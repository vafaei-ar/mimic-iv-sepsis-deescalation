#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd


def _latest_final(run_dir: Path) -> Path:
    candidates = sorted((run_dir / "inference_reruns").glob("final_vital_corrected_final_*"))
    if not candidates:
        raise FileNotFoundError("No final_vital_corrected_final_* inference rerun found")
    return candidates[-1]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _candidate_bootstrap_inventory(final_dir: Path) -> list[dict]:
    rows: list[dict] = []
    for path in sorted(final_dir.rglob("*.csv")):
        low = path.name.lower()
        if "bootstrap" not in low and "progressive" not in low and "primary" not in low:
            continue
        item = {
            "path": str(path.relative_to(final_dir)),
            "size_bytes": int(path.stat().st_size),
            "sha256": _sha256(path),
        }
        try:
            df = pd.read_csv(path)
            item["rows"] = int(len(df))
            item["columns"] = [str(c) for c in df.columns]
            for col in ["bootstrap_success", "n_success", "n_requested", "replicate", "rep", "model", "analysis", "outcome"]:
                if col in df.columns:
                    vals = df[col].dropna()
                    item[f"{col}_unique_count"] = int(vals.nunique())
                    if len(vals) and vals.nunique() <= 12:
                        item[f"{col}_values"] = [str(v) for v in vals.unique().tolist()]
        except Exception as exc:
            item["read_error"] = type(exc).__name__
        rows.append(item)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a safe aggregate publication-freeze bundle from the final corrected MIMIC inference run.")
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/publication_freeze/latest"))
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    final_dir = _latest_final(run_dir)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    primary_path = final_dir / "tables" / "primary_secondary_outcomes.csv"
    progressive_path = final_dir / "tables" / "progressive_adjustment.csv"
    weighting_dir = final_dir / "final_weighting"
    weighting_point_path = weighting_dir / "final_weighting_point_estimates.csv"
    weighting_ci_path = weighting_dir / "final_weighting_bootstrap_ci.csv"
    weighting_diag_path = weighting_dir / "final_weighting_bootstrap_diagnostics.csv"

    required = [primary_path, progressive_path, weighting_point_path, weighting_ci_path, weighting_diag_path]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required aggregate outputs: " + ", ".join(missing))

    primary = pd.read_csv(primary_path)
    progressive = pd.read_csv(progressive_path)
    # The three weighting files are required publication artifacts and are copied below.
    # We only need to verify that they are readable here; their contents are not used to
    # derive the primary-vs-M4 CI audit.
    pd.read_csv(weighting_point_path)
    pd.read_csv(weighting_ci_path)
    pd.read_csv(weighting_diag_path)

    mort = primary.loc[primary["analysis"] == "30-day mortality"].iloc[0]
    afd = primary.loc[primary["analysis"] == "antibiotic-free days"].iloc[0]
    m4 = progressive.loc[progressive["model"].str.startswith("M4")].iloc[0]

    primary_rd = float(mort["risk_difference"])
    m4_rd = float(m4["risk_difference"])
    primary_ci = [float(mort["lower_95"]), float(mort["upper_95"])]
    m4_ci = [float(m4["rd_lower_95"]), float(m4["rd_upper_95"])]
    point_delta = primary_rd - m4_rd
    lower_delta = primary_ci[0] - m4_ci[0]
    upper_delta = primary_ci[1] - m4_ci[1]

    m4_success = None
    for col in ["bootstrap_success", "n_bootstrap_success", "bootstrap_n", "n_success"]:
        if col in progressive.columns and pd.notna(m4[col]):
            m4_success = int(m4[col])
            break

    ci_audit = {
        "purpose": "Audit only; no model, cohort, exposure, outcome, covariate, or weighting changes.",
        "primary_source_file": str(primary_path.relative_to(final_dir)),
        "m4_source_file": str(progressive_path.relative_to(final_dir)),
        "primary_rd": primary_rd,
        "m4_rd": m4_rd,
        "point_estimate_delta": point_delta,
        "same_point_estimate_within_1e-12": bool(abs(point_delta) <= 1e-12),
        "primary_ci95": primary_ci,
        "m4_ci95": m4_ci,
        "ci_endpoint_deltas_primary_minus_m4": [lower_delta, upper_delta],
        "same_ci_within_1e-12": bool(abs(lower_delta) <= 1e-12 and abs(upper_delta) <= 1e-12),
        "primary_bootstrap_success": int(mort["bootstrap_success"]),
        "m4_bootstrap_success": m4_success,
        "primary_file_sha256": _sha256(primary_path),
        "m4_file_sha256": _sha256(progressive_path),
        "bootstrap_candidate_inventory": _candidate_bootstrap_inventory(final_dir),
        "interpretation_rule": "If point estimates are identical but CIs differ, treat the difference as bootstrap-procedure/provenance variation until the generating paths are shown to be identical. Do not choose a CI based on favorability.",
    }

    summary = {
        "source_run": str(final_dir.relative_to(Path.cwd())) if final_dir.is_relative_to(Path.cwd()) else str(final_dir),
        "primary_30d_mortality": {
            "risk_deescalated": float(mort["risk_deescalated_stopped"]),
            "risk_continued": float(mort["risk_continued"]),
            "risk_difference": primary_rd,
            "risk_difference_ci95": primary_ci,
            "risk_ratio": float(mort["risk_ratio"]),
            "bootstrap_success": int(mort["bootstrap_success"]),
        },
        "antibiotic_free_days": {
            "mean_difference": float(afd["mean_difference"]),
            "ci95": [float(afd["lower_95"]), float(afd["upper_95"])],
        },
        "progressive_m4": {
            "risk_difference": m4_rd,
            "risk_difference_ci95": m4_ci,
            "risk_ratio": float(m4["risk_ratio"]),
            "max_post_smd": float(m4["max_post_smd"]),
            "balance_warning_gt_0_10": bool(float(m4["max_post_smd"]) > 0.10),
        },
        "ci_audit": ci_audit,
        "interpretation_note": "Fully adjusted mortality is compatible with no clear benefit or harm; antibiotic burden is lower after de-escalation. Residual IPTW imbalance remains above 0.10 and should be reported transparently.",
        "data_safety": "Aggregate derived outputs only; no raw or row-level MIMIC-IV data included.",
    }

    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    shutil.copy2(primary_path, out / "primary_secondary_outcomes.csv")
    shutil.copy2(progressive_path, out / "progressive_adjustment.csv")
    shutil.copy2(weighting_point_path, out / "final_weighting_point_estimates.csv")
    shutil.copy2(weighting_ci_path, out / "final_weighting_bootstrap_ci.csv")
    shutil.copy2(weighting_diag_path, out / "final_weighting_bootstrap_diagnostics.csv")
    print(json.dumps(summary, indent=2))
    print(f"Safe publication-freeze outputs: {out}")


if __name__ == "__main__":
    main()

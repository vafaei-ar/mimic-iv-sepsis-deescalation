#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd


def _latest_final(run_dir: Path) -> Path:
    candidates = sorted((run_dir / "inference_reruns").glob("final_vital_corrected_final_*"))
    if not candidates:
        raise FileNotFoundError("No final_vital_corrected_final_* inference rerun found")
    return candidates[-1]


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
    weighting_point = pd.read_csv(weighting_point_path)
    weighting_ci = pd.read_csv(weighting_ci_path)
    weighting_diag = pd.read_csv(weighting_diag_path)

    mort = primary.loc[primary["analysis"] == "30-day mortality"].iloc[0]
    afd = primary.loc[primary["analysis"] == "antibiotic-free days"].iloc[0]
    m4 = progressive.loc[progressive["model"].str.startswith("M4")].iloc[0]

    summary = {
        "source_run": str(final_dir.relative_to(Path.cwd())) if final_dir.is_relative_to(Path.cwd()) else str(final_dir),
        "primary_30d_mortality": {
            "risk_deescalated": float(mort["risk_deescalated_stopped"]),
            "risk_continued": float(mort["risk_continued"]),
            "risk_difference": float(mort["risk_difference"]),
            "risk_difference_ci95": [float(mort["lower_95"]), float(mort["upper_95"])],
            "risk_ratio": float(mort["risk_ratio"]),
            "bootstrap_success": int(mort["bootstrap_success"]),
        },
        "antibiotic_free_days": {
            "mean_difference": float(afd["mean_difference"]),
            "ci95": [float(afd["lower_95"]), float(afd["upper_95"])],
        },
        "progressive_m4": {
            "risk_difference": float(m4["risk_difference"]),
            "risk_difference_ci95": [float(m4["rd_lower_95"]), float(m4["rd_upper_95"])],
            "risk_ratio": float(m4["risk_ratio"]),
            "max_post_smd": float(m4["max_post_smd"]),
            "balance_warning_gt_0_10": bool(float(m4["max_post_smd"]) > 0.10),
        },
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

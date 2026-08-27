#!/usr/bin/env python3
"""Collect frozen aggregate publication source tables into one sanitized JSON bundle.

This task reads only project-local aggregate outputs already produced by the frozen
MIMIC and PSU pipelines. It does not read patient-level data and does not recompute
any estimand.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "publication_integration" / "source_snapshot"

MIMIC_BASE = ROOT / "outputs" / "mimic" / "mimic_iv_v5_7_final_20260820T003506Z" / "inference_reruns" / "final_vital_corrected_final_20260825T010041Z"
PSU_BASE = ROOT / "outputs" / "psu_publication_parity" / "latest"

SOURCES = {
    "mimic_primary_secondary": MIMIC_BASE / "tables" / "primary_secondary_outcomes.csv",
    "mimic_progressive_adjustment": MIMIC_BASE / "tables" / "progressive_adjustment.csv",
    "mimic_final_weighting_point_estimates": MIMIC_BASE / "final_weighting" / "final_weighting_point_estimates.csv",
    "mimic_final_weighting_bootstrap_ci": MIMIC_BASE / "final_weighting" / "final_weighting_bootstrap_ci.csv",
    "psu_parity_report": PSU_BASE / "parity_report.json",
    "psu_point_estimates": PSU_BASE / "point_estimates" / "point_estimates.csv",
    "psu_bootstrap_ci": PSU_BASE / "bootstrap" / "bootstrap_ci.csv",
    "psu_ps_balance_summary": PSU_BASE / "ps_balance" / "summary.json",
    "psu_robustness_bootstrap_summary": PSU_BASE / "robustness_bootstrap" / "robustness_bootstrap_summary.csv",
    "psu_outcome_overall_summary": PSU_BASE / "outcome_freeze" / "outcome_overall_summary.csv",
}


def read_aggregate(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path)
        return {
            "columns": list(df.columns),
            "n_rows": int(len(df)),
            "records": df.where(pd.notna(df), None).to_dict(orient="records"),
        }
    if path.suffix.lower() == ".json":
        return json.loads(path.read_text())
    raise ValueError(f"Unsupported aggregate source: {path}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    bundle = {
        "purpose": "Frozen aggregate publication source snapshot; no new analysis.",
        "mimic_source_root": str(MIMIC_BASE.relative_to(ROOT)),
        "psu_source_root": str(PSU_BASE.relative_to(ROOT)),
        "sources": {},
        "data_safety": "Aggregate outputs only; no patient-level MIMIC or PSU data included.",
    }
    rows = []
    for name, path in SOURCES.items():
        rel = str(path.relative_to(ROOT))
        content = read_aggregate(path)
        bundle["sources"][name] = {"path": rel, "content": content}
        rows.append({"source": name, "path": rel, "suffix": path.suffix.lower(), "exists": True})

    (OUT / "publication_source_snapshot.json").write_text(json.dumps(bundle, indent=2, default=str))
    pd.DataFrame(rows).to_csv(OUT / "publication_source_files.csv", index=False)
    print(json.dumps({"n_sources": len(rows), "output_dir": str(OUT), "status": "ok"}, indent=2))


if __name__ == "__main__":
    main()

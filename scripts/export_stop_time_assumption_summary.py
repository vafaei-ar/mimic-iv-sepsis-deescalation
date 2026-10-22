#!/usr/bin/env python3
"""Export only existing, aggregate MIMIC stop-time sensitivity results.

This report audits previously generated deterministic point estimates.
No cohort reconstruction, model refitting, or restricted row-level output.
"""
from __future__ import annotations

import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/publication_integration/stop_time_audit"
FILENAME = "stop_time_assumption_sensitivity.csv"
NEEDED = {"stop_fill_hours", "n", "risk_treated", "risk_control", "risk_difference", "risk_ratio"}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    candidates = sorted((ROOT / "outputs/mimic").glob("**/" + FILENAME))
    sources = []
    for path in candidates:
        if not path.is_file():
            continue
        df = pd.read_csv(path)
        if not NEEDED.issubset(df.columns) or len(df) > 12:
            raise RuntimeError("Unexpected stop-time sensitivity summary format")
        if not set(df["stop_fill_hours"].dropna().astype(int)).issubset({12, 24, 48}):
            raise RuntimeError("Unexpected stop-fill hour values")
        for row in df.to_dict(orient="records"):
            sources.append({"source_path": str(path.relative_to(ROOT)),
                "stop_fill_hours": int(row["stop_fill_hours"]), "n": int(row["n"]),
                "risk_treated": float(row["risk_treated"]),
                "risk_control": float(row["risk_control"]),
                "risk_difference": float(row["risk_difference"]),
                "risk_ratio": float(row["risk_ratio"])})
    summary = {"found": bool(sources), "n_source_files": len(candidates),
               "point_estimates_only": True, "source_rows": sources}
    (OUT / "stop_time_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps({"found": summary["found"], "n_source_files": len(candidates),
                      "n_rows": len(sources)}, indent=2))


if __name__ == "__main__":
    main()

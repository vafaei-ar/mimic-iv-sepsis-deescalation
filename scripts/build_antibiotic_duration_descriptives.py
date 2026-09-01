#!/usr/bin/env python3
"""Build aggregate antibiotic-duration descriptors for clinical manuscript reporting.

This task does not refit any model or alter any frozen estimand. For MIMIC-IV it reads
the final vital-corrected analytic cohort and summarizes the already-constructed
``antibiotic_days`` outcome, which counts calendar days with observed systemic
antibiotic coverage during the 30-day period beginning at the 96-hour landmark.
For Penn State it reads the frozen aggregate outcome summary produced by the
publication-parity pipeline, where ``antibiotic_days_30d`` is the corresponding
calendar-day PRESCRIBING proxy after the 96-hour landmark.

Only aggregate descriptors are exported. No row-level clinical data leave the local
machine.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MIMIC = ROOT / "outputs/mimic/mimic_iv_v5_7_final_20260820T003506Z/audits/vital_repair/analysis_cohort_vital_corrected.csv"
PSU = ROOT / "outputs/psu_publication_parity/latest/outcome_freeze/outcome_overall_summary.csv"
OUT = ROOT / "outputs/publication_integration/clinical_duration"
REVIEWER_SUMMARY = ROOT / "outputs/publication_integration/reviewer_support/reviewer_support_final_summary.md"

EXPECTED_MIMIC_N = 9589
EXPECTED_MIMIC_DEESC = 1863
EXPECTED_MIMIC_CONT = 7726
EXPECTED_PSU_N = 19841


def summarize(series: pd.Series) -> dict[str, float | int]:
    x = pd.to_numeric(series, errors="coerce").dropna()
    return {
        "n": int(len(x)),
        "mean_days": float(x.mean()),
        "sd_days": float(x.std(ddof=1)),
        "p25_days": float(x.quantile(0.25)),
        "median_days": float(x.quantile(0.50)),
        "p75_days": float(x.quantile(0.75)),
        "min_days": float(x.min()),
        "max_days": float(x.max()),
    }


def main() -> None:
    if not MIMIC.exists():
        raise FileNotFoundError(MIMIC)
    if not PSU.exists():
        raise FileNotFoundError(PSU)
    OUT.mkdir(parents=True, exist_ok=True)

    d = pd.read_csv(MIMIC, low_memory=False)
    if len(d) != EXPECTED_MIMIC_N:
        raise RuntimeError(f"MIMIC cohort parity failure: {len(d)} != {EXPECTED_MIMIC_N}")
    if int((d["A"] == 1).sum()) != EXPECTED_MIMIC_DEESC or int((d["A"] == 0).sum()) != EXPECTED_MIMIC_CONT:
        raise RuntimeError("MIMIC treatment-count parity failure")
    if "antibiotic_days" not in d.columns:
        raise RuntimeError("Final MIMIC cohort does not contain antibiotic_days")

    rows: list[dict] = []
    for group, mask in [
        ("Overall", pd.Series(True, index=d.index)),
        ("De-escalation or stopping", d["A"] == 1),
        ("Continued broad-spectrum", d["A"] == 0),
    ]:
        row = {
            "dataset": "MIMIC-IV",
            "group": group,
            "measure": "Observed systemic antibiotic days during 30-day post-landmark period",
            "measurement_note": "Calendar days with observed systemic antibiotic coverage from the 96-hour landmark; not necessarily one uninterrupted course.",
        }
        row.update(summarize(d.loc[mask, "antibiotic_days"]))
        rows.append(row)

    psu = pd.read_csv(PSU)
    hit = psu.loc[psu["outcome"] == "antibiotic_days_30d"]
    if len(hit) != 1:
        raise RuntimeError("Expected one PSU antibiotic_days_30d row in frozen aggregate outcome summary")
    p = hit.iloc[0]
    if int(p["n"]) != EXPECTED_PSU_N:
        raise RuntimeError(f"PSU cohort parity failure: {int(p['n'])} != {EXPECTED_PSU_N}")
    rows.append({
        "dataset": "Penn State",
        "group": "Overall",
        "measure": "Observed systemic antibiotic days during 30-day post-landmark period",
        "measurement_note": "Calendar days with the frozen systemic-antibiotic PRESCRIBING proxy from the 96-hour landmark; ordered treatment, not verified administration.",
        "n": int(p["n"]),
        "mean_days": float(p["mean"]),
        "sd_days": float(p["sd"]),
        "p25_days": float(p["p25"]),
        "median_days": float(p["median"]),
        "p75_days": float(p["p75"]),
        "min_days": float(p["min"]),
        "max_days": float(p["max"]),
    })

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "antibiotic_duration_descriptives.csv", index=False)
    metadata = {
        "purpose": "Clinical descriptive response to coauthor request for mean/median antibiotic duration.",
        "primary_science_changed": False,
        "model_refit": False,
        "bootstrap_rerun": False,
        "mimic_source": str(MIMIC.relative_to(ROOT)),
        "psu_source": str(PSU.relative_to(ROOT)),
        "interpretation_guardrail": "These are observed antibiotic-exposure days after the 96-hour landmark, not necessarily uninterrupted total course duration.",
    }
    (OUT / "antibiotic_duration_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    if REVIEWER_SUMMARY.exists():
        display = out[["dataset", "group", "n", "mean_days", "sd_days", "p25_days", "median_days", "p75_days"]].copy()
        lines = [
            "\n## Clinically interpretable antibiotic-duration descriptors\n",
            "These values count calendar days with observed systemic antibiotic exposure during the 30-day period after the 96-hour landmark; they are not necessarily uninterrupted total course durations.\n",
            "| Dataset | Group | n | Mean days | SD | P25 | Median | P75 |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
        for r in display.itertuples(index=False):
            lines.append(f"| {r.dataset} | {r.group} | {int(r.n):,} | {r.mean_days:.2f} | {r.sd_days:.2f} | {r.p25_days:.1f} | {r.median_days:.1f} | {r.p75_days:.1f} |")
        with REVIEWER_SUMMARY.open("a", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    print(out.to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build MIMIC-IV and Penn State baseline-characteristics publication artifacts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

import build_mimic_baseline_characteristics as mimic
import build_psu_baseline_characteristics as psu

OUT = Path("outputs/publication_integration/baseline_characteristics")


def combine_markdown(out: Path) -> Path:
    mimic_md = (out / "mimic_baseline_characteristics.md").read_text(encoding="utf-8")
    psu_md = (out / "psu_baseline_characteristics.md").read_text(encoding="utf-8")
    combined = [
        "# ESM baseline characteristics",
        "",
        "## Panel A. MIMIC-IV primary cohort",
        "",
        mimic_md.replace("# Candidate baseline characteristics table\n\n", "", 1).strip(),
        "",
        "## Panel B. Penn State modified external-replication cohort",
        "",
        psu_md.replace("# Penn State baseline characteristics\n\n", "", 1).strip(),
        "",
        "Cross-dataset note: the panels intentionally report the variables available and frozen within each data source rather than implying exact covariate harmonization. Penn State lacks several MIMIC-specific ICU-context, urine-output, and microbiology-intensity measures.",
    ]
    path = out / "baseline_characteristics_combined.md"
    path.write_text("\n".join(combined) + "\n", encoding="utf-8")
    return path


def _cohort_provenance(mimic_df: pd.DataFrame, psu_df: pd.DataFrame) -> dict:
    """Return aggregate-only reporting diagnostics for reviewer-facing methods text."""
    out: dict[str, object] = {}

    if "subject_id" in mimic_df.columns:
        counts = mimic_df.groupby("subject_id", dropna=True).size()
        out.update(
            {
                "mimic_unique_subjects": int(counts.size),
                "mimic_subjects_with_multiple_eligible_admissions": int((counts > 1).sum()),
                "mimic_admissions_from_subjects_with_multiple_eligible_admissions": int(
                    counts.loc[counts > 1].sum()
                ),
                "mimic_max_eligible_admissions_per_subject": int(counts.max()) if len(counts) else 0,
            }
        )
    if "admittime" in mimic_df.columns:
        dt = pd.to_datetime(mimic_df["admittime"], errors="coerce").dropna()
        if len(dt):
            out["mimic_admission_date_min"] = dt.min().date().isoformat()
            out["mimic_admission_date_max"] = dt.max().date().isoformat()

    if "patid" in psu_df.columns:
        counts = psu_df.groupby("patid", dropna=True).size()
        out.update(
            {
                "psu_unique_patients": int(counts.size),
                "psu_patients_with_multiple_eligible_encounters": int((counts > 1).sum()),
                "psu_encounters_from_patients_with_multiple_eligible_encounters": int(
                    counts.loc[counts > 1].sum()
                ),
                "psu_max_eligible_encounters_per_patient": int(counts.max()) if len(counts) else 0,
            }
        )
    if "admit_date" in psu_df.columns:
        dt = pd.to_datetime(psu_df["admit_date"], errors="coerce").dropna()
        if len(dt):
            out["psu_admission_date_min"] = dt.min().date().isoformat()
            out["psu_admission_date_max"] = dt.max().date().isoformat()

    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("psu_data_root", type=Path)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    mimic.main()
    mimic_df = pd.read_csv(mimic.COHORT, low_memory=False)

    df, balance, _, fitmeta = psu.reconstruct(args.psu_data_root)
    detailed, formatted, metadata = psu.build_table(df, balance, fitmeta)
    detailed.to_csv(OUT / "psu_baseline_characteristics_detailed.csv", index=False)
    formatted.to_csv(OUT / "psu_baseline_characteristics_formatted.csv", index=False)
    psu.write_markdown(formatted, metadata, OUT / "psu_baseline_characteristics.md")
    (OUT / "psu_baseline_characteristics_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    combined = combine_markdown(OUT)
    summary = {
        "mimic_cohort_n": mimic.EXPECTED_N,
        "psu_cohort_n": psu.EXPECTED_N,
        **_cohort_provenance(mimic_df, df),
        "combined_markdown": str(combined),
        "privacy": "Aggregate/sanitized outputs only; no row-level MIMIC or PSU data exported.",
    }
    (OUT / "baseline_characteristics_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

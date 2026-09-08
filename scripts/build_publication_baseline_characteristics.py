#!/usr/bin/env python3
"""Build MIMIC-IV and Penn State baseline-characteristics publication artifacts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
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


def _observed_continuous_summary(x: pd.Series) -> dict[str, float | int]:
    z = pd.to_numeric(x, errors="coerce")
    observed = z.dropna()
    return {
        "n": int(len(z)),
        "nonmissing_n": int(observed.size),
        "missing_n": int(z.isna().sum()),
        "missing_percent": float(100 * z.isna().mean()),
        "mean": float(observed.mean()) if observed.size else np.nan,
        "sd": float(observed.std(ddof=1)) if observed.size > 1 else np.nan,
        "median": float(observed.median()) if observed.size else np.nan,
        "q1": float(observed.quantile(0.25)) if observed.size else np.nan,
        "q3": float(observed.quantile(0.75)) if observed.size else np.nan,
    }


def _format_continuous(summary: dict[str, float | int], digits: int) -> str:
    if not np.isfinite(float(summary["median"])):
        return "NA"
    fmt = f"{{:.{digits}f}}"
    return (
        f"{fmt.format(float(summary['median']))} "
        f"[{fmt.format(float(summary['q1']))}, {fmt.format(float(summary['q3']))}]"
    )


def _restore_observed_mimic_descriptives(mimic_df: pd.DataFrame) -> None:
    """Remove median-imputed values from descriptive summaries without changing SMDs.

    The frozen analytic cohort stores several continuous trajectory variables after
    median imputation, together with explicit ``*_missing`` indicators. The propensity
    model and SMD diagnostics must continue to use that frozen analysis representation,
    but a clinician-facing baseline table should summarize the measurements that were
    actually observed. This function updates only the descriptive cells and missingness
    counts in the already-built aggregate table. No weights, balance values, outcomes,
    or treatment effects are recomputed.
    """
    detailed_path = OUT / "mimic_baseline_characteristics_detailed.csv"
    formatted_path = OUT / "mimic_baseline_characteristics_formatted.csv"
    metadata_path = OUT / "mimic_baseline_characteristics_metadata.json"

    detailed = pd.read_csv(detailed_path)
    formatted = pd.read_csv(formatted_path)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    group1_col = f"De-escalation or stopping (n={mimic.EXPECTED_DEESC:,})"
    group0_col = f"Continued broad-spectrum (n={mimic.EXPECTED_CONT:,})"

    corrected: list[str] = []
    for _, var, label, kind, digits in mimic.CHARACTERISTICS:
        if kind != "continuous":
            continue
        missing_col = f"{var}_missing"
        if missing_col not in mimic_df.columns:
            continue

        summaries = {}
        for a, prefix in [(1, "deescalated_or_stopped"), (0, "continued_broad")]:
            sub = mimic_df.loc[mimic_df["A"] == a, [var, missing_col]].copy()
            values = pd.to_numeric(sub[var], errors="coerce").mask(
                pd.to_numeric(sub[missing_col], errors="coerce").fillna(0).astype(int) == 1
            )
            summaries[prefix] = _observed_continuous_summary(values)

        s1 = summaries["deescalated_or_stopped"]
        s0 = summaries["continued_broad"]
        display1 = _format_continuous(s1, digits)
        display0 = _format_continuous(s0, digits)

        mask = detailed["variable"].eq(var)
        if int(mask.sum()) != 1:
            raise RuntimeError(f"Expected one detailed MIMIC baseline row for {var}")
        idx = detailed.index[mask][0]
        detailed.loc[idx, "deescalated_or_stopped_display"] = display1
        detailed.loc[idx, "continued_broad_display"] = display0
        for prefix, summary in summaries.items():
            detailed.loc[idx, f"{prefix}_nonmissing_n"] = summary["nonmissing_n"]
            detailed.loc[idx, f"{prefix}_missing_n"] = summary["missing_n"]
            detailed.loc[idx, f"{prefix}_missing_percent"] = summary["missing_percent"]
            for stat in ["mean", "sd", "median", "q1", "q3"]:
                detailed.loc[idx, f"{prefix}_{stat}"] = summary[stat]

        fmask = formatted["characteristic"].eq(label)
        if int(fmask.sum()) != 1:
            raise RuntimeError(f"Expected one formatted MIMIC baseline row for {label}")
        fidx = formatted.index[fmask][0]
        formatted.loc[fidx, group1_col] = display1
        formatted.loc[fidx, group0_col] = display0
        corrected.append(var)

    detailed.to_csv(detailed_path, index=False)
    formatted.to_csv(formatted_path, index=False)
    metadata["summary_convention"] = (
        "Continuous variables: median [IQR]. For continuous variables with stored original "
        "missingness indicators, descriptive summaries exclude values inserted by the frozen "
        "median-imputation step; derived analysis scores retain their frozen construction. "
        "Categorical variables: n (%) with missing binary values coded as zero to match the "
        "frozen PS analysis convention."
    )
    metadata["descriptive_imputation_correction"] = {
        "applied": True,
        "variables": corrected,
        "scientific_effect": "Descriptive summaries only; propensity weights, SMDs, outcomes, and treatment effects unchanged.",
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    mimic.write_markdown(formatted, metadata, OUT / "mimic_baseline_characteristics.md")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("psu_data_root", type=Path)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    mimic.main()
    mimic_df = pd.read_csv(mimic.COHORT, low_memory=False)
    _restore_observed_mimic_descriptives(mimic_df)

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

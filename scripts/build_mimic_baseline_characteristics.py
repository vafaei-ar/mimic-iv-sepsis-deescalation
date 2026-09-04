#!/usr/bin/env python3
"""Build aggregate MIMIC-IV baseline characteristics for manuscript reporting.

The table is descriptive only. It reads the frozen vital-corrected primary MIMIC
analytic cohort, reproduces the frozen stabilized-IPTW weights for balance
reporting, and exports aggregate summaries only. No outcome or treatment-effect
estimate is computed.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from sepsis_deescalation.specification import CANDIDATE_PS_VARS
from sepsis_deescalation.stats import balance_table, fit_stabilized_iptw

BASE_RUN = Path("outputs/mimic/mimic_iv_v5_7_final_20260820T003506Z")
COHORT = BASE_RUN / "audits/vital_repair/analysis_cohort_vital_corrected.csv"
OUT = Path("outputs/publication_integration/baseline_characteristics")

EXPECTED_N = 9589
EXPECTED_DEESC = 1863
EXPECTED_CONT = 7726
EXPECTED_MAX_POST_SMD = 0.13261
EXPECTED_WORST_POST_SMD = "temperature_48_72h"

# Curated clinician-facing subset of pre-decision variables. This is not a new
# propensity specification; the complete frozen PS specification remains in
# sepsis_deescalation.specification.CANDIDATE_PS_VARS.
CHARACTERISTICS = [
    # Demographics
    ("Demographics", "age", "Age, years", "continuous", 1),
    ("Demographics", "sex_male", "Male sex", "binary", 1),
    ("Demographics", "race_white", "White race", "binary", 1),
    # Comorbidity and ICU context
    ("Comorbidity and ICU context", "comorb", "Any recorded comorbidity", "binary", 1),
    ("Comorbidity and ICU context", "heart_failure", "Heart failure", "binary", 1),
    ("Comorbidity and ICU context", "chronic_kidney", "Chronic kidney disease", "binary", 1),
    (
        "Comorbidity and ICU context",
        "hours_admit_to_icu",
        "Admission-to-ICU interval, h",
        "continuous",
        1,
    ),
    ("Comorbidity and ICU context", "micu", "Medical ICU", "binary", 1),
    ("Comorbidity and ICU context", "sicu", "Surgical ICU", "binary", 1),
    ("Comorbidity and ICU context", "cardiac_icu", "Cardiac ICU", "binary", 1),
    ("Comorbidity and ICU context", "neuro_icu", "Neurologic ICU", "binary", 1),
    # Severity and organ support before the treatment decision
    (
        "Severity and organ support",
        "vent_proc",
        "Mechanical ventilation, pre-72 h",
        "binary",
        1,
    ),
    (
        "Severity and organ support",
        "vasopressor_any_0_24h",
        "Vasopressor exposure, 0-24 h",
        "binary",
        1,
    ),
    (
        "Severity and organ support",
        "vasopressor_any_48_72h",
        "Vasopressor exposure, 48-72 h",
        "binary",
        1,
    ),
    (
        "Severity and organ support",
        "lactate_early_last_0_24h",
        "Last lactate, 0-24 h, mmol/L",
        "continuous",
        1,
    ),
    (
        "Severity and organ support",
        "lactate_late_last_48_72h",
        "Last lactate, 48-72 h, mmol/L",
        "continuous",
        1,
    ),
    (
        "Severity and organ support",
        "creatinine_late_last_48_72h",
        "Last creatinine, 48-72 h, mg/dL",
        "continuous",
        2,
    ),
    (
        "Severity and organ support",
        "wbc_late_last_48_72h",
        "Last WBC count, 48-72 h",
        "continuous",
        1,
    ),
    (
        "Severity and organ support",
        "platelet_late_worst_48_72h",
        "Worst platelet count, 48-72 h",
        "continuous",
        0,
    ),
    (
        "Severity and organ support",
        "bilirubin_late_worst_48_72h",
        "Worst bilirubin, 48-72 h, mg/dL",
        "continuous",
        1,
    ),
    (
        "Severity and organ support",
        "heart_rate_48_72h",
        "Maximum heart rate, 48-72 h, beats/min",
        "continuous",
        0,
    ),
    (
        "Severity and organ support",
        "temperature_48_72h",
        "Temperature, 48-72 h, degC",
        "continuous",
        1,
    ),
    (
        "Severity and organ support",
        "urine_output_ml_48_72h",
        "Urine output, 48-72 h, mL",
        "continuous",
        0,
    ),
    (
        "Severity and organ support",
        "sofa_like_0_24h",
        "SOFA-like score, 0-24 h",
        "continuous",
        1,
    ),
    (
        "Severity and organ support",
        "sofa_like_48_72h",
        "SOFA-like score, 48-72 h",
        "continuous",
        1,
    ),
    # Treatment intensity before the decision
    (
        "Antibiotic treatment intensity",
        "systemic_abx_agents_pre72",
        "Systemic antibiotic agents, pre-72 h",
        "continuous",
        0,
    ),
    (
        "Antibiotic treatment intensity",
        "broad_abx_hours_pre72",
        "Broad-spectrum antibiotic hours, pre-72 h",
        "continuous",
        1,
    ),
    (
        "Antibiotic treatment intensity",
        "broad_abx_agents_pre72",
        "Broad-spectrum agents, pre-72 h",
        "continuous",
        0,
    ),
    (
        "Antibiotic treatment intensity",
        "anti_mrsa_pre72",
        "Anti-MRSA therapy, pre-72 h",
        "binary",
        1,
    ),
    (
        "Antibiotic treatment intensity",
        "antipseudomonal_pre72",
        "Antipseudomonal therapy, pre-72 h",
        "binary",
        1,
    ),
    (
        "Antibiotic treatment intensity",
        "carbapenem_pre72",
        "Carbapenem use, pre-72 h",
        "binary",
        1,
    ),
    # Microbiology and diagnostic sampling. These are not validated infection-site labels.
    (
        "Microbiology and diagnostic sampling",
        "distinct_specimen_types_pre72",
        "Distinct microbiology specimen types, pre-72 h",
        "continuous",
        0,
    ),
    (
        "Microbiology and diagnostic sampling",
        "blood_culture_pre72",
        "Blood culture sampled, pre-72 h",
        "binary",
        1,
    ),
    (
        "Microbiology and diagnostic sampling",
        "respiratory_culture_pre72",
        "Respiratory culture sampled, pre-72 h",
        "binary",
        1,
    ),
    (
        "Microbiology and diagnostic sampling",
        "urine_culture_pre72",
        "Urine culture sampled, pre-72 h",
        "binary",
        1,
    ),
    (
        "Microbiology and diagnostic sampling",
        "sterile_fluid_culture_pre72",
        "Sterile-fluid culture sampled, pre-72 h",
        "binary",
        1,
    ),
    (
        "Microbiology and diagnostic sampling",
        "repeat_micro_48_72h",
        "Repeat microbiology sampling, 48-72 h",
        "binary",
        1,
    ),
]


def _continuous_summary(x: pd.Series) -> dict[str, float | int]:
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


def _binary_summary(x: pd.Series) -> dict[str, float | int]:
    raw = pd.to_numeric(x, errors="coerce")
    # The frozen MIMIC PS implementation codes missing binary covariates as zero.
    z = raw.fillna(0).clip(0, 1)
    positive_n = int(z.sum())
    return {
        "n": int(len(z)),
        "nonmissing_n": int(raw.notna().sum()),
        "missing_n": int(raw.isna().sum()),
        "missing_percent": float(100 * raw.isna().mean()),
        "positive_n": positive_n,
        "positive_percent": float(100 * positive_n / len(z)) if len(z) else np.nan,
    }


def _format_continuous(summary: dict[str, float | int], digits: int) -> str:
    if not np.isfinite(float(summary["median"])):
        return "NA"
    fmt = f"{{:.{digits}f}}"
    return (
        f"{fmt.format(float(summary['median']))} "
        f"[{fmt.format(float(summary['q1']))}, {fmt.format(float(summary['q3']))}]"
    )


def _format_binary(summary: dict[str, float | int]) -> str:
    return f"{int(summary['positive_n']):,} ({float(summary['positive_percent']):.1f}%)"


def build_table(d: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    missing = [var for _, var, _, _, _ in CHARACTERISTICS if var not in d.columns]
    if missing:
        raise RuntimeError(f"Frozen MIMIC cohort is missing expected descriptive variables: {missing}")
    if len(d) != EXPECTED_N:
        raise RuntimeError(f"Frozen MIMIC cohort n={len(d)}, expected {EXPECTED_N}")
    n1 = int((d["A"] == 1).sum())
    n0 = int((d["A"] == 0).sum())
    if (n1, n0) != (EXPECTED_DEESC, EXPECTED_CONT):
        raise RuntimeError(
            "Frozen treatment counts do not match publication lock: "
            f"observed A=1 {n1}, A=0 {n0}"
        )

    # Reproduce the frozen weights only for post-weighting balance diagnostics.
    # No outcome is read or analyzed in this script.
    w, _, _ = fit_stabilized_iptw(d, CANDIDATE_PS_VARS)
    all_balance = balance_table(w, CANDIDATE_PS_VARS)
    max_row = all_balance.sort_values("after", ascending=False).iloc[0]
    max_post = float(max_row["after"])
    worst_post = str(max_row["variable"])
    if abs(max_post - EXPECTED_MAX_POST_SMD) > 5e-4:
        raise RuntimeError(
            f"Primary balance parity failed: max post-SMD {max_post:.6f}, "
            f"expected about {EXPECTED_MAX_POST_SMD:.5f}"
        )
    if worst_post != EXPECTED_WORST_POST_SMD:
        raise RuntimeError(
            f"Primary balance parity failed: worst variable {worst_post}, "
            f"expected {EXPECTED_WORST_POST_SMD}"
        )

    balance = all_balance.set_index("variable")
    detailed_rows: list[dict] = []
    formatted_rows: list[dict] = []

    for section, var, label, kind, digits in CHARACTERISTICS:
        g1 = d.loc[d["A"] == 1, var]
        g0 = d.loc[d["A"] == 0, var]
        if kind == "binary":
            s1 = _binary_summary(g1)
            s0 = _binary_summary(g0)
            display1 = _format_binary(s1)
            display0 = _format_binary(s0)
        else:
            s1 = _continuous_summary(g1)
            s0 = _continuous_summary(g0)
            display1 = _format_continuous(s1, digits)
            display0 = _format_continuous(s0, digits)

        if var in balance.index:
            smd_before = float(balance.loc[var, "before"])
            smd_after = float(balance.loc[var, "after"])
        else:
            smd_before = np.nan
            smd_after = np.nan

        detailed_rows.append(
            {
                "section": section,
                "variable": var,
                "characteristic": label,
                "type": kind,
                "deescalated_or_stopped_display": display1,
                "continued_broad_display": display0,
                "deescalated_or_stopped_n": s1["n"],
                "continued_broad_n": s0["n"],
                "deescalated_or_stopped_nonmissing_n": s1["nonmissing_n"],
                "continued_broad_nonmissing_n": s0["nonmissing_n"],
                "deescalated_or_stopped_missing_n": s1["missing_n"],
                "continued_broad_missing_n": s0["missing_n"],
                "deescalated_or_stopped_missing_percent": s1["missing_percent"],
                "continued_broad_missing_percent": s0["missing_percent"],
                "absolute_smd_before_weighting": smd_before,
                "absolute_smd_after_weighting": smd_after,
                **{
                    f"deescalated_or_stopped_{k}": v
                    for k, v in s1.items()
                    if k not in {"n", "nonmissing_n", "missing_n", "missing_percent"}
                },
                **{
                    f"continued_broad_{k}": v
                    for k, v in s0.items()
                    if k not in {"n", "nonmissing_n", "missing_n", "missing_percent"}
                },
            }
        )
        formatted_rows.append(
            {
                "section": section,
                "characteristic": label,
                f"De-escalation or stopping (n={EXPECTED_DEESC:,})": display1,
                f"Continued broad-spectrum (n={EXPECTED_CONT:,})": display0,
                "Absolute SMD before weighting": (
                    f"{smd_before:.3f}" if np.isfinite(smd_before) else ""
                ),
                "Absolute SMD after weighting": (
                    f"{smd_after:.3f}" if np.isfinite(smd_after) else ""
                ),
            }
        )

    metadata = {
        "purpose": "Aggregate baseline characteristics for manuscript/ESM reporting.",
        "source": str(COHORT),
        "cohort_n": EXPECTED_N,
        "deescalation_or_stopping_n": EXPECTED_DEESC,
        "continued_broad_n": EXPECTED_CONT,
        "summary_convention": "Continuous variables: median [IQR] from observed values; categorical variables: n (%) with missing binary values coded as zero to match the frozen PS analysis convention.",
        "smd_convention": "Absolute SMDs before and after the frozen primary stabilized IPTW; continuous missingness follows the prespecified PS preparation for balance calculations.",
        "max_post_weighting_absolute_smd": max_post,
        "worst_post_weighting_variable": worst_post,
        "infection_site_note": "The frozen MIMIC analytic cohort does not contain a validated infection-site variable. Microbiology specimen-category indicators are shown as diagnostic-sampling descriptors and should not be labeled as infection site.",
        "outcomes_analyzed": False,
        "row_level_artifacts_exported": False,
    }
    return pd.DataFrame(detailed_rows), pd.DataFrame(formatted_rows), metadata


def write_markdown(formatted: pd.DataFrame, metadata: dict, path: Path) -> None:
    group1 = f"De-escalation or stopping (n={EXPECTED_DEESC:,})"
    group0 = f"Continued broad-spectrum (n={EXPECTED_CONT:,})"
    lines = [
        "# Candidate baseline characteristics table",
        "",
        "Continuous variables are median [IQR]; categorical variables are n (%). "
        "No hypothesis-test p-values are shown.",
        "",
        f"| Characteristic | {group1} | {group0} | Absolute SMD before | Absolute SMD after |",
        "|---|---:|---:|---:|---:|",
    ]
    last_section = None
    for _, row in formatted.iterrows():
        section = str(row["section"])
        if section != last_section:
            lines.append(f"| **{section}** |  |  |  |  |")
            last_section = section
        lines.append(
            "| {characteristic} | {g1} | {g0} | {pre} | {post} |".format(
                characteristic=row["characteristic"],
                g1=row[group1],
                g0=row[group0],
                pre=row["Absolute SMD before weighting"],
                post=row["Absolute SMD after weighting"],
            )
        )
    lines.extend(
        [
            "",
            "Notes: Absolute standardized mean differences are reported before and after "
            "the primary stabilized inverse probability treatment weighting. Binary "
            "missing values follow the prespecified analysis convention of coding missing "
            "as zero; continuous descriptive summaries use observed values, while the SMD "
            "calculation follows the frozen propensity-score preparation described in Methods.",
            "",
            "Microbiology specimen-category indicators describe diagnostic sampling and are "
            "not validated infection-site labels.",
            "",
            f"Primary balance check: maximum post-weighting absolute SMD "
            f"{metadata['max_post_weighting_absolute_smd']:.3f} "
            f"({metadata['worst_post_weighting_variable']}).",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    if not COHORT.exists():
        raise FileNotFoundError(f"Missing frozen corrected MIMIC cohort: {COHORT}")
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(COHORT, low_memory=False)
    detailed, formatted, metadata = build_table(d)
    detailed.to_csv(OUT / "mimic_baseline_characteristics_detailed.csv", index=False)
    formatted.to_csv(OUT / "mimic_baseline_characteristics_formatted.csv", index=False)
    write_markdown(formatted, metadata, OUT / "mimic_baseline_characteristics.md")
    (OUT / "mimic_baseline_characteristics_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

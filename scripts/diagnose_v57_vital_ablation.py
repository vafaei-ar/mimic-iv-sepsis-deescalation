#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from sepsis_deescalation.features import add_sofa_like
from sepsis_deescalation.specification import CANDIDATE_PS_VARS
from sepsis_deescalation.stats import (
    balance_table,
    effective_sample_size,
    fit_stabilized_iptw,
    risks,
)


TEMP_COLS = ["temp_max_pre72", "temp_max_pre72_missing"]
GCS_COLS = [
    "gcs_total_0_24h",
    "gcs_total_0_24h_missing",
    "gcs_total_48_72h",
    "gcs_total_48_72h_missing",
]
FIO2_COLS = [
    "fio2_0_24h",
    "fio2_0_24h_missing",
    "fio2_48_72h",
    "fio2_48_72h_missing",
]
SOFA_COLS = [
    "cv_score_0_24h",
    "renal_score_0_24h",
    "coag_score_0_24h",
    "liver_score_0_24h",
    "neuro_score_0_24h",
    "resp_score_0_24h",
    "sofa_like_0_24h",
    "cv_score_48_72h",
    "renal_score_48_72h",
    "coag_score_48_72h",
    "liver_score_48_72h",
    "neuro_score_48_72h",
    "resp_score_48_72h",
    "sofa_like_48_72h",
    "sofa_like_change_pre72",
    "sofa_like_improved_pre72",
]


def _copy_cols(dst: pd.DataFrame, src: pd.DataFrame, cols: list[str]) -> None:
    for col in cols:
        if col in src.columns:
            dst[col] = src[col].to_numpy()


def _scenario(
    original: pd.DataFrame,
    corrected: pd.DataFrame,
    *,
    use_temp: bool,
    use_gcs: bool,
    use_fio2: bool,
) -> pd.DataFrame:
    d = original.copy()
    if use_temp:
        _copy_cols(d, corrected, TEMP_COLS)
    if use_gcs:
        _copy_cols(d, corrected, GCS_COLS)
        d = add_sofa_like(d)
    if use_fio2:
        _copy_cols(d, corrected, FIO2_COLS)
    return d


def _summarize(name: str, d: pd.DataFrame, ps_vars: list[str]) -> dict:
    w, _, info = fit_stabilized_iptw(d, ps_vars)
    rt, rc, rd, rr = risks(w, "death_by_horizon", "SW_A")
    bal = balance_table(w, ps_vars)
    if len(bal):
        worst = bal.sort_values("after", ascending=False).iloc[0]
        max_smd = float(worst["after"])
        worst_var = str(worst["variable"])
    else:
        max_smd = np.nan
        worst_var = ""
    wt_t = pd.to_numeric(w.loc[w["A"] == 1, "SW_A"], errors="coerce").dropna()
    wt_c = pd.to_numeric(w.loc[w["A"] == 0, "SW_A"], errors="coerce").dropna()
    ps = pd.to_numeric(w["ps_den"], errors="coerce")
    sw = pd.to_numeric(w["SW_A"], errors="coerce")
    return {
        "analysis": name,
        "risk_deescalated_stopped": rt,
        "risk_continued": rc,
        "risk_difference": rd,
        "risk_ratio": rr,
        "max_post_smd": max_smd,
        "worst_covariate": worst_var,
        "ess_deescalated_stopped": effective_sample_size(wt_t),
        "ess_continued": effective_sample_size(wt_c),
        "max_weight": float(sw.max()),
        "ps_min": float(ps.min()),
        "ps_max": float(ps.max()),
        "den_method": str(info.get("den_method", "")),
        "n_ps_vars_used": len(info.get("used_vars", [])),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Point-estimate ablation diagnostic for repaired v5.7 temperature/GCS/FiO2 covariates.")
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--corrected-cohort", type=Path, default=None)
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    original_path = run_dir / "analysis_cohort_weighted.csv"
    corrected_path = args.corrected_cohort or (run_dir / "audits" / "vital_repair" / "analysis_cohort_vital_corrected.csv")
    if not original_path.exists():
        raise FileNotFoundError(original_path)
    if not corrected_path.exists():
        raise FileNotFoundError(corrected_path)

    original = pd.read_csv(original_path, low_memory=False)
    corrected = pd.read_csv(corrected_path, low_memory=False)
    if len(original) != len(corrected):
        raise ValueError("Original and corrected cohorts have different row counts.")
    if not original["stay_id"].reset_index(drop=True).equals(corrected["stay_id"].reset_index(drop=True)):
        raise ValueError("Original and corrected cohorts are not aligned by stay_id.")

    scenarios = [
        ("A_original_v57", _scenario(original, corrected, use_temp=False, use_gcs=False, use_fio2=False), CANDIDATE_PS_VARS),
        ("B_temperature_only", _scenario(original, corrected, use_temp=True, use_gcs=False, use_fio2=False), CANDIDATE_PS_VARS),
        ("C_temperature_plus_gcs", _scenario(original, corrected, use_temp=True, use_gcs=True, use_fio2=False), CANDIDATE_PS_VARS),
        ("D_temperature_plus_fio2", _scenario(original, corrected, use_temp=True, use_gcs=False, use_fio2=True), CANDIDATE_PS_VARS),
        ("E_all_corrected", _scenario(original, corrected, use_temp=True, use_gcs=True, use_fio2=True), CANDIDATE_PS_VARS),
    ]

    all_corrected = scenarios[-1][1]
    no_fio2 = [v for v in CANDIDATE_PS_VARS if not v.startswith("fio2_")]
    no_gcs = [v for v in CANDIDATE_PS_VARS if not v.startswith("gcs_total_")]
    no_both = [v for v in no_fio2 if not v.startswith("gcs_total_")]
    scenarios.extend([
        ("F_all_corrected_exclude_fio2_ps", all_corrected, no_fio2),
        ("G_all_corrected_exclude_gcs_ps", all_corrected, no_gcs),
        ("H_all_corrected_exclude_gcs_fio2_ps", all_corrected, no_both),
    ])

    rows = [_summarize(name, data, vars_) for name, data, vars_ in scenarios]
    out = pd.DataFrame(rows)
    out_dir = run_dir / "audits" / "vital_ablation"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / "vital_ablation_summary.csv", index=False)
    print(out.to_string(index=False))
    print(f"Outputs: {out_dir}")


if __name__ == "__main__":
    main()

from __future__ import annotations

import numpy as np
import pandas as pd


def add_mortality(cohort: pd.DataFrame, time0_col: str, horizon_days: int) -> pd.DataFrame:
    d = cohort.copy()
    death_dt = pd.to_datetime(d["deathtime"].fillna(d["dod"]), errors="coerce")
    time0 = pd.to_datetime(d[time0_col])
    horizon = time0 + pd.to_timedelta(horizon_days, unit="D")
    d["death_dt"] = death_dt
    d["death_by_horizon"] = (death_dt.notna() & (death_dt >= time0) & (death_dt <= horizon)).astype(int)
    death_days = (death_dt - time0).dt.total_seconds() / 86400.0
    d["followup_time"] = float(horizon_days)
    d.loc[d["death_by_horizon"] == 1, "followup_time"] = death_days.loc[d["death_by_horizon"] == 1]
    d["followup_time"] = d["followup_time"].clip(0, horizon_days)
    return d


def add_hospital_free_days(cohort: pd.DataFrame, time0_col: str, horizon_days: int) -> pd.DataFrame:
    d = cohort.copy()
    los_after = ((pd.to_datetime(d["dischtime"]) - pd.to_datetime(d[time0_col])).dt.total_seconds() / 86400.0).clip(0, horizon_days)
    d["hospital_free_days"] = (horizon_days - los_after).clip(0, horizon_days)
    d.loc[d["death_by_horizon"] == 1, "hospital_free_days"] = 0
    return d


def antibiotic_days(
    cohort: pd.DataFrame,
    rx: pd.DataFrame,
    time0_col: str,
    horizon_days: int,
) -> pd.Series:
    d = cohort.copy()
    coverage = pd.Series(0, index=d.index, dtype=int)
    if len(d) == 0 or len(rx) == 0:
        return coverage
    meds = rx[["hadm_id", "coverage_start", "coverage_stop"]].dropna().copy()
    meds = meds.loc[meds["hadm_id"].isin(set(d["hadm_id"]))]
    base = d[["hadm_id", time0_col]].copy()
    base["_row_index"] = d.index
    for day in range(horizon_days):
        win = base.copy()
        win["window_start"] = pd.to_datetime(win[time0_col]) + pd.to_timedelta(day, unit="D")
        win["window_end"] = pd.to_datetime(win[time0_col]) + pd.to_timedelta(day + 1, unit="D")
        ov = meds.merge(win[["hadm_id", "_row_index", "window_start", "window_end"]], on="hadm_id", how="inner")
        ov = ov.loc[(ov["coverage_start"] < ov["window_end"]) & (ov["coverage_stop"] > ov["window_start"])]
        if len(ov):
            coverage.loc[ov["_row_index"].unique()] += 1
    return coverage


def add_antibiotic_free_days(cohort: pd.DataFrame, rx: pd.DataFrame, time0_col: str, horizon_days: int) -> pd.DataFrame:
    d = cohort.copy()
    d["antibiotic_days"] = antibiotic_days(d, rx, time0_col, horizon_days).reindex(d.index).fillna(0).astype(int)
    d["antibiotic_free_days"] = (horizon_days - d["antibiotic_days"]).clip(0, horizon_days)
    d.loc[d["death_by_horizon"] == 1, "antibiotic_free_days"] = 0
    return d


def add_antibiotic_burden(
    cohort: pd.DataFrame,
    all_rx: pd.DataFrame,
    broad_rx: pd.DataFrame,
    time0_col: str,
    horizon_days: int,
) -> pd.DataFrame:
    d = cohort.copy()
    d["days_alive_30d"] = d["followup_time"].clip(0.1, horizon_days)
    d["normalized_antibiotic_exposure_30d"] = (d["antibiotic_days"] / d["days_alive_30d"]).clip(0, 1)
    d["broad_antibiotic_days_30d"] = antibiotic_days(d, broad_rx, time0_col, horizon_days).reindex(d.index).fillna(0).astype(int)
    d["normalized_broad_antibiotic_exposure_30d"] = (d["broad_antibiotic_days_30d"] / d["days_alive_30d"]).clip(0, 1)

    # Exploratory outcome only. It is vulnerable to differential inpatient observation time.
    meds = all_rx[["hadm_id", "coverage_start", "coverage_stop"]].dropna().copy()
    meds = meds.loc[meds["hadm_id"].isin(set(d["hadm_id"]))]
    base = d[["hadm_id", time0_col]].copy(); base["_row_index"] = d.index
    matrix = {idx: np.zeros(horizon_days, dtype=np.int8) for idx in d.index}
    for day in range(horizon_days):
        win = base.copy(); win["window_start"] = pd.to_datetime(win[time0_col]) + pd.to_timedelta(day, unit="D"); win["window_end"] = pd.to_datetime(win[time0_col]) + pd.to_timedelta(day + 1, unit="D")
        ov = meds.merge(win[["hadm_id", "_row_index", "window_start", "window_end"]], on="hadm_id", how="inner")
        ov = ov.loc[(ov["coverage_start"] < ov["window_end"]) & (ov["coverage_stop"] > ov["window_start"])]
        for idx in ov["_row_index"].unique(): matrix[idx][day] = 1
    flags = []
    for idx in d.index:
        arr = matrix[idx]; flag = any(arr[start:start + 3].sum() == 3 for start in range(7, max(7, horizon_days - 2)))
        flags.append(int(flag))
    d["late_recurrent_or_persistent_abx_course_30d"] = flags
    d["observable_hospital_days_after_landmark"] = ((pd.to_datetime(d["dischtime"]) - pd.to_datetime(d[time0_col])).dt.total_seconds() / 86400.0).clip(0, horizon_days)
    return d

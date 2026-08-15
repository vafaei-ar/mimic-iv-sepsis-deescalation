from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

# Frozen MIMIC v5.5 medication phenotype. Keep changes explicit and versioned.
BROAD_PATTERN = (
    r"vancomycin|piperacillin|tazobactam|zosyn|cefepime|ceftazidime|meropenem|"
    r"imipenem|aztreonam|linezolid|daptomycin|ceftolozane|avibactam"
)
NON_BROAD_PATTERN = (
    r"ceftriaxone|cefazolin|ampicillin|amoxicillin|doxycycline|azithromycin|"
    r"metronidazole|clindamycin|cephalexin|ciprofloxacin|levofloxacin|gentamicin|tobramycin"
)
NON_SYSTEMIC_TEXT = (
    r"oral liquid|enema|ophth|ophthalmic|otic|cream|ointment|topical|lock|irrigation|"
    r"inhal|intrathecal|intravitreal|desensitization|graded challenge|pharmacy to dose"
)
NON_SYSTEMIC_ROUTE = (
    r"rectal|enema|topical|ophth|otic|inhal|lock|irrigation|intravitreal|intrathecal|^it$|^ip$|^im$"
)

ANTI_MRSA_PATTERN = r"vancomycin|linezolid|daptomycin"
ANTIPSEUDOMONAL_PATTERN = (
    r"piperacillin|tazobactam|zosyn|cefepime|ceftazidime|meropenem|imipenem|"
    r"aztreonam|ceftolozane|avibactam"
)
CARBAPENEM_PATTERN = r"meropenem|imipenem"
ANAEROBIC_PATTERN = r"piperacillin|tazobactam|zosyn|meropenem|imipenem|metronidazole|clindamycin"


def prepare_coverage(rx: pd.DataFrame, stop_fill_hours: int = 24) -> pd.DataFrame:
    d = rx.copy()
    d["coverage_start"] = pd.to_datetime(d["starttime"], errors="coerce")
    d["coverage_stop"] = pd.to_datetime(d["stoptime"], errors="coerce")
    missing = d["coverage_stop"].isna()
    d.loc[missing, "coverage_stop"] = d.loc[missing, "coverage_start"] + pd.Timedelta(hours=stop_fill_hours)
    invalid = d["coverage_stop"] <= d["coverage_start"]
    d.loc[invalid, "coverage_stop"] = d.loc[invalid, "coverage_start"] + pd.Timedelta(hours=1)
    d["missing_stop_assumed"] = missing.astype(int)
    return d


def _text_route(df: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    if "drug_lower" in df:
        text = df["drug_lower"].fillna("")
    elif "medication_lower" in df:
        text = df["medication_lower"].fillna("")
    else:
        text = pd.Series("", index=df.index, dtype="object")
    if "route_lower" in df:
        route = df["route_lower"].fillna("").str.strip()
    elif "route" in df:
        route = df["route"].fillna("").str.lower().str.strip()
    else:
        route = pd.Series("", index=df.index, dtype="object")
    return text, route


def systemic_antibiotic_mask(df: pd.DataFrame, broad: bool) -> pd.Series:
    text, route = _text_route(df)
    broad_med = text.str.contains(BROAD_PATTERN, na=False, regex=True)
    non_broad_med = text.str.contains(NON_BROAD_PATTERN, na=False, regex=True)
    bad_text = text.str.contains(NON_SYSTEMIC_TEXT, na=False, regex=True)
    bad_route = route.str.contains(NON_SYSTEMIC_ROUTE, na=False, regex=True)

    if broad:
        clear_iv = route.eq("iv") | route.str.contains("intravenous", na=False)
        return broad_med & clear_iv & ~bad_text & ~bad_route

    any_med = broad_med | non_broad_med
    oral_enteral = route.str.contains(r"po|oral|enteral|ng|og|tube", na=False, regex=True)
    oral_vanco = text.str.contains("vancomycin", na=False) & oral_enteral
    return any_med & ~bad_text & ~bad_route & ~oral_vanco


def raw_broad_mask(df: pd.DataFrame) -> pd.Series:
    text, _ = _text_route(df)
    return text.str.contains(BROAD_PATTERN, na=False, regex=True)


def overlap_rows(
    events: pd.DataFrame,
    id_col: str,
    start_col: str,
    stop_col: str,
    windows: pd.DataFrame,
) -> pd.DataFrame:
    d = events.merge(windows[[id_col, "window_start", "window_end"]], on=id_col, how="inner")
    return d.loc[(d[start_col] < d["window_end"]) & (d[stop_col] > d["window_start"])].copy()


def overlap_hours(row: pd.Series) -> float:
    start = max(row["coverage_start"], row["window_start"])
    stop = min(row["coverage_stop"], row["window_end"])
    return max(0.0, (stop - start).total_seconds() / 3600.0)


def definition_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"domain": "broad-spectrum primary exposure", "rule": "drug text", "details": BROAD_PATTERN},
            {"domain": "broad-spectrum primary exposure", "rule": "route", "details": "IV/intravenous required"},
            {"domain": "broad-spectrum primary exposure", "rule": "excluded text", "details": NON_SYSTEMIC_TEXT},
            {"domain": "broad-spectrum primary exposure", "rule": "excluded route", "details": NON_SYSTEMIC_ROUTE},
            {"domain": "any systemic antibiotic", "rule": "drug text", "details": BROAD_PATTERN + "|" + NON_BROAD_PATTERN},
            {"domain": "any systemic antibiotic", "rule": "oral/enteral non-broad", "details": "allowed"},
            {"domain": "any systemic antibiotic", "rule": "oral/enteral vancomycin", "details": "excluded"},
        ]
    )


def classify_treatment(
    cohort: pd.DataFrame,
    broad_rx: pd.DataFrame,
    all_rx: pd.DataFrame,
    decision_col: str = "decision_time",
    landmark_col: str = "analysis_time0",
) -> pd.DataFrame:
    d = cohort.copy()
    windows = d[["hadm_id", decision_col, landmark_col]].rename(
        columns={decision_col: "window_start", landmark_col: "window_end"}
    )
    broad_post = overlap_rows(broad_rx, "hadm_id", "coverage_start", "coverage_stop", windows)
    any_post = overlap_rows(all_rx, "hadm_id", "coverage_start", "coverage_stop", windows)

    continued = set(broad_post["hadm_id"].dropna().astype(int))
    any_systemic = set(any_post["hadm_id"].dropna().astype(int))
    d["A"] = (~d["hadm_id"].astype(int).isin(continued)).astype(int)
    d["has_any_systemic_72_96"] = d["hadm_id"].astype(int).isin(any_systemic).astype(int)
    d["deescalation_type"] = "continued_broad"
    d.loc[(d["A"] == 1) & (d["has_any_systemic_72_96"] == 1), "deescalation_type"] = "narrowed_or_non_broad_only"
    d.loc[(d["A"] == 1) & (d["has_any_systemic_72_96"] == 0), "deescalation_type"] = "stopped_all_observed_systemic_antibiotics"

    if len(broad_post):
        broad_post = broad_post.copy()
        broad_post["overlap_hours"] = broad_post.apply(overlap_hours, axis=1)
        overlap_sum = broad_post.groupby("hadm_id")["overlap_hours"].sum()
        d["broad_overlap_hours_72_96"] = d["hadm_id"].map(overlap_sum).fillna(0.0)
    else:
        d["broad_overlap_hours_72_96"] = 0.0
    return d


def predecision_antibiotic_intensity(
    cohort: pd.DataFrame,
    all_rx: pd.DataFrame,
    broad_rx: pd.DataFrame,
) -> pd.DataFrame:
    d = cohort.copy()
    windows = d[["hadm_id", "first_broad_time", "decision_time"]].rename(
        columns={"first_broad_time": "window_start", "decision_time": "window_end"}
    )
    broad = overlap_rows(broad_rx, "hadm_id", "coverage_start", "coverage_stop", windows)
    any_abx = overlap_rows(all_rx, "hadm_id", "coverage_start", "coverage_stop", windows)

    if len(broad):
        broad = broad.copy()
        broad["overlap_hours"] = broad.apply(overlap_hours, axis=1)
        hours = broad.groupby("hadm_id")["overlap_hours"].sum()
        agents = broad.groupby("hadm_id")["drug_lower"].nunique()
    else:
        hours = pd.Series(dtype=float)
        agents = pd.Series(dtype=float)
    d["broad_hours_pre72"] = d["hadm_id"].map(hours).fillna(0.0)
    d["broad_agents_pre72"] = d["hadm_id"].map(agents).fillna(0.0)

    text = any_abx.get("drug_lower", pd.Series("", index=any_abx.index)).fillna("")
    for name, pattern in {
        "anti_mrsa_pre72": ANTI_MRSA_PATTERN,
        "antipseudomonal_pre72": ANTIPSEUDOMONAL_PATTERN,
        "carbapenem_pre72": CARBAPENEM_PATTERN,
        "anaerobic_pre72": ANAEROBIC_PATTERN,
    }.items():
        ids = set(any_abx.loc[text.str.contains(pattern, na=False, regex=True), "hadm_id"].astype(int))
        d[name] = d["hadm_id"].astype(int).isin(ids).astype(int)
    return d

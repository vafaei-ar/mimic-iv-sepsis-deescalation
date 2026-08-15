from __future__ import annotations

import numpy as np
import pandas as pd

from .antibiotics import overlap_rows
from .mimic_io import read_csv, read_csv_filtered


def fill_numeric_with_median(d: pd.DataFrame, col: str, default: float = 0.0) -> pd.DataFrame:
    d[col] = pd.to_numeric(d[col], errors="coerce")
    d[f"{col}_missing"] = d[col].isna().astype(int)
    med = d[col].median()
    if pd.isna(med):
        med = default
    d[col] = d[col].fillna(med)
    return d


def select_d_items(source, concept_patterns: dict[str, str]) -> pd.DataFrame:
    items = read_csv(source, "icu/d_items.csv.gz", usecols=["itemid", "label"])
    items["label_lower"] = items["label"].fillna("").str.lower()
    parts = []
    for concept, pattern in concept_patterns.items():
        sub = items.loc[items["label_lower"].str.contains(pattern, na=False, regex=True)].copy()
        sub["concept"] = concept
        parts.append(sub)
    if not parts:
        return pd.DataFrame(columns=["itemid", "label", "concept"])
    return pd.concat(parts, ignore_index=True).drop_duplicates(["itemid", "concept"])[["itemid", "label", "concept"]]


def lab_items(source, labels: list[str]) -> pd.DataFrame:
    items = read_csv(source, "hosp/d_labitems.csv.gz", usecols=["itemid", "label", "fluid", "category"])
    items["label_lower"] = items["label"].fillna("").str.lower().str.strip()
    wanted = {x.lower().strip() for x in labels}
    return items.loc[items["label_lower"].isin(wanted)].copy()


def add_baseline_covariates(cohort: pd.DataFrame, dx: pd.DataFrame, proc: pd.DataFrame) -> pd.DataFrame:
    d = cohort.copy()
    dxs = dx.loc[dx["hadm_id"].isin(d["hadm_id"])].copy()
    px = proc.loc[proc["hadm_id"].isin(d["hadm_id"])].copy()
    chronic = r"diabetes|chronic kidney|renal failure|heart failure|copd|chronic obstructive|liver disease|cirrhosis|malignan|cancer"
    d["age"] = pd.to_numeric(d["anchor_age"], errors="coerce")
    d["sex_male"] = (d["gender"] == "M").astype(int)
    d["race_white"] = d["race"].fillna("").str.lower().str.contains("white").astype(int)
    d["comorb"] = d["hadm_id"].isin(set(dxs.loc[dxs["long_title_lower"].str.contains(chronic, na=False), "hadm_id"])).astype(int)
    d["heart_failure"] = d["hadm_id"].isin(set(dxs.loc[dxs["long_title_lower"].str.contains("heart failure", na=False), "hadm_id"])).astype(int)
    d["chronic_kidney"] = d["hadm_id"].isin(set(dxs.loc[dxs["long_title_lower"].str.contains("chronic kidney|renal failure", na=False), "hadm_id"])).astype(int)
    d["vent_proc"] = d["hadm_id"].isin(set(px.loc[px["long_title_lower"].str.contains("mechanical ventilation|respiratory ventilation|intubation", na=False), "hadm_id"])).astype(int)
    careunit = d["first_careunit"].fillna("").str.lower()
    d["micu"] = careunit.str.contains("medical").astype(int)
    d["sicu"] = careunit.str.contains("surgical").astype(int)
    d["cardiac_icu"] = careunit.str.contains("coronary|cardiac|cardiovascular").astype(int)
    d["neuro_icu"] = careunit.str.contains("neuro").astype(int)
    d["hours_admit_to_icu"] = ((pd.to_datetime(d["intime"]) - pd.to_datetime(d["admittime"])).dt.total_seconds() / 3600).clip(-24, 24 * 30).fillna(0)
    return d


def add_baseline_labs(cohort: pd.DataFrame, source) -> pd.DataFrame:
    d = cohort.copy()
    labels = ["Lactate", "Creatinine", "White Blood Cells"]
    items = lab_items(source, labels)
    ids = set(d["hadm_id"].dropna().astype(int))
    itemids = set(items["itemid"].astype(int))
    labs = read_csv_filtered(
        source,
        "hosp/labevents.csv.gz",
        usecols=["subject_id", "hadm_id", "itemid", "charttime", "valuenum"],
        parse_dates=["charttime"],
        filter_func=lambda c: c["hadm_id"].isin(ids) & c["itemid"].isin(itemids),
    )
    labs = labs.merge(items[["itemid", "label"]], on="itemid", how="left").merge(
        d[["hadm_id", "intime", "decision_time"]], on="hadm_id", how="inner"
    )
    labs = labs.loc[(labs["charttime"] >= labs["intime"]) & (labs["charttime"] <= labs["decision_time"]) & labs["valuenum"].notna()].copy()
    for label in labels:
        base = label.lower().replace(" ", "_")
        sub = labs.loc[labs["label"].str.lower().str.strip() == label.lower()].sort_values(["hadm_id", "charttime"])
        if len(sub):
            last = sub.drop_duplicates("hadm_id", keep="last").set_index("hadm_id")["valuenum"]
            first = sub.drop_duplicates("hadm_id", keep="first").set_index("hadm_id")["valuenum"]
            d[f"{base}_last_pre72"] = d["hadm_id"].map(last)
            d[f"{base}_change_pre72"] = d["hadm_id"].map(last - first)
        else:
            d[f"{base}_last_pre72"] = np.nan
            d[f"{base}_change_pre72"] = np.nan
        d[f"{base}_missing_pre72"] = d[f"{base}_last_pre72"].isna().astype(int)
        d[f"{base}_last_pre72"] = d[f"{base}_last_pre72"].fillna(d[f"{base}_last_pre72"].median() if d[f"{base}_last_pre72"].notna().any() else 0)
        d[f"{base}_change_pre72"] = d[f"{base}_change_pre72"].fillna(0)
    return d


def add_baseline_vitals(cohort: pd.DataFrame, source) -> pd.DataFrame:
    d = cohort.copy()
    patterns = {
        "heart_rate": r"^heart rate$",
        "resp_rate": r"respiratory rate",
        "spo2": r"o2 saturation|spo2|oxygen saturation",
        "map": r"blood pressure mean|arterial blood pressure mean|non invasive blood pressure mean",
        "temperature": r"temperature",
    }
    items = select_d_items(source, patterns)
    stay_ids = set(d["stay_id"].dropna().astype(int))
    itemids = set(items["itemid"].astype(int))
    ce = read_csv_filtered(
        source,
        "icu/chartevents.csv.gz",
        usecols=["subject_id", "hadm_id", "stay_id", "charttime", "itemid", "valuenum"],
        parse_dates=["charttime"],
        filter_func=lambda c: c["stay_id"].isin(stay_ids) & c["itemid"].isin(itemids),
        chunksize=500_000,
    )
    ce = ce.merge(items, on="itemid", how="inner").merge(d[["stay_id", "intime", "decision_time"]], on="stay_id", how="inner")
    ce = ce.loc[(ce["charttime"] >= ce["intime"]) & (ce["charttime"] <= ce["decision_time"]) & ce["valuenum"].notna()].copy()
    ranges = {"heart_rate": (20, 250), "resp_rate": (1, 90), "spo2": (20, 100), "map": (20, 200)}
    for concept, (lo, hi) in ranges.items():
        ce.loc[(ce["concept"] == concept) & ((ce["valuenum"] < lo) | (ce["valuenum"] > hi)), "valuenum"] = np.nan
    tmask = ce["concept"] == "temperature"
    ce.loc[tmask & (ce["valuenum"] > 70), "valuenum"] = (ce.loc[tmask & (ce["valuenum"] > 70), "valuenum"] - 32) * 5 / 9
    ce.loc[tmask & ((ce["valuenum"] < 25) | (ce["valuenum"] > 45)), "valuenum"] = np.nan
    specs = {"heart_rate": ("max", "hr_max_pre72"), "resp_rate": ("max", "rr_max_pre72"), "spo2": ("min", "spo2_min_pre72"), "map": ("min", "map_min_pre72"), "temperature": ("max", "temp_max_pre72")}
    for concept, (fn, out) in specs.items():
        sub = ce.loc[ce["concept"] == concept]
        vals = sub.groupby("stay_id")["valuenum"].agg(fn) if len(sub) else pd.Series(dtype=float)
        d[out] = d["stay_id"].map(vals)
        d[f"{out}_missing"] = d[out].isna().astype(int)
        d[out] = d[out].fillna(d[out].median() if d[out].notna().any() else 0)
    return d


def prepare_vasopressors(cohort: pd.DataFrame, inp: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = cohort.copy()
    vaso = inp.loc[inp["label_lower"].str.contains("norepinephrine|phenylephrine|vasopressin|epinephrine|dopamine", na=False) & inp["starttime"].notna()].copy()
    vaso["vaso_start"] = pd.to_datetime(vaso["starttime"])
    vaso["vaso_stop"] = pd.to_datetime(vaso["endtime"])
    missing = vaso["vaso_stop"].isna()
    vaso.loc[missing, "vaso_stop"] = vaso.loc[missing, "vaso_start"] + pd.Timedelta(hours=1)
    bad = vaso["vaso_stop"] <= vaso["vaso_start"]
    vaso.loc[bad, "vaso_stop"] = vaso.loc[bad, "vaso_start"] + pd.Timedelta(hours=1)
    pre = d[["stay_id", "intime", "decision_time"]].rename(columns={"intime": "window_start", "decision_time": "window_end"})
    ov = overlap_rows(vaso, "stay_id", "vaso_start", "vaso_stop", pre)
    d["vasopressor_any_pre72"] = d["stay_id"].isin(set(ov["stay_id"])).astype(int)
    if len(ov):
        ov["overlap_start"] = ov[["vaso_start", "window_start"]].max(axis=1)
        ov["overlap_stop"] = ov[["vaso_stop", "window_end"]].min(axis=1)
        ov["vaso_hours"] = (ov["overlap_stop"] - ov["overlap_start"]).dt.total_seconds() / 3600
        d["vasopressor_hours_pre72"] = d["stay_id"].map(ov.groupby("stay_id")["vaso_hours"].sum()).fillna(0).clip(0, 72)
    else:
        d["vasopressor_hours_pre72"] = 0.0
    late = d[["stay_id", "decision_time"]].rename(columns={"decision_time": "window_end"})
    late["window_start"] = late["window_end"] - pd.Timedelta(hours=6)
    late_ov = overlap_rows(vaso, "stay_id", "vaso_start", "vaso_stop", late)
    d["vasopressor_overlap_last6h_pre72"] = d["stay_id"].isin(set(late_ov["stay_id"])).astype(int)
    return d, vaso


def add_vasopressor_windows(d: pd.DataFrame, vaso: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    early = d[["stay_id", "first_broad_time"]].copy(); early["window_start"] = early["first_broad_time"]; early["window_end"] = early["first_broad_time"] + pd.Timedelta(hours=24)
    late = d[["stay_id", "decision_time"]].copy(); late["window_start"] = late["decision_time"] - pd.Timedelta(hours=24); late["window_end"] = late["decision_time"]
    v0 = overlap_rows(vaso, "stay_id", "vaso_start", "vaso_stop", early) if len(vaso) else pd.DataFrame()
    v1 = overlap_rows(vaso, "stay_id", "vaso_start", "vaso_stop", late) if len(vaso) else pd.DataFrame()
    d["vasopressor_any_0_24h"] = d["stay_id"].isin(set(v0.get("stay_id", []))).astype(int)
    d["vasopressor_any_48_72h"] = d["stay_id"].isin(set(v1.get("stay_id", []))).astype(int)
    d["vasopressor_stopped_before_72h"] = ((d["vasopressor_any_0_24h"] == 1) & (d["vasopressor_any_48_72h"] == 0)).astype(int)
    return d


def add_lab_trajectories(d: pd.DataFrame, source) -> pd.DataFrame:
    d = d.copy()
    labels = ["Lactate", "Creatinine", "White Blood Cells", "Platelet Count", "Bilirubin, Total"]
    mapping = {"Lactate": "lactate", "Creatinine": "creatinine", "White Blood Cells": "wbc", "Platelet Count": "platelet", "Bilirubin, Total": "bilirubin"}
    items = lab_items(source, labels)
    ids, itemids = set(d["hadm_id"].dropna().astype(int)), set(items["itemid"].astype(int))
    labs = read_csv_filtered(source, "hosp/labevents.csv.gz", usecols=["subject_id", "hadm_id", "itemid", "charttime", "valuenum"], parse_dates=["charttime"], filter_func=lambda c: c["hadm_id"].isin(ids) & c["itemid"].isin(itemids))
    labs = labs.merge(items[["itemid", "label"]], on="itemid").merge(d[["hadm_id", "first_broad_time", "decision_time"]], on="hadm_id")
    labs = labs.loc[labs["valuenum"].notna()].copy()
    labs["early"] = (labs["charttime"] >= labs["first_broad_time"]) & (labs["charttime"] <= labs["first_broad_time"] + pd.Timedelta(hours=24))
    labs["late"] = (labs["charttime"] >= labs["decision_time"] - pd.Timedelta(hours=24)) & (labs["charttime"] <= labs["decision_time"])
    plaus = {"Lactate": (0.1, 30), "Creatinine": (0.1, 30), "White Blood Cells": (0.1, 200), "Platelet Count": (1, 2000), "Bilirubin, Total": (0.01, 80)}
    for label, (lo, hi) in plaus.items():
        m = labs["label"].str.lower().str.strip() == label.lower(); labs.loc[m & ((labs["valuenum"] < lo) | (labs["valuenum"] > hi)), "valuenum"] = np.nan
    for label, base in mapping.items():
        sub = labs.loc[labs["label"].str.lower().str.strip() == label.lower()].copy()
        for flag, suffix in [("early", "0_24h"), ("late", "48_72h")]:
            win = sub.loc[sub[flag]].sort_values(["hadm_id", "charttime"])
            last_name = f"{base}_{'early_last' if flag == 'early' else 'late_last'}_{suffix}"
            worst_name = f"{base}_{'early_worst' if flag == 'early' else 'late_worst'}_{suffix}"
            d[last_name] = d["hadm_id"].map(win.drop_duplicates("hadm_id", keep="last").set_index("hadm_id")["valuenum"] if len(win) else pd.Series(dtype=float))
            agg = "min" if base == "platelet" else "max"
            d[worst_name] = d["hadm_id"].map(win.groupby("hadm_id")["valuenum"].agg(agg) if len(win) else pd.Series(dtype=float))
        d[f"{base}_change_early_to_late"] = d[f"{base}_late_last_48_72h"] - d[f"{base}_early_last_0_24h"]
        for col in [f"{base}_early_last_0_24h", f"{base}_early_worst_0_24h", f"{base}_late_last_48_72h", f"{base}_late_worst_48_72h", f"{base}_change_early_to_late"]:
            d = fill_numeric_with_median(d, col)
    d["lactate_rising_pre72"] = (d["lactate_change_early_to_late"] > 0).astype(int)
    d["wbc_rising_pre72"] = (d["wbc_change_early_to_late"] > 0).astype(int)
    d["creatinine_rising_pre72"] = (d["creatinine_change_early_to_late"] > 0).astype(int)
    d["platelet_falling_pre72"] = (d["platelet_change_early_to_late"] < 0).astype(int)
    d["bilirubin_rising_pre72"] = (d["bilirubin_change_early_to_late"] > 0).astype(int)
    return d


def add_vital_trajectories(d: pd.DataFrame, source) -> pd.DataFrame:
    d = d.copy()
    patterns = {"heart_rate": r"^heart rate$", "resp_rate": r"respiratory rate", "spo2": r"o2 saturation|spo2|oxygen saturation", "map": r"blood pressure mean|arterial blood pressure mean|non invasive blood pressure mean", "temperature": r"temperature", "gcs_total": r"gcs.*total|glasgow coma scale.*total", "fio2": r"fraction inspired oxygen|fio2"}
    items = select_d_items(source, patterns)
    stay_ids, itemids = set(d["stay_id"].dropna().astype(int)), set(items["itemid"].astype(int))
    ce = read_csv_filtered(source, "icu/chartevents.csv.gz", usecols=["subject_id", "hadm_id", "stay_id", "charttime", "itemid", "valuenum"], parse_dates=["charttime"], filter_func=lambda c: c["stay_id"].isin(stay_ids) & c["itemid"].isin(itemids), chunksize=500_000)
    ce = ce.merge(items, on="itemid").merge(d[["stay_id", "first_broad_time", "decision_time"]], on="stay_id")
    ce = ce.loc[ce["valuenum"].notna()].copy()
    ce["early"] = (ce["charttime"] >= ce["first_broad_time"]) & (ce["charttime"] <= ce["first_broad_time"] + pd.Timedelta(hours=24))
    ce["late"] = (ce["charttime"] >= ce["decision_time"] - pd.Timedelta(hours=24)) & (ce["charttime"] <= ce["decision_time"])
    ce["last12"] = (ce["charttime"] >= ce["decision_time"] - pd.Timedelta(hours=12)) & (ce["charttime"] <= ce["decision_time"])
    ranges = {"heart_rate": (20, 250), "resp_rate": (1, 90), "spo2": (20, 100), "map": (20, 200), "gcs_total": (3, 15), "fio2": (0.2, 100)}
    for concept, (lo, hi) in ranges.items():
        ce.loc[(ce["concept"] == concept) & ((ce["valuenum"] < lo) | (ce["valuenum"] > hi)), "valuenum"] = np.nan
    tm = ce["concept"] == "temperature"; ce.loc[tm & (ce["valuenum"] > 70), "valuenum"] = (ce.loc[tm & (ce["valuenum"] > 70), "valuenum"] - 32) * 5 / 9; ce.loc[tm & ((ce["valuenum"] < 25) | (ce["valuenum"] > 45)), "valuenum"] = np.nan
    fm = ce["concept"] == "fio2"; ce.loc[fm & (ce["valuenum"] > 1.5), "valuenum"] = ce.loc[fm & (ce["valuenum"] > 1.5), "valuenum"] / 100
    aggregation = {"heart_rate": "max", "resp_rate": "max", "spo2": "min", "map": "min", "temperature": "max", "gcs_total": "min", "fio2": "max"}
    for concept, fn in aggregation.items():
        sub = ce.loc[ce["concept"] == concept]
        for flag, suffix in [("early", "0_24h"), ("late", "48_72h")]:
            win = sub.loc[sub[flag]]
            d[f"{concept}_{suffix}"] = d["stay_id"].map(win.groupby("stay_id")["valuenum"].agg(fn) if len(win) else pd.Series(dtype=float))
            d = fill_numeric_with_median(d, f"{concept}_{suffix}")
    fever_ids = set(ce.loc[(ce["concept"] == "temperature") & ce["last12"] & (ce["valuenum"] >= 38), "stay_id"])
    d["fever_last12h_pre72"] = d["stay_id"].isin(fever_ids).astype(int)
    d["map_improved_pre72"] = (d["map_48_72h"] > d["map_0_24h"]).astype(int)
    d["spo2_improved_pre72"] = (d["spo2_48_72h"] > d["spo2_0_24h"]).astype(int)
    d["temp_improved_pre72"] = (d["temperature_48_72h"] < d["temperature_0_24h"]).astype(int)
    d["rr_improved_pre72"] = (d["resp_rate_48_72h"] < d["resp_rate_0_24h"]).astype(int)
    return d


def add_urine_output(d: pd.DataFrame, source) -> pd.DataFrame:
    d = d.copy()
    items = read_csv(source, "icu/d_items.csv.gz", usecols=["itemid", "label"]); items["label_lower"] = items["label"].fillna("").str.lower(); items = items.loc[items["label_lower"].str.contains("urine|foley|void|urinary", na=False)]
    stay_ids, itemids = set(d["stay_id"].dropna().astype(int)), set(items["itemid"].astype(int))
    out = read_csv_filtered(source, "icu/outputevents.csv.gz", usecols=["subject_id", "hadm_id", "stay_id", "charttime", "itemid", "value"], parse_dates=["charttime"], filter_func=lambda c: c["stay_id"].isin(stay_ids) & c["itemid"].isin(itemids))
    out["value_num"] = pd.to_numeric(out["value"], errors="coerce"); out = out.loc[(out["value_num"] >= 0) & (out["value_num"] <= 10000)].merge(d[["stay_id", "first_broad_time", "decision_time"]], on="stay_id")
    early = out.loc[(out["charttime"] >= out["first_broad_time"]) & (out["charttime"] <= out["first_broad_time"] + pd.Timedelta(hours=24))].groupby("stay_id")["value_num"].sum()
    late = out.loc[(out["charttime"] >= out["decision_time"] - pd.Timedelta(hours=24)) & (out["charttime"] <= out["decision_time"])].groupby("stay_id")["value_num"].sum()
    d["urine_output_ml_0_24h"] = d["stay_id"].map(early); d["urine_output_ml_48_72h"] = d["stay_id"].map(late); d["urine_output_change_pre72"] = d["urine_output_ml_48_72h"] - d["urine_output_ml_0_24h"]
    for col in ["urine_output_ml_0_24h", "urine_output_ml_48_72h", "urine_output_change_pre72"]: d = fill_numeric_with_median(d, col)
    d["low_urine_output_48_72h"] = (d["urine_output_ml_48_72h"] < 500).astype(int)
    return d


def add_steroid_bmi(d: pd.DataFrame, rx: pd.DataFrame, source) -> pd.DataFrame:
    d = d.copy(); pre = d[["hadm_id", "first_broad_time", "decision_time"]].copy(); pre["window_start"] = pre["first_broad_time"]; pre["window_end"] = pre["decision_time"]
    steroid = rx.loc[rx["drug_lower"].str.contains("hydrocortisone|methylprednisolone|prednisone|prednisolone|dexamethasone", na=False) & ~rx["route_lower"].str.contains("topical|ophth|otic|inhal", na=False)]
    ov = overlap_rows(steroid, "hadm_id", "coverage_start", "coverage_stop", pre) if len(steroid) else pd.DataFrame()
    d["steroid_any_pre72"] = d["hadm_id"].isin(set(ov.get("hadm_id", []))).astype(int); d["hydrocortisone_any_pre72"] = d["hadm_id"].isin(set(ov.loc[ov.get("drug_lower", pd.Series(dtype=str)).str.contains("hydrocortisone", na=False), "hadm_id"] if len(ov) else [])).astype(int)
    try:
        items = select_d_items(source, {"weight_kg": r"admission weight|daily weight|weight", "height_cm": r"height"}); stay_ids, itemids = set(d["stay_id"].dropna().astype(int)), set(items["itemid"].astype(int))
        ce = read_csv_filtered(source, "icu/chartevents.csv.gz", usecols=["subject_id", "hadm_id", "stay_id", "charttime", "itemid", "valuenum"], parse_dates=["charttime"], filter_func=lambda c: c["stay_id"].isin(stay_ids) & c["itemid"].isin(itemids), chunksize=500_000).merge(items, on="itemid").merge(d[["stay_id", "intime", "decision_time"]], on="stay_id")
        ce = ce.loc[(ce["charttime"] >= ce["intime"]) & (ce["charttime"] <= ce["decision_time"]) & ce["valuenum"].notna()]
        h = ce.loc[ce["concept"] == "height_cm"].copy(); h.loc[(h["valuenum"] > 40) & (h["valuenum"] < 90), "valuenum"] *= 2.54; h = h.loc[(h["valuenum"] >= 120) & (h["valuenum"] <= 230)].sort_values(["stay_id", "charttime"]).drop_duplicates("stay_id").set_index("stay_id")["valuenum"]
        w = ce.loc[(ce["concept"] == "weight_kg") & (ce["valuenum"] >= 25) & (ce["valuenum"] <= 300)].sort_values(["stay_id", "charttime"]).drop_duplicates("stay_id").set_index("stay_id")["valuenum"]
        d["height_cm_pre72"] = d["stay_id"].map(h); d["weight_kg_pre72"] = d["stay_id"].map(w); d["bmi_pre72"] = d["weight_kg_pre72"] / ((d["height_cm_pre72"] / 100) ** 2)
    except Exception:
        d["bmi_pre72"] = np.nan
    return fill_numeric_with_median(d, "bmi_pre72", 30.0)


def add_sofa_like(d: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    creat = lambda x: np.select([x < 1.2, x < 2, x < 3.5, x < 5, x >= 5], [0, 1, 2, 3, 4], default=0)
    bili = lambda x: np.select([x < 1.2, x < 2, x < 6, x < 12, x >= 12], [0, 1, 2, 3, 4], default=0)
    plate = lambda x: np.select([x >= 150, x >= 100, x >= 50, x >= 20, x < 20], [0, 1, 2, 3, 4], default=0)
    gcs = lambda x: np.select([x >= 15, x >= 13, x >= 10, x >= 6, x < 6], [0, 1, 2, 3, 4], default=0)
    for suffix, phase in [("0_24h", "early"), ("48_72h", "late")]:
        d[f"cv_score_{suffix}"] = np.where(d[f"vasopressor_any_{suffix}"] == 1, 3, np.where(d[f"map_{suffix}"] < 70, 1, 0))
        creat_col = "creatinine_early_last_0_24h" if suffix == "0_24h" else "creatinine_late_last_48_72h"
        urine_col = "urine_output_ml_0_24h" if suffix == "0_24h" else "urine_output_ml_48_72h"
        plate_col = "platelet_early_worst_0_24h" if suffix == "0_24h" else "platelet_late_worst_48_72h"
        bili_col = "bilirubin_early_worst_0_24h" if suffix == "0_24h" else "bilirubin_late_worst_48_72h"
        d[f"renal_score_{suffix}"] = np.maximum(creat(d[creat_col]), np.where(d[urine_col] < 500, 3, 0)); d[f"coag_score_{suffix}"] = plate(d[plate_col]); d[f"liver_score_{suffix}"] = bili(d[bili_col]); d[f"neuro_score_{suffix}"] = gcs(d[f"gcs_total_{suffix}"]); d[f"resp_score_{suffix}"] = np.where((d[f"spo2_{suffix}"] < 90) | (d["vent_proc"] == 1), 2, np.where(d[f"spo2_{suffix}"] < 94, 1, 0))
        components = [f"cv_score_{suffix}", f"renal_score_{suffix}", f"coag_score_{suffix}", f"liver_score_{suffix}", f"neuro_score_{suffix}", f"resp_score_{suffix}"]
        d[f"sofa_like_{suffix}"] = d[components].sum(axis=1)
    d["sofa_like_change_pre72"] = d["sofa_like_48_72h"] - d["sofa_like_0_24h"]; d["sofa_like_improved_pre72"] = (d["sofa_like_change_pre72"] < 0).astype(int)
    return d


def add_severity_score(d: pd.DataFrame) -> pd.DataFrame:
    d = d.copy(); terms = []
    for col in ["lactate_last_pre72", "creatinine_last_pre72", "white_blood_cells_last_pre72", "hr_max_pre72", "rr_max_pre72", "map_min_pre72", "spo2_min_pre72", "vasopressor_hours_pre72"]:
        x = pd.to_numeric(d[col], errors="coerce"); sd = x.std(ddof=0)
        if sd and np.isfinite(sd): terms.append((-(x - x.mean()) / sd if col in ["map_min_pre72", "spo2_min_pre72"] else (x - x.mean()) / sd).fillna(0))
    d["severity_pre72"] = np.sum(terms, axis=0) if terms else 0.0; d["severity_pre72"] = d["severity_pre72"] + 0.3 * d["vent_proc"] + 0.2 * d["comorb"]
    return d

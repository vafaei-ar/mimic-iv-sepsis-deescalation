from __future__ import annotations

"""Audited MIMIC vital-sign corrections used by the final v5.7 analysis.

The original scripted pipeline used broad D_ITEMS label matching for temperature,
GCS total, and FiO2. A targeted real-data audit found three problems:

* baseline temperature could mix Fahrenheit and Celsius before aggregation;
* the direct GCS-total extraction was not a trustworthy routine measurement;
* direct FiO2 matching was sparse/ambiguous.

The publication analysis repaired those variables from validated item definitions before
refitting the propensity score. This module moves that repair into reusable pipeline code
so a clean from-scratch run produces the same corrected feature semantics rather than
requiring a one-off post-processing script.

Direct GCS and FiO2 are still excluded from the primary propensity model by
``specification.PS_EXCLUDED_VITAL_DIRECT_TERMS``. Reconstructed GCS is retained because
it contributes to the SOFA-like trajectory that is calculated later in the pipeline.
"""

import numpy as np
import pandas as pd

from .mimic_io import read_csv, read_csv_filtered


AUDITED_TEMPERATURE_LABELS = {
    "temperature fahrenheit",
    "temperature celsius",
    "blood temperature cco (c)",
    "cerebral temperature (c)",
}
AUDITED_FIO2_LABEL_PATTERN = r"inspired o2 fraction|fraction inspired oxygen"
GCS_COMPONENT_MAP = {
    "gcs - eye opening": "eye",
    "gcs - verbal response": "verbal",
    "gcs - motor response": "motor",
}


def temperature_to_celsius(values: pd.Series, labels: pd.Series) -> pd.Series:
    """Normalize audited temperature readings before patient-level aggregation.

    Label semantics are preferred. The value >60 fallback handles generic temperature
    labels whose numeric scale is clearly Fahrenheit. Conversion is deliberately done
    at the reading level; converting a patient-level maximum after aggregation can be
    wrong when one stay contains both Celsius and Fahrenheit rows.
    """
    v = pd.to_numeric(values, errors="coerce")
    lab = labels.fillna("").astype(str).str.lower()
    is_f = lab.str.contains("fahrenheit|temperaturef", regex=True, na=False)
    is_c = lab.str.contains(r"celsius|\(c\)", regex=True, na=False)
    generic_f = (~is_f & ~is_c) & (v > 60)
    out = v.copy()
    out.loc[is_f | generic_f] = (v.loc[is_f | generic_f] - 32.0) * 5.0 / 9.0
    return out


def _median_impute_with_flag(d: pd.DataFrame, col: str) -> None:
    """Match the historical feature-layer imputation used before PS preparation."""
    x = pd.to_numeric(d[col], errors="coerce")
    d[f"{col}_missing"] = x.isna().astype(int)
    med = x.median()
    d[col] = x.fillna(0.0 if pd.isna(med) else med)


def _select_audited_items(source) -> pd.DataFrame:
    items = read_csv(source, "icu/d_items.csv.gz", usecols=["itemid", "label", "unitname"])
    items["label_lower"] = items["label"].fillna("").astype(str).str.lower().str.strip()

    temp = items.loc[items["label_lower"].isin(AUDITED_TEMPERATURE_LABELS)].copy()
    temp["domain"] = "temperature"

    fio2 = items.loc[items["label_lower"].str.fullmatch(AUDITED_FIO2_LABEL_PATTERN, na=False)].copy()
    fio2["domain"] = "fio2"

    gcs = items.loc[items["label_lower"].isin(GCS_COMPONENT_MAP)].copy()
    gcs["domain"] = "gcs"
    gcs["gcs_component"] = gcs["label_lower"].map(GCS_COMPONENT_MAP)

    if temp.empty:
        raise RuntimeError("Audited MIMIC temperature items were not found in d_items")
    if gcs.empty:
        raise RuntimeError("Audited MIMIC GCS component items were not found in d_items")
    if fio2.empty:
        raise RuntimeError("Audited MIMIC FiO2 items were not found in d_items")

    return pd.concat([temp, fio2, gcs], ignore_index=True, sort=False)


def _load_events(source, cohort: pd.DataFrame, selected: pd.DataFrame, chunksize: int) -> pd.DataFrame:
    stay_ids = set(cohort["stay_id"].dropna().astype(int))
    itemids = set(selected["itemid"].dropna().astype(int))
    events = read_csv_filtered(
        source,
        "icu/chartevents.csv.gz",
        usecols=["stay_id", "charttime", "itemid", "valuenum"],
        parse_dates=["charttime"],
        filter_func=lambda c: c["stay_id"].isin(stay_ids) & c["itemid"].isin(itemids),
        chunksize=chunksize,
    )
    events = events.merge(
        selected[[c for c in ["itemid", "label", "unitname", "domain", "gcs_component"] if c in selected]],
        on="itemid",
        how="inner",
    )
    windows = cohort[["stay_id", "intime", "first_broad_time", "decision_time"]].drop_duplicates("stay_id")
    events = events.merge(windows, on="stay_id", how="inner")
    events["charttime"] = pd.to_datetime(events["charttime"], errors="coerce")
    return events.loc[events["valuenum"].notna() & events["charttime"].notna()].copy()


def _correct_baseline_temperature(d: pd.DataFrame, events: pd.DataFrame) -> dict:
    # This reproduces the accepted v5.7 repair: baseline temperature spans ICU
    # intime through the 72-h decision, is converted reading-by-reading to C,
    # then takes the within-stay maximum.
    t = events.loc[events["domain"] == "temperature"].copy()
    t = t.loc[(t["charttime"] >= t["intime"]) & (t["charttime"] <= t["decision_time"])].copy()
    t["temperature_c"] = temperature_to_celsius(t["valuenum"], t["label"])
    t = t.loc[t["temperature_c"].between(30.0, 45.0, inclusive="both")].copy()
    vals = t.groupby("stay_id")["temperature_c"].max()
    d["temp_max_pre72"] = d["stay_id"].map(vals)
    _median_impute_with_flag(d, "temp_max_pre72")
    return {
        "n_events_used": int(len(t)),
        "n_stays_observed": int(vals.index.nunique()),
        "unit": "degrees_C_after_reading_level_normalization",
    }


def _reconstruct_gcs(d: pd.DataFrame, events: pd.DataFrame) -> dict:
    # A GCS total is accepted only when eye, verbal, and motor components occur
    # at the same stay/timestamp. We do not sum components from different times.
    g = events.loc[events["domain"] == "gcs"].copy()
    ranges = {"eye": (1, 4), "verbal": (1, 5), "motor": (1, 6)}
    valid_parts = []
    for component, (lo, hi) in ranges.items():
        s = g.loc[g["gcs_component"] == component].copy()
        s["valuenum"] = pd.to_numeric(s["valuenum"], errors="coerce")
        valid_parts.append(s.loc[s["valuenum"].between(lo, hi, inclusive="both")])
    g = pd.concat(valid_parts, ignore_index=True) if valid_parts else g.iloc[0:0].copy()

    p = (
        g.groupby(["stay_id", "charttime", "gcs_component"], as_index=False)["valuenum"]
        .max()
        .pivot(index=["stay_id", "charttime"], columns="gcs_component", values="valuenum")
        .reset_index()
    )
    for component in ["eye", "verbal", "motor"]:
        if component not in p:
            p[component] = np.nan
    p = p.dropna(subset=["eye", "verbal", "motor"]).copy()
    p["gcs_total"] = p["eye"] + p["verbal"] + p["motor"]
    p = p.loc[p["gcs_total"].between(3, 15, inclusive="both")].copy()
    p = p.merge(
        d[["stay_id", "first_broad_time", "decision_time"]].drop_duplicates("stay_id"),
        on="stay_id",
        how="inner",
    )

    early = p.loc[
        (p["charttime"] >= p["first_broad_time"])
        & (p["charttime"] <= p["first_broad_time"] + pd.Timedelta(hours=24))
    ]
    late = p.loc[
        (p["charttime"] >= p["decision_time"] - pd.Timedelta(hours=24))
        & (p["charttime"] <= p["decision_time"])
    ]
    early_vals = early.groupby("stay_id")["gcs_total"].min()
    late_vals = late.groupby("stay_id")["gcs_total"].min()
    d["gcs_total_0_24h"] = d["stay_id"].map(early_vals)
    d["gcs_total_48_72h"] = d["stay_id"].map(late_vals)
    _median_impute_with_flag(d, "gcs_total_0_24h")
    _median_impute_with_flag(d, "gcs_total_48_72h")
    return {
        "n_complete_triplets": int(len(p)),
        "n_stays_early": int(early_vals.index.nunique()),
        "n_stays_late": int(late_vals.index.nunique()),
    }


def _correct_fio2(d: pd.DataFrame, events: pd.DataFrame) -> dict:
    # MIMIC FiO2 appears both as fraction and percentage. Normalize to fraction
    # before aggregation. The direct FiO2 term remains excluded from the primary
    # PS, but keeping the corrected feature makes the analytic cohort reproducible.
    f = events.loc[events["domain"] == "fio2"].copy()
    v = pd.to_numeric(f["valuenum"], errors="coerce")
    f["fio2_fraction"] = np.where(v > 1.5, v / 100.0, v)
    f = f.loc[f["fio2_fraction"].between(0.20, 1.0, inclusive="both")].copy()
    early = f.loc[
        (f["charttime"] >= f["first_broad_time"])
        & (f["charttime"] <= f["first_broad_time"] + pd.Timedelta(hours=24))
    ]
    late = f.loc[
        (f["charttime"] >= f["decision_time"] - pd.Timedelta(hours=24))
        & (f["charttime"] <= f["decision_time"])
    ]
    early_vals = early.groupby("stay_id")["fio2_fraction"].max()
    late_vals = late.groupby("stay_id")["fio2_fraction"].max()
    d["fio2_0_24h"] = d["stay_id"].map(early_vals)
    d["fio2_48_72h"] = d["stay_id"].map(late_vals)
    _median_impute_with_flag(d, "fio2_0_24h")
    _median_impute_with_flag(d, "fio2_48_72h")
    return {
        "n_events_used": int(len(f)),
        "n_stays_early": int(early_vals.index.nunique()),
        "n_stays_late": int(late_vals.index.nunique()),
        "unit": "fraction_0_to_1",
    }


def apply_audited_vital_corrections(
    cohort: pd.DataFrame,
    source,
    chunksize: int = 500_000,
) -> tuple[pd.DataFrame, dict]:
    """Apply the final audited vital corrections before SOFA-like calculation/PS fitting."""
    d = cohort.copy()
    selected = _select_audited_items(source)
    events = _load_events(source, d, selected, chunksize)
    summary = {
        "temperature": _correct_baseline_temperature(d, events),
        "gcs": _reconstruct_gcs(d, events),
        "fio2": _correct_fio2(d, events),
        "guardrail": (
            "Audit-driven measurement correction only. Treatment, outcomes, cohort eligibility, "
            "and propensity estimand are unchanged. Direct GCS/FiO2 remain excluded from the primary PS."
        ),
    }
    return d, summary

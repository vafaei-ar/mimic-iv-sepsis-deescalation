#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from sepsis_deescalation.config import resolve_mimic_source, load_config
from sepsis_deescalation.features import add_sofa_like
from sepsis_deescalation.mimic_io import read_csv, read_csv_filtered
from sepsis_deescalation.specification import CANDIDATE_PS_VARS
from sepsis_deescalation.stats import balance_table, fit_stabilized_iptw, risks


def _temperature_c(values: pd.Series, labels: pd.Series) -> pd.Series:
    v = pd.to_numeric(values, errors="coerce")
    lab = labels.fillna("").astype(str).str.lower()
    is_f = lab.str.contains("fahrenheit|temperaturef", regex=True, na=False)
    is_c = lab.str.contains("celsius|\(c\)", regex=True, na=False)
    generic_f = (~is_f & ~is_c) & (v > 60)
    out = v.copy()
    out.loc[is_f | generic_f] = (v.loc[is_f | generic_f] - 32.0) * 5.0 / 9.0
    return out


def _median_impute_with_flag(d: pd.DataFrame, col: str) -> None:
    x = pd.to_numeric(d[col], errors="coerce")
    d[f"{col}_missing"] = x.isna().astype(int)
    med = x.median()
    d[col] = x.fillna(0.0 if pd.isna(med) else med)


def _select_items(source: Path) -> tuple[pd.DataFrame, dict[str, set[int]]]:
    items = read_csv(source, "icu/d_items.csv.gz", usecols=["itemid", "label", "unitname"])
    items["label_lower"] = items["label"].fillna("").astype(str).str.lower().str.strip()

    temp = items.loc[
        items["label_lower"].isin({
            "temperature fahrenheit",
            "temperature celsius",
            "blood temperature cco (c)",
            "cerebral temperature (c)",
        })
    ].copy()

    fio2 = items.loc[
        items["label_lower"].str.fullmatch(r"inspired o2 fraction|fraction inspired oxygen", na=False)
    ].copy()

    component_map = {
        "gcs - eye opening": "eye",
        "gcs - verbal response": "verbal",
        "gcs - motor response": "motor",
    }
    gcs = items.loc[items["label_lower"].isin(component_map)].copy()
    gcs["gcs_component"] = gcs["label_lower"].map(component_map)

    selected = pd.concat(
        [temp.assign(domain="temperature"), fio2.assign(domain="fio2"), gcs.assign(domain="gcs")],
        ignore_index=True,
        sort=False,
    )
    ids = {
        "temperature": set(temp["itemid"].astype(int)),
        "fio2": set(fio2["itemid"].astype(int)),
        "gcs": set(gcs["itemid"].astype(int)),
    }
    if not ids["temperature"]:
        raise RuntimeError("No true temperature items were found in d_items.")
    if not ids["gcs"]:
        raise RuntimeError("No routine GCS component items were found in d_items.")
    if not ids["fio2"]:
        raise RuntimeError(
            "No routine bedside FiO2 item matched 'Inspired O2 Fraction'. Inspect d_items before proceeding."
        )
    return selected, ids


def _load_events(source: Path, cohort: pd.DataFrame, selected: pd.DataFrame, chunksize: int) -> pd.DataFrame:
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
    events = events.loc[events["valuenum"].notna() & events["charttime"].notna()].copy()
    return events


def _repair_temperature(d: pd.DataFrame, events: pd.DataFrame) -> dict:
    t = events.loc[events["domain"] == "temperature"].copy()
    t = t.loc[(t["charttime"] >= t["intime"]) & (t["charttime"] <= t["decision_time"])].copy()
    t["temperature_c"] = _temperature_c(t["valuenum"], t["label"])
    t = t.loc[t["temperature_c"].between(30.0, 45.0, inclusive="both")].copy()
    vals = t.groupby("stay_id")["temperature_c"].max()
    old = pd.to_numeric(d.get("temp_max_pre72"), errors="coerce").copy()
    d["temp_max_pre72"] = d["stay_id"].map(vals)
    _median_impute_with_flag(d, "temp_max_pre72")
    old_c = np.where(old > 60, (old - 32.0) * 5.0 / 9.0, old)
    diff = np.abs(pd.to_numeric(d["temp_max_pre72"]) - old_c)
    return {
        "n_events_used": int(len(t)),
        "n_stays_observed": int(vals.index.nunique()),
        "n_stays_changed_gt_0_1c": int(np.nansum(diff > 0.1)),
        "n_stays_changed_gt_0_5c": int(np.nansum(diff > 0.5)),
        "max_abs_change_c": float(np.nanmax(diff)) if np.isfinite(diff).any() else None,
    }


def _repair_gcs(d: pd.DataFrame, events: pd.DataFrame) -> dict:
    g = events.loc[events["domain"] == "gcs"].copy()
    ranges = {"eye": (1, 4), "verbal": (1, 5), "motor": (1, 6)}
    valid_parts = []
    for comp, (lo, hi) in ranges.items():
        s = g.loc[g["gcs_component"] == comp].copy()
        s["valuenum"] = pd.to_numeric(s["valuenum"], errors="coerce")
        s = s.loc[s["valuenum"].between(lo, hi, inclusive="both")]
        valid_parts.append(s)
    g = pd.concat(valid_parts, ignore_index=True) if valid_parts else g.iloc[0:0].copy()

    p = (
        g.groupby(["stay_id", "charttime", "gcs_component"], as_index=False)["valuenum"]
        .max()
        .pivot(index=["stay_id", "charttime"], columns="gcs_component", values="valuenum")
        .reset_index()
    )
    for comp in ["eye", "verbal", "motor"]:
        if comp not in p:
            p[comp] = np.nan
    p = p.dropna(subset=["eye", "verbal", "motor"]).copy()
    p["gcs_total"] = p["eye"] + p["verbal"] + p["motor"]
    p = p.loc[p["gcs_total"].between(3, 15, inclusive="both")].copy()
    p = p.merge(d[["stay_id", "first_broad_time", "decision_time"]].drop_duplicates("stay_id"), on="stay_id", how="inner")

    early = p.loc[
        (p["charttime"] >= p["first_broad_time"]) &
        (p["charttime"] <= p["first_broad_time"] + pd.Timedelta(hours=24))
    ]
    late = p.loc[
        (p["charttime"] >= p["decision_time"] - pd.Timedelta(hours=24)) &
        (p["charttime"] <= p["decision_time"])
    ]
    early_vals = early.groupby("stay_id")["gcs_total"].min()
    late_vals = late.groupby("stay_id")["gcs_total"].min()
    d["gcs_total_0_24h"] = d["stay_id"].map(early_vals)
    d["gcs_total_48_72h"] = d["stay_id"].map(late_vals)
    _median_impute_with_flag(d, "gcs_total_0_24h")
    _median_impute_with_flag(d, "gcs_total_48_72h")
    return {
        "n_complete_triplets": int(len(p)),
        "n_stays_any_complete_triplet": int(p["stay_id"].nunique()),
        "n_stays_early": int(early_vals.index.nunique()),
        "n_stays_late": int(late_vals.index.nunique()),
    }


def _repair_fio2(d: pd.DataFrame, events: pd.DataFrame) -> dict:
    f = events.loc[events["domain"] == "fio2"].copy()
    v = pd.to_numeric(f["valuenum"], errors="coerce")
    f["fio2_fraction"] = np.where(v > 1.5, v / 100.0, v)
    f = f.loc[f["fio2_fraction"].between(0.20, 1.0, inclusive="both")].copy()
    early = f.loc[
        (f["charttime"] >= f["first_broad_time"]) &
        (f["charttime"] <= f["first_broad_time"] + pd.Timedelta(hours=24))
    ]
    late = f.loc[
        (f["charttime"] >= f["decision_time"] - pd.Timedelta(hours=24)) &
        (f["charttime"] <= f["decision_time"])
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
        "median_imputed_early": float(pd.to_numeric(d["fio2_0_24h"]).median()),
        "median_imputed_late": float(pd.to_numeric(d["fio2_48_72h"]).median()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Repair only the audited MIMIC v5.7 vital covariates, without rebuilding cohort membership/exposure/outcomes."
    )
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--config", default="config/mimic.yaml")
    parser.add_argument("--mimic-source", type=Path, default=None)
    parser.add_argument("--chunksize", type=int, default=500_000)
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    cohort_path = run_dir / "analysis_cohort_weighted.csv"
    if not cohort_path.exists():
        raise FileNotFoundError(cohort_path)
    d = pd.read_csv(cohort_path, low_memory=False)
    original = d.copy()
    for col in ["intime", "first_broad_time", "decision_time"]:
        d[col] = pd.to_datetime(d[col], errors="coerce")

    cfg = load_config(args.config)
    source = args.mimic_source.expanduser().resolve() if args.mimic_source else resolve_mimic_source(cfg)
    selected, _ = _select_items(source)
    events = _load_events(source, d, selected, args.chunksize)

    out = run_dir / "audits" / "vital_repair"
    out.mkdir(parents=True, exist_ok=True)
    selected.to_csv(out / "selected_repair_items.csv", index=False)

    summary = {
        "temperature": _repair_temperature(d, events),
        "gcs": _repair_gcs(d, events),
        "fio2": _repair_fio2(d, events),
    }

    d = add_sofa_like(d)

    old_w, _, _ = fit_stabilized_iptw(original, CANDIDATE_PS_VARS)
    new_w, _, _ = fit_stabilized_iptw(d, CANDIDATE_PS_VARS)
    old_rt, old_rc, old_rd, old_rr = risks(old_w, "death_by_horizon", "SW_A")
    new_rt, new_rc, new_rd, new_rr = risks(new_w, "death_by_horizon", "SW_A")
    old_bal = balance_table(old_w, CANDIDATE_PS_VARS)
    new_bal = balance_table(new_w, CANDIDATE_PS_VARS)

    comparison = pd.DataFrame([
        {
            "analysis": "original_v57",
            "risk_deescalated_stopped": old_rt,
            "risk_continued": old_rc,
            "risk_difference": old_rd,
            "risk_ratio": old_rr,
            "max_post_smd": float(old_bal["after"].max()) if len(old_bal) else np.nan,
        },
        {
            "analysis": "vital_corrected",
            "risk_deescalated_stopped": new_rt,
            "risk_continued": new_rc,
            "risk_difference": new_rd,
            "risk_ratio": new_rr,
            "max_post_smd": float(new_bal["after"].max()) if len(new_bal) else np.nan,
        },
    ])
    comparison.to_csv(out / "point_estimate_comparison.csv", index=False)
    new_bal.to_csv(out / "balance_vital_corrected.csv", index=False)

    # Save the repaired cohort with the repaired propensity-score weights, not
    # the stale v5.7 SW_A/ps columns inherited from the original CSV.
    for col in ["ps_den", "ps_num", "SW_A"]:
        d[col] = new_w[col].to_numpy()
    d.to_csv(out / "analysis_cohort_vital_corrected.csv", index=False)
    (out / "vital_repair_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print(comparison.to_string(index=False))
    print(f"Corrected cohort: {out / 'analysis_cohort_vital_corrected.csv'}")


if __name__ == "__main__":
    main()

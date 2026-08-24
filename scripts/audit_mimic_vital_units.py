#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _read_csv(source: Path, rel: str, **kwargs) -> pd.DataFrame:
    path = source / rel
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, **kwargs)


def _source_from_manifest(run_dir: Path) -> Path | None:
    manifest = run_dir / "run_manifest.json"
    if not manifest.exists():
        return None
    obj = json.loads(manifest.read_text(encoding="utf-8"))
    value = obj.get("source")
    if not value:
        return None
    path = Path(value).expanduser()
    return path if path.exists() else None


def _summarize_items(events: pd.DataFrame, label_col: str = "label") -> pd.DataFrame:
    if events.empty:
        return pd.DataFrame(columns=["itemid", "label", "n", "n_stays", "min", "p01", "p05", "median", "p95", "p99", "max"])
    rows = []
    for (itemid, label), sub in events.groupby(["itemid", label_col], dropna=False):
        v = pd.to_numeric(sub["valuenum"], errors="coerce").dropna()
        if v.empty:
            continue
        rows.append({
            "itemid": itemid,
            "label": label,
            "n": int(len(v)),
            "n_stays": int(sub.loc[v.index, "stay_id"].nunique()),
            "min": float(v.min()),
            "p01": float(v.quantile(0.01)),
            "p05": float(v.quantile(0.05)),
            "median": float(v.median()),
            "p95": float(v.quantile(0.95)),
            "p99": float(v.quantile(0.99)),
            "max": float(v.max()),
        })
    return pd.DataFrame(rows).sort_values(["n", "label"], ascending=[False, True])


def _temperature_c(value: pd.Series, label: pd.Series) -> pd.Series:
    v = pd.to_numeric(value, errors="coerce")
    lab = label.fillna("").astype(str).str.lower()
    is_f = lab.str.contains("fahrenheit|deg f|°f", regex=True, na=False)
    is_c = lab.str.contains("celsius|deg c|°c", regex=True, na=False)
    # For generic temperature labels, values >60 are Fahrenheit in MIMIC; lower
    # values are Celsius. Explicit label information takes precedence.
    generic_f = (~is_f & ~is_c) & (v > 60)
    out = v.copy()
    out.loc[is_f | generic_f] = (v.loc[is_f | generic_f] - 32.0) * 5.0 / 9.0
    return out


def _fio2_fraction(value: pd.Series) -> pd.Series:
    v = pd.to_numeric(value, errors="coerce")
    out = v.copy()
    pct = (v > 1.5) & (v <= 100)
    out.loc[pct] = v.loc[pct] / 100.0
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit MIMIC temperature, FiO2, and GCS item/value semantics before freezing the v5.7 analysis."
    )
    parser.add_argument("run_dir", type=Path, help="Completed MIMIC run directory")
    parser.add_argument("--mimic-source", type=Path, default=None, help="MIMIC-IV root; defaults to run_manifest.json source")
    parser.add_argument("--chunksize", type=int, default=500_000)
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    cohort_path = run_dir / "analysis_cohort_weighted.csv"
    if not cohort_path.exists():
        raise FileNotFoundError(cohort_path)
    cohort = pd.read_csv(cohort_path, low_memory=False)
    for col in ["stay_id", "intime", "decision_time"]:
        if col not in cohort:
            raise ValueError(f"analysis cohort missing required column: {col}")
    cohort["intime"] = pd.to_datetime(cohort["intime"], errors="coerce")
    cohort["decision_time"] = pd.to_datetime(cohort["decision_time"], errors="coerce")

    source = args.mimic_source.expanduser().resolve() if args.mimic_source else _source_from_manifest(run_dir)
    if source is None or not source.exists():
        raise FileNotFoundError("Could not resolve MIMIC source; pass --mimic-source")

    out = run_dir / "audits" / "vital_unit_review"
    out.mkdir(parents=True, exist_ok=True)

    items = _read_csv(source, "icu/d_items.csv.gz", usecols=lambda c: c in {"itemid", "label", "unitname"})
    items["label_lower"] = items["label"].fillna("").astype(str).str.lower()
    temp_items = items.loc[items["label_lower"].str.contains(r"temperature", regex=True, na=False)].copy()
    fio2_items = items.loc[items["label_lower"].str.contains(r"fio2|fraction.*inspired.*oxygen", regex=True, na=False)].copy()
    gcs_items = items.loc[items["label_lower"].str.contains(r"gcs|glasgow", regex=True, na=False)].copy()

    selected = pd.concat([
        temp_items.assign(domain="temperature"),
        fio2_items.assign(domain="fio2"),
        gcs_items.assign(domain="gcs"),
    ], ignore_index=True).drop_duplicates(["itemid", "domain"])
    selected.to_csv(out / "selected_d_items.csv", index=False)

    item_to_domain = selected.groupby("itemid")["domain"].apply(list).to_dict()
    wanted_itemids = set(int(x) for x in selected["itemid"].dropna().astype(int))
    wanted_stays = set(int(x) for x in cohort["stay_id"].dropna().astype(int))

    event_parts = []
    chartevents_path = source / "icu/chartevents.csv.gz"
    usecols = ["stay_id", "charttime", "itemid", "valuenum"]
    for chunk in pd.read_csv(chartevents_path, usecols=usecols, parse_dates=["charttime"], chunksize=args.chunksize):
        hit = chunk.loc[chunk["stay_id"].isin(wanted_stays) & chunk["itemid"].isin(wanted_itemids) & chunk["valuenum"].notna()].copy()
        if not hit.empty:
            event_parts.append(hit)
    events = pd.concat(event_parts, ignore_index=True) if event_parts else pd.DataFrame(columns=usecols)
    events = events.merge(selected[["itemid", "label", "unitname", "domain"]], on="itemid", how="inner")
    windows = cohort[["stay_id", "intime", "decision_time"]].drop_duplicates("stay_id")
    events = events.merge(windows, on="stay_id", how="inner")
    events = events.loc[(events["charttime"] >= events["intime"]) & (events["charttime"] <= events["decision_time"])].copy()

    for domain in ["temperature", "fio2", "gcs"]:
        sub = events.loc[events["domain"] == domain].copy()
        _summarize_items(sub).to_csv(out / f"{domain}_item_summary.csv", index=False)

    # Temperature unit audit and corrected patient-level aggregate comparison.
    temp = events.loc[events["domain"] == "temperature"].copy()
    temp["temperature_c"] = _temperature_c(temp["valuenum"], temp["label"])
    temp["plausible_c"] = temp["temperature_c"].between(30.0, 45.0, inclusive="both")
    temp_plaus = temp.loc[temp["plausible_c"]].copy()
    corrected_max_c = temp_plaus.groupby("stay_id")["temperature_c"].max()

    comp = cohort[["stay_id", "temp_max_pre72"]].drop_duplicates("stay_id").copy()
    comp["current_temp_max_raw"] = pd.to_numeric(comp["temp_max_pre72"], errors="coerce")
    comp["current_temp_max_as_c"] = np.where(
        comp["current_temp_max_raw"] > 60,
        (comp["current_temp_max_raw"] - 32.0) * 5.0 / 9.0,
        comp["current_temp_max_raw"],
    )
    comp["corrected_temp_max_c"] = comp["stay_id"].map(corrected_max_c)
    comp["abs_difference_c"] = (comp["corrected_temp_max_c"] - comp["current_temp_max_as_c"]).abs()
    # Patient-level comparison is intentionally local and remains under outputs/.
    comp.to_csv(out / "temperature_patient_comparison.csv", index=False)

    # FiO2 normalization audit. This does not alter the production feature yet.
    fio2 = events.loc[events["domain"] == "fio2"].copy()
    fio2["fio2_fraction"] = _fio2_fraction(fio2["valuenum"])
    fio2["plausible_fraction"] = fio2["fio2_fraction"].between(0.20, 1.0, inclusive="both")
    fio2.groupby(["itemid", "label"], dropna=False).agg(
        n=("valuenum", "size"),
        n_stays=("stay_id", "nunique"),
        raw_min=("valuenum", "min"),
        raw_median=("valuenum", "median"),
        raw_max=("valuenum", "max"),
        normalized_min=("fio2_fraction", "min"),
        normalized_median=("fio2_fraction", "median"),
        normalized_max=("fio2_fraction", "max"),
        plausible_fraction=("plausible_fraction", "mean"),
    ).reset_index().sort_values("n", ascending=False).to_csv(out / "fio2_normalization_summary.csv", index=False)

    # GCS coverage audit. Total-score and component labels are kept separate so
    # we can decide whether the current feature should use a total score, a sum
    # of components, or be removed from the harmonized model.
    gcs = events.loc[events["domain"] == "gcs"].copy()
    gcs["is_total_label"] = gcs["label"].fillna("").str.lower().str.contains(r"total|gcs total|glasgow coma scale total", regex=True, na=False)
    gcs["is_component_label"] = gcs["label"].fillna("").str.lower().str.contains(r"eye|verbal|motor", regex=True, na=False)
    gcs.groupby(["itemid", "label", "is_total_label", "is_component_label"], dropna=False).agg(
        n=("valuenum", "size"),
        n_stays=("stay_id", "nunique"),
        min=("valuenum", "min"),
        median=("valuenum", "median"),
        max=("valuenum", "max"),
    ).reset_index().sort_values("n", ascending=False).to_csv(out / "gcs_label_coverage.csv", index=False)

    summary = {
        "n_cohort": int(len(cohort)),
        "n_temperature_events": int(len(temp)),
        "n_temperature_events_plausible_after_normalization": int(len(temp_plaus)),
        "n_stays_with_corrected_temperature": int(corrected_max_c.index.nunique()),
        "n_current_temp_max_below_60": int((pd.to_numeric(cohort["temp_max_pre72"], errors="coerce") < 60).sum()),
        "n_stays_temperature_difference_gt_0_1c": int((comp["abs_difference_c"] > 0.1).sum()),
        "n_stays_temperature_difference_gt_0_5c": int((comp["abs_difference_c"] > 0.5).sum()),
        "max_temperature_difference_c": float(comp["abs_difference_c"].max()) if comp["abs_difference_c"].notna().any() else None,
        "n_fio2_events": int(len(fio2)),
        "n_fio2_stays": int(fio2["stay_id"].nunique()),
        "n_gcs_events": int(len(gcs)),
        "n_gcs_stays": int(gcs["stay_id"].nunique()),
        "n_gcs_total_events": int(gcs["is_total_label"].sum()),
        "n_gcs_component_events": int(gcs["is_component_label"].sum()),
    }
    (out / "vital_unit_review_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()

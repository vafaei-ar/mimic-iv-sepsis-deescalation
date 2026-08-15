from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .config import load_config, make_run_paths
from .harmonization import default_crosswalk
from .pcornet_antibiotics import add_antibiotic_classification, antibiotic_mapping_audit
from .pcornet_io import combine_date_time, normalize_columns, read_table, resolve_path
from .provenance import copy_config, write_run_manifest, zip_run

SCREEN_PATTERN = r"screen|surveillance|mrsa|vre|rectal|nares|nasal|stool ova|parasite|viral screen"


def _write(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def _rename_from_mapping(df: pd.DataFrame, mapping: dict[str, str], keys: list[str]) -> pd.DataFrame:
    d = normalize_columns(df)
    rename = {}
    for canonical in keys:
        source = mapping.get(canonical)
        if source and source.lower() in d.columns:
            rename[source.lower()] = canonical
    return d.rename(columns=rename)


def _medication_table(cfg: dict, source_name: str) -> pd.DataFrame:
    meds = cfg["medications"]
    if source_name == "prescribing":
        mapping = meds["prescribing_columns"]
        raw = read_table(resolve_path(cfg, "prescribing"))
        keys = ["patid", "encounterid", "code", "raw_code", "raw_name", "route", "start_date", "start_time", "end_date"]
        d = _rename_from_mapping(raw, mapping, keys)
        if "code" not in d:
            d["code"] = d.get("raw_code", "")
        if "raw_name" not in d:
            raise ValueError("PRESCRIBING requires a configured raw medication name column.")
        d["coverage_start"] = combine_date_time(d["start_date"], d.get("start_time"))
        d["coverage_stop"] = pd.to_datetime(d.get("end_date"), errors="coerce")
    elif source_name == "med_admin":
        mapping = meds["med_admin_columns"]
        raw = read_table(resolve_path(cfg, "med_admin"))
        keys = ["patid", "encounterid", "code", "code_type", "raw_name", "route", "raw_route", "start_date", "start_time", "stop_date", "stop_time"]
        d = _rename_from_mapping(raw, mapping, keys)
        if "raw_name" not in d:
            raise ValueError("MED_ADMIN requires a configured raw medication name column.")
        d["coverage_start"] = combine_date_time(d["start_date"], d.get("start_time"))
        d["coverage_stop"] = combine_date_time(d["stop_date"], d.get("stop_time")) if "stop_date" in d else pd.NaT
        if "route" not in d and "raw_route" in d:
            d["route"] = d["raw_route"]
    else:
        raise ValueError(f"Unknown medication source: {source_name}")

    for col in ["patid", "encounterid", "code", "raw_name", "route"]:
        if col not in d:
            d[col] = ""
    missing_stop = d["coverage_stop"].isna()
    fill_hours = int(meds.get("missing_stop_fill_hours", 24))
    d.loc[missing_stop, "coverage_stop"] = d.loc[missing_stop, "coverage_start"] + pd.Timedelta(hours=fill_hours)
    bad = d["coverage_stop"] <= d["coverage_start"]
    d.loc[bad, "coverage_stop"] = d.loc[bad, "coverage_start"] + pd.Timedelta(hours=1)
    d["missing_stop_assumed"] = missing_stop.astype(int)
    return add_antibiotic_classification(d, "code", "raw_name", "route")


def _icu_table(cfg: dict) -> pd.DataFrame:
    raw = read_table(resolve_path(cfg, "icu_stays"))
    d = _rename_from_mapping(raw, cfg["icu"], ["patid", "encounterid", "start_datetime", "end_datetime", "careunit"])
    required = ["patid", "encounterid", "start_datetime", "end_datetime"]
    missing = [c for c in required if c not in d]
    if missing:
        raise ValueError(f"Local ICU table is missing mapped columns: {missing}")
    d["icu_start"] = pd.to_datetime(d["start_datetime"], errors="coerce")
    d["icu_end"] = pd.to_datetime(d["end_datetime"], errors="coerce")
    d = d.loc[d["icu_start"].notna()].sort_values(["encounterid", "icu_start"]).drop_duplicates("encounterid", keep="first")
    return d


def _base_encounters(cfg: dict, icu: pd.DataFrame) -> pd.DataFrame:
    enc = normalize_columns(read_table(resolve_path(cfg, "encounter")))
    dem = normalize_columns(read_table(resolve_path(cfg, "demographic")))
    required_enc = ["patid", "encounterid", "admit_date", "discharge_date", "enc_type"]
    required_dem = ["patid", "birth_date", "sex", "race"]
    for name, d, req in [("ENCOUNTER", enc, required_enc), ("DEMOGRAPHIC", dem, required_dem)]:
        missing = [c for c in req if c not in d]
        if missing:
            raise ValueError(f"{name} missing columns: {missing}")
    e = enc.loc[enc["enc_type"].isin(["IP", "EI"])].copy()
    e["admit_dt"] = combine_date_time(e["admit_date"], e["admit_time"] if "admit_time" in e else None)
    if "discharge_time" in e:
        e["discharge_dt"] = combine_date_time(e["discharge_date"], e["discharge_time"])
    else:
        e["discharge_dt"] = pd.to_datetime(e["discharge_date"], errors="coerce") + pd.Timedelta(hours=23, minutes=59)
        e["discharge_time_imputed_to_end_of_day"] = 1
    e = e.merge(dem[[c for c in ["patid", "birth_date", "sex", "race", "hispanic"] if c in dem]], on="patid", how="left")
    e = e.merge(icu[[c for c in ["patid", "encounterid", "icu_start", "icu_end", "careunit"] if c in icu]], on=["patid", "encounterid"], how="inner")
    e["birth_date"] = pd.to_datetime(e["birth_date"], errors="coerce")
    e["age"] = (e["icu_start"] - e["birth_date"]).dt.total_seconds() / (365.25 * 86400)
    return e.loc[e["age"] >= 18].copy()


def _overlap(events: pd.DataFrame, cohort: pd.DataFrame, start_col: str, end_col: str) -> pd.DataFrame:
    w = cohort[["encounterid", start_col, end_col]].rename(columns={start_col: "window_start", end_col: "window_end"})
    x = events.merge(w, on="encounterid", how="inner")
    return x.loc[(x["coverage_start"] < x["window_end"]) & (x["coverage_stop"] > x["window_start"])].copy()


def _microbiology(cfg: dict, cohort: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    micro_cfg = cfg["microbiology"]
    source = micro_cfg.get("source", "local")
    if source == "local":
        raw = read_table(resolve_path(cfg, "microbiology_local"))
        d = _rename_from_mapping(raw, micro_cfg["local_columns"], ["patid", "encounterid", "specimen_datetime", "result_datetime", "specimen", "test_name", "organism"])
        d["specimen_time"] = pd.to_datetime(d["specimen_datetime"], errors="coerce")
        d["result_available_time"] = pd.to_datetime(d["result_datetime"], errors="coerce")
    else:
        raw = read_table(resolve_path(cfg, "lab_result_cm"))
        d = _rename_from_mapping(raw, micro_cfg["lab_result_cm_columns"], ["patid", "encounterid", "specimen_date", "specimen_time", "result_date", "result_time", "specimen", "test_name", "loinc", "result_qual", "raw_result", "result_snomed"])
        d["specimen_time"] = combine_date_time(d["specimen_date"], d.get("specimen_time"))
        d["result_available_time"] = combine_date_time(d["result_date"], d.get("result_time"))
        # This fallback is intentionally conservative and must be audited before the external analysis is frozen.
        d["organism"] = d.get("raw_result", "").fillna("")
    for c in ["specimen", "test_name", "organism"]:
        if c not in d:
            d[c] = ""
        d[c] = d[c].fillna("").astype(str)
    d["combined"] = d["specimen"].str.lower() + " " + d["test_name"].str.lower()
    d["clinical_micro"] = (~d["combined"].str.contains(SCREEN_PATTERN, na=False, regex=True)).astype(int)
    d["positive_organism"] = ((d["clinical_micro"] == 1) & d["organism"].str.strip().ne("")).astype(int)
    d = d.merge(cohort[["encounterid", "first_broad_time", "decision_time"]], on="encounterid", how="inner")
    win = d.loc[d["specimen_time"].notna() & (d["specimen_time"] >= d["first_broad_time"] - pd.Timedelta(hours=24)) & (d["specimen_time"] <= d["decision_time"]) & (d["clinical_micro"] == 1)].copy()
    sampled = set(win["encounterid"])
    positive = set(win.loc[(win["positive_organism"] == 1) & win["result_available_time"].notna() & (win["result_available_time"] <= win["decision_time"]), "encounterid"])
    eligible = cohort.loc[cohort["encounterid"].isin(sampled - positive)].copy()
    return eligible, win


def _vasopressor_filter(cfg: dict, cohort: pd.DataFrame) -> pd.DataFrame:
    ma = _medication_table(cfg, "med_admin")
    patterns = "|".join(cfg.get("vasopressors", {}).get("name_patterns", []))
    if not patterns:
        raise ValueError("No vasopressor name patterns configured.")
    vaso = ma.loc[ma["raw_name"].str.lower().str.contains(patterns, na=False, regex=True)].copy()
    w = cohort[["encounterid", "decision_time"]].copy()
    w["window_start"] = w["decision_time"] - pd.Timedelta(hours=int(cfg["clock"]["vasopressor_free_lookback_hours"]))
    w["window_end"] = w["decision_time"]
    x = vaso.merge(w[["encounterid", "window_start", "window_end"]], on="encounterid", how="inner")
    active = set(x.loc[(x["coverage_start"] < x["window_end"]) & (x["coverage_stop"] > x["window_start"]), "encounterid"])
    return cohort.loc[~cohort["encounterid"].isin(active)].copy()


def _classify(cohort: pd.DataFrame, meds: pd.DataFrame) -> pd.DataFrame:
    d = cohort.copy()
    broad = meds.loc[meds["is_broad_spectrum"] == 1].copy()
    systemic = meds.loc[meds["is_systemic_antibiotic"] == 1].copy()
    d["window_48_start"] = d["first_broad_time"] + pd.Timedelta(hours=48)
    d["window_72"] = d["decision_time"]
    pre = _overlap(broad, d, "window_48_start", "window_72")
    d = d.loc[d["encounterid"].isin(set(pre["encounterid"]))].copy()
    d["window_96"] = d["analysis_time0"]
    post_b = _overlap(broad, d, "decision_time", "window_96")
    post_a = _overlap(systemic, d, "decision_time", "window_96")
    continued = set(post_b["encounterid"]); any_abx = set(post_a["encounterid"])
    d["A"] = (~d["encounterid"].isin(continued)).astype(int)
    d["deescalation_type"] = "continued_broad"
    d.loc[(d["A"] == 1) & d["encounterid"].isin(any_abx), "deescalation_type"] = "narrowed_or_non_broad_only"
    d.loc[(d["A"] == 1) & ~d["encounterid"].isin(any_abx), "deescalation_type"] = "stopped_all_observed_systemic_antibiotics"
    return d


def run_pcornet(config_path: str | Path) -> Path:
    cfg = load_config(config_path)
    paths = make_run_paths(cfg["output_root"], cfg["site"], cfg["analysis_version"])
    copy_config(cfg, paths.run_dir)
    write_run_manifest(paths.run_dir, cfg, source="PSU/PCORnet local datamart")
    _write(default_crosswalk(), paths.tables / "harmonization_crosswalk.csv")

    icu = _icu_table(cfg)
    base = _base_encounters(cfg, icu)
    primary_source = cfg["medications"].get("primary_source", "prescribing")
    meds = _medication_table(cfg, primary_source)
    _write(antibiotic_mapping_audit(meds, "code", "raw_name", "route"), paths.audits / f"{primary_source}_antibiotic_mapping.csv")

    broad = meds.loc[meds["is_broad_spectrum"] == 1].merge(base[["encounterid", "icu_start"]], on="encounterid", how="inner")
    before = int(cfg["clock"].get("qualifying_start_hours_before_icu", 6)); after = int(cfg["clock"].get("qualifying_start_hours_after_icu", 24))
    broad = broad.loc[(broad["coverage_start"] >= broad["icu_start"] - pd.Timedelta(hours=before)) & (broad["coverage_start"] <= broad["icu_start"] + pd.Timedelta(hours=after))]
    first = broad.sort_values(["encounterid", "coverage_start"]).drop_duplicates("encounterid")[["encounterid", "coverage_start", "broad_drug_group"]].rename(columns={"coverage_start": "first_broad_time", "broad_drug_group": "first_broad_drug"})
    cohort = base.merge(first, on="encounterid", how="inner")
    cohort["decision_time"] = cohort["first_broad_time"] + pd.Timedelta(hours=int(cfg["clock"]["decision_hours"]))
    cohort["analysis_time0"] = cohort["first_broad_time"] + pd.Timedelta(hours=int(cfg["clock"]["classification_end_hours"]))
    cohort = cohort.loc[cohort["discharge_dt"].notna() & (cohort["discharge_dt"] > cohort["analysis_time0"])].copy()
    n_before_micro = len(cohort)
    cohort, micro_win = _microbiology(cfg, cohort)
    n_after_micro = len(cohort)
    cohort = _vasopressor_filter(cfg, cohort)
    n_after_vaso = len(cohort)
    cohort = _classify(cohort, meds)

    flow = pd.DataFrame([
        {"step": "Adult IP/EI encounters with timestamped ICU episode", "n": len(base)},
        {"step": "Qualifying early broad-spectrum exposure; hospitalized through 96 h", "n": n_before_micro},
        {"step": "Clinical microbiology sampled and no positive result available by 72 h", "n": n_after_micro},
        {"step": "No active vasopressor overlap during 66-72 h", "n": n_after_vaso},
        {"step": "Broad-spectrum coverage during 48-72 h; final analytic cohort", "n": len(cohort)},
        {"step": "De-escalated/stopped", "n": int(cohort["A"].sum()) if len(cohort) else 0},
        {"step": "Continued broad-spectrum", "n": int((cohort["A"] == 0).sum()) if len(cohort) else 0},
    ])
    _write(flow, paths.tables / "cohort_flow.csv")
    _write(micro_win.head(2000), paths.audits / "microbiology_window_audit_sample.csv")

    if cfg["medications"].get("administration_sensitivity", True):
        admin = _medication_table(cfg, "med_admin")
        common = cohort.copy()
        admin_class = _classify(common, admin)
        comparison = cohort[["encounterid", "A", "deescalation_type"]].merge(
            admin_class[["encounterid", "A", "deescalation_type"]], on="encounterid", how="left", suffixes=("_order", "_admin")
        )
        _write(pd.crosstab(comparison["A_order"], comparison["A_admin"], dropna=False).reset_index(), paths.tables / "order_vs_administration_reclassification.csv")
        _write(comparison, paths.audits / "order_vs_administration_patient_level_local.csv")

    # Death is extracted for crude validation diagnostics. Full adjusted external inference is intentionally
    # blocked until PSU covariate/source mappings pass the harmonization audit.
    death = normalize_columns(read_table(resolve_path(cfg, "death")))
    if "death_date" in death:
        death["death_date"] = pd.to_datetime(death["death_date"], errors="coerce")
        first_death = death.groupby("patid")["death_date"].min()
        cohort["death_date"] = cohort["patid"].map(first_death)
        cohort["death_by_30d"] = (cohort["death_date"].notna() & (cohort["death_date"] >= cohort["analysis_time0"]) & (cohort["death_date"] <= cohort["analysis_time0"] + pd.Timedelta(days=int(cfg["clock"]["horizon_days"])))).astype(int)
        crude = cohort.groupby("A")["death_by_30d"].agg(["count", "sum", "mean"]).reset_index()
        _write(crude, paths.tables / "crude_mortality_diagnostic.csv")

    if cfg.get("outputs", {}).get("save_patient_level_cohort", True):
        cohort.to_csv(paths.run_dir / "external_cohort_unadjusted_local.csv", index=False)
    status = {
        "stage": "cohort_exposure_and_measurement_validation",
        "n": len(cohort),
        "adjusted_external_effect_estimate_ready": False,
        "next_requirement": "Freeze PSU lab/vital/SOFA/comorbidity mappings after inspecting actual datamart output; then reuse common IPTW/statistical modules.",
    }
    (paths.run_dir / "external_run_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    write_run_manifest(paths.run_dir, cfg, source="PSU/PCORnet local datamart", extra={"external_status": status})
    if cfg.get("outputs", {}).get("zip_run_directory", True):
        zip_run(paths.run_dir)
    return paths.run_dir

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .antibiotics import (
    ANAEROBIC_PATTERN,
    ANTI_MRSA_PATTERN,
    ANTIPSEUDOMONAL_PATTERN,
    CARBAPENEM_PATTERN,
    classify_treatment,
    definition_table,
    overlap_rows,
    prepare_coverage,
    raw_broad_mask,
    systemic_antibiotic_mask,
)
from .config import load_config, make_run_paths, resolve_mimic_source
from .features import (
    add_baseline_covariates,
    add_baseline_labs,
    add_baseline_vitals,
    add_lab_trajectories,
    add_severity_score,
    add_sofa_like,
    add_steroid_bmi,
    add_urine_output,
    add_vasopressor_windows,
    add_vital_trajectories,
    prepare_vasopressors,
)
from .microbiology import audit_counts, eligible_microbiology
from .mimic_io import diagnoses, inputevents, master_admissions, microbiology, prescriptions, procedures
from .outcomes import add_antibiotic_burden, add_antibiotic_free_days, add_hospital_free_days, add_mortality
from .provenance import copy_config, write_run_manifest, zip_run
from .specification import BINARY_VARS, CANDIDATE_PS_VARS, CONTINUOUS_VARS, PROGRESSIVE_MODELS
from .stats import (
    balance_table,
    bootstrap_iptw_ci,
    fit_stabilized_iptw,
    mean_difference,
    risks,
    weight_summary,
    weighted_mean,
)
from .tables import cohort_flow_row, table1

LOG = logging.getLogger(__name__)


def _write(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def _results_row(name: str, values: tuple, kind: str = "risk") -> pd.DataFrame:
    if kind == "risk":
        rt, rc, rd, rr = values
        return pd.DataFrame([{"analysis": name, "risk_deescalated_stopped": rt, "risk_continued": rc, "risk_difference": rd, "risk_ratio": rr}])
    mt, mc, md = values
    return pd.DataFrame([{"analysis": name, "mean_deescalated_stopped": mt, "mean_continued": mc, "mean_difference": md}])


def _add_antibiotic_intensity(d: pd.DataFrame, all_rx: pd.DataFrame, broad_rx: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    windows = d[["hadm_id", "first_broad_time", "decision_time"]].copy()
    windows["window_start"] = windows["first_broad_time"]
    windows["window_end"] = windows["decision_time"]
    any_pre = overlap_rows(all_rx, "hadm_id", "coverage_start", "coverage_stop", windows)
    broad_pre = overlap_rows(broad_rx, "hadm_id", "coverage_start", "coverage_stop", windows)

    def summarize(events: pd.DataFrame, prefix: str) -> None:
        if not len(events):
            d[f"{prefix}_abx_hours_pre72"] = 0.0
            d[f"{prefix}_abx_agents_pre72"] = 0.0
            return
        e = events.copy()
        e["overlap_start"] = e[["coverage_start", "window_start"]].max(axis=1)
        e["overlap_stop"] = e[["coverage_stop", "window_end"]].min(axis=1)
        e["hours"] = ((e["overlap_stop"] - e["overlap_start"]).dt.total_seconds() / 3600.0).clip(lower=0)
        e["drug_simple"] = e["drug_lower"].str.replace(r"\s+", " ", regex=True).str.strip()
        d[f"{prefix}_abx_hours_pre72"] = d["hadm_id"].map(e.groupby("hadm_id")["hours"].sum()).fillna(0.0)
        d[f"{prefix}_abx_agents_pre72"] = d["hadm_id"].map(e.groupby("hadm_id")["drug_simple"].nunique()).fillna(0.0)

    summarize(any_pre, "systemic")
    summarize(broad_pre, "broad")
    text = broad_pre.get("drug_lower", pd.Series("", index=broad_pre.index)).fillna("")
    for col, pattern in {
        "anti_mrsa_pre72": ANTI_MRSA_PATTERN,
        "antipseudomonal_pre72": ANTIPSEUDOMONAL_PATTERN,
        "carbapenem_pre72": CARBAPENEM_PATTERN,
    }.items():
        ids = set(broad_pre.loc[text.str.contains(pattern, na=False, regex=True), "hadm_id"].astype(int))
        d[col] = d["hadm_id"].astype(int).isin(ids).astype(int)
    text_any = any_pre.get("drug_lower", pd.Series("", index=any_pre.index)).fillna("")
    anaerobe_ids = set(any_pre.loc[text_any.str.contains(ANAEROBIC_PATTERN, na=False, regex=True), "hadm_id"].astype(int))
    d["anaerobic_coverage_pre72"] = d["hadm_id"].astype(int).isin(anaerobe_ids).astype(int)
    return d


def _add_diagnostic_intensity(d: pd.DataFrame, micro_win: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    mw = micro_win.loc[micro_win["hadm_id"].isin(set(d["hadm_id"]))].copy()
    mw["text"] = mw["spec_type_desc_lower"].fillna("") + " " + mw["test_name_lower"].fillna("")
    d["micro_records_pre72"] = d["hadm_id"].map(mw.groupby("hadm_id").size()).fillna(0.0)
    d["clinical_micro_records_pre72"] = d["hadm_id"].map(mw.loc[mw["clinical_micro"] == 1].groupby("hadm_id").size()).fillna(0.0)
    d["strict_culture_records_pre72"] = d["hadm_id"].map(mw.loc[mw["strict_culture_record"] == 1].groupby("hadm_id").size()).fillna(0.0)
    d["distinct_specimen_types_pre72"] = d["hadm_id"].map(mw.groupby("hadm_id")["spec_type_desc_lower"].nunique()).fillna(0.0)
    late_key = d[["hadm_id", "decision_time"]].rename(columns={"decision_time": "decision_time_for_late_window"})
    late_key["late_start"] = late_key["decision_time_for_late_window"] - pd.Timedelta(hours=24)
    tmp = mw.merge(late_key, on="hadm_id", how="left")
    late = tmp.loc[(tmp["specimen_time"] >= tmp["late_start"]) & (tmp["specimen_time"] <= tmp["decision_time_for_late_window"])]
    d["repeat_micro_48_72h"] = d["hadm_id"].isin(set(late["hadm_id"])).astype(int)
    for name, pattern in {
        "blood_culture_pre72": "blood",
        "respiratory_culture_pre72": "sputum|respiratory|bronch|lavage|tracheal",
        "urine_culture_pre72": "urine",
        "sterile_fluid_culture_pre72": "csf|cerebrospinal|pleural|peritoneal|synovial|fluid|tissue|wound",
    }.items():
        ids = set(mw.loc[mw["text"].str.contains(pattern, na=False, regex=True), "hadm_id"])
        d[name] = d["hadm_id"].isin(ids).astype(int)
    return d


def _primary_and_secondary_results(
    cohort_w: pd.DataFrame,
    cfg: dict,
    paths,
) -> dict[str, pd.DataFrame]:
    reps = int(cfg["bootstrap"]["primary_reps"])
    seed = int(cfg["bootstrap"]["seed"])
    outputs: dict[str, pd.DataFrame] = {}
    outcomes = [
        ("30-day mortality", "death_by_horizon", "risk"),
        ("hospital-free days", "hospital_free_days", "mean"),
        ("antibiotic-free days", "antibiotic_free_days", "mean"),
        ("normalized systemic antibiotic exposure", "normalized_antibiotic_exposure_30d", "mean"),
        ("normalized broad-spectrum exposure", "normalized_broad_antibiotic_exposure_30d", "mean"),
        ("late recurrent/persistent antibiotic-course use", "late_recurrent_or_persistent_abx_course_30d", "risk"),
    ]
    rows = []
    for i, (label, outcome, kind) in enumerate(outcomes):
        point = risks(cohort_w, outcome, "SW_A") if kind == "risk" else mean_difference(cohort_w, outcome, "SW_A")
        result = _results_row(label, point, kind)
        ci, boot = bootstrap_iptw_ci(cohort_w, CANDIDATE_PS_VARS, outcome, kind, reps, seed + i)
        _write(boot, paths.diagnostics / f"bootstrap_{outcome}.csv")
        _write(ci.assign(analysis=label), paths.tables / f"ci_{outcome}.csv")
        result["kind"] = kind
        if len(ci):
            primary_est = "risk_difference" if kind == "risk" else "mean_difference"
            hit = ci.loc[ci["estimand"] == primary_est]
            if len(hit):
                result["lower_95"] = float(hit.iloc[0]["lower_95"])
                result["upper_95"] = float(hit.iloc[0]["upper_95"])
                result["bootstrap_success"] = int(hit.iloc[0]["n_success"])
        rows.append(result)
    outputs["outcomes"] = pd.concat(rows, ignore_index=True)
    return outputs


def _progressive_adjustment(cohort: pd.DataFrame, cfg: dict, paths) -> pd.DataFrame:
    reps = int(cfg["bootstrap"]["progressive_reps"])
    seed = int(cfg["bootstrap"]["seed"])
    rows = []
    for i, spec in enumerate(PROGRESSIVE_MODELS):
        w, _, diag = fit_stabilized_iptw(cohort, spec["vars"])
        rt, rc, rd, rr = risks(w, "death_by_horizon", "SW_A")
        bal = balance_table(w, spec["vars"])
        ci, boot = bootstrap_iptw_ci(w, spec["vars"], "death_by_horizon", "risk", reps, seed + 100 + i)
        _write(boot, paths.diagnostics / f"progressive_bootstrap_M{i+1}.csv")
        rd_ci = ci.loc[ci["estimand"] == "risk_difference"] if len(ci) else pd.DataFrame()
        rr_ci = ci.loc[ci["estimand"] == "risk_ratio"] if len(ci) else pd.DataFrame()
        rows.append({
            "model": spec["model"],
            "n_vars_requested": len(spec["vars"]),
            "n_vars_used": len(diag["used_vars"]),
            "risk_deescalated_stopped": rt,
            "risk_continued": rc,
            "risk_difference": rd,
            "rd_lower_95": float(rd_ci.iloc[0]["lower_95"]) if len(rd_ci) else np.nan,
            "rd_upper_95": float(rd_ci.iloc[0]["upper_95"]) if len(rd_ci) else np.nan,
            "risk_ratio": rr,
            "rr_lower_95": float(rr_ci.iloc[0]["lower_95"]) if len(rr_ci) else np.nan,
            "rr_upper_95": float(rr_ci.iloc[0]["upper_95"]) if len(rr_ci) else np.nan,
            "max_post_smd": float(bal["after"].max()) if len(bal) else np.nan,
        })
    return pd.DataFrame(rows)


def _simple_sensitivity(
    df: pd.DataFrame,
    label: str,
    cfg: dict,
    reps: int | None = None,
    seed_offset: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if len(df) < 100 or df["A"].nunique() != 2:
        return pd.DataFrame([{"analysis": label, "n": len(df), "status": "not_run"}]), pd.DataFrame()
    w, _, _ = fit_stabilized_iptw(df, CANDIDATE_PS_VARS)
    rt, rc, rd, rr = risks(w, "death_by_horizon", "SW_A")
    out = pd.DataFrame([{"analysis": label, "n": len(w), "treated": int(w["A"].sum()), "control": int((w["A"] == 0).sum()), "risk_treated": rt, "risk_control": rc, "risk_difference": rd, "risk_ratio": rr, "status": "ok"}])
    if reps is None:
        reps = int(cfg["bootstrap"]["primary_reps"])
    ci, boot = bootstrap_iptw_ci(w, CANDIDATE_PS_VARS, "death_by_horizon", "risk", reps, int(cfg["bootstrap"]["seed"]) + seed_offset)
    if len(ci):
        rd_ci = ci.loc[ci["estimand"] == "risk_difference"]
        rr_ci = ci.loc[ci["estimand"] == "risk_ratio"]
        if len(rd_ci): out["rd_lower_95"] = float(rd_ci.iloc[0]["lower_95"]); out["rd_upper_95"] = float(rd_ci.iloc[0]["upper_95"])
        if len(rr_ci): out["rr_lower_95"] = float(rr_ci.iloc[0]["lower_95"]); out["rr_upper_95"] = float(rr_ci.iloc[0]["upper_95"])
    return out, boot


def run_mimic(config_path: str | Path) -> Path:
    cfg = load_config(config_path)
    source = resolve_mimic_source(cfg)
    paths = make_run_paths(cfg["output_root"], cfg["site"], cfg["analysis_version"])
    copy_config(cfg, paths.run_dir)
    write_run_manifest(paths.run_dir, cfg, source)
    LOG.info("MIMIC source: %s", source)
    LOG.info("Run directory: %s", paths.run_dir)

    flow = []
    master = master_admissions(source)
    cohort0 = master.loc[master["intime"].notna() & master["outtime"].notna() & (master["anchor_age"] >= 18)].copy()
    flow.append(cohort_flow_row("Adult ICU admissions with hospital data", len(cohort0)))

    rx_raw = prescriptions(source, cohort0["hadm_id"].dropna().astype(int).unique())
    fill_hours = int(cfg["medications"]["missing_stop_fill_hours"])
    rx = prepare_coverage(rx_raw, fill_hours)
    rx_broad = rx.loc[systemic_antibiotic_mask(rx, broad=True) & rx["coverage_start"].notna()].copy()
    rx_abx = rx.loc[systemic_antibiotic_mask(rx, broad=False) & rx["coverage_start"].notna()].copy()

    early = rx_broad.merge(cohort0[["hadm_id", "intime"]], on="hadm_id", how="inner")
    early = early.loc[(early["coverage_start"] >= early["intime"] - pd.Timedelta(hours=6)) & (early["coverage_start"] <= early["intime"] + pd.Timedelta(hours=24))]
    first = early.sort_values(["hadm_id", "coverage_start"]).drop_duplicates("hadm_id")[["hadm_id", "coverage_start"]].rename(columns={"coverage_start": "first_broad_time"})
    cohort = cohort0.merge(first, on="hadm_id", how="inner")
    cohort["decision_time"] = cohort["first_broad_time"] + pd.Timedelta(hours=int(cfg["clock"]["decision_hours"]))
    cohort["analysis_time0"] = cohort["first_broad_time"] + pd.Timedelta(hours=int(cfg["clock"]["classification_end_hours"]))
    death_dt = pd.to_datetime(cohort["deathtime"].fillna(cohort["dod"]), errors="coerce")
    cohort = cohort.loc[(death_dt.isna() | (death_dt > cohort["analysis_time0"])) & (cohort["dischtime"] > cohort["analysis_time0"])].copy()
    flow.append(cohort_flow_row("Early systemic IV broad-spectrum exposure; alive and hospitalized through 96 h", len(cohort)))

    micro_raw = microbiology(source, cohort["hadm_id"].dropna().astype(int).unique())
    cohort, micro_sets, micro_win = eligible_microbiology(cohort, micro_raw)
    flow.append(cohort_flow_row("Clinical microbiology sampled and no positive result available by 72 h", len(cohort)))
    _write(audit_counts(micro_sets), paths.audits / "microbiology_counts.csv")
    positive_missing = micro_win.loc[(micro_win["positive_organism_row"] == 1) & micro_win["result_available_time"].isna()].copy()
    _write(positive_missing.head(500), paths.audits / "positive_organism_rows_missing_result_time.csv")

    dx = diagnoses(source, cohort["hadm_id"].dropna().astype(int).unique())
    px = procedures(source, cohort["hadm_id"].dropna().astype(int).unique())
    inp = inputevents(source, cohort["stay_id"].dropna().astype(int).unique())
    cohort = add_baseline_covariates(cohort, dx, px)
    cohort, vaso = prepare_vasopressors(cohort, inp)
    cohort = cohort.loc[cohort["vasopressor_overlap_last6h_pre72"] == 0].copy()
    flow.append(cohort_flow_row("No active vasopressor overlap during 66-72 h", len(cohort)))
    cohort = add_baseline_labs(cohort, source)
    cohort = add_baseline_vitals(cohort, source)
    cohort = add_severity_score(cohort)

    pre = cohort[["hadm_id", "decision_time"]].copy(); pre["window_start"] = pre["decision_time"] - pd.Timedelta(hours=24); pre["window_end"] = pre["decision_time"]
    pre_broad = overlap_rows(rx_broad, "hadm_id", "coverage_start", "coverage_stop", pre)
    cohort = cohort.loc[cohort["hadm_id"].isin(set(pre_broad["hadm_id"]))].copy()
    flow.append(cohort_flow_row("Systemic IV broad-spectrum coverage during 48-72 h", len(cohort)))
    cohort_pre_exposure = cohort.copy()

    cohort = classify_treatment(cohort, rx_broad, rx_abx)
    horizon = int(cfg["clock"]["horizon_days"])
    cohort = add_mortality(cohort, "analysis_time0", horizon)
    cohort = add_hospital_free_days(cohort, "analysis_time0", horizon)
    cohort = add_antibiotic_free_days(cohort, rx_abx, "analysis_time0", horizon)

    cohort = add_lab_trajectories(cohort, source)
    cohort = add_vasopressor_windows(cohort, vaso)
    cohort = add_vital_trajectories(cohort, source)
    cohort = add_urine_output(cohort, source)
    cohort = _add_antibiotic_intensity(cohort, rx_abx, rx_broad)
    cohort = _add_diagnostic_intensity(cohort, micro_win)
    cohort = add_steroid_bmi(cohort, rx, source)
    cohort = add_sofa_like(cohort)
    cohort = add_antibiotic_burden(cohort, rx_abx, rx_broad, "analysis_time0", horizon)
    cohort["icu_discharge_within_24h_after_landmark"] = (pd.to_datetime(cohort["outtime"]) <= cohort["analysis_time0"] + pd.Timedelta(hours=24)).astype(int)
    cohort["hospital_discharge_within_24h_after_landmark"] = (pd.to_datetime(cohort["dischtime"]) <= cohort["analysis_time0"] + pd.Timedelta(hours=24)).astype(int)
    cohort["near_discharge_24h_after_landmark"] = ((cohort["icu_discharge_within_24h_after_landmark"] == 1) | (cohort["hospital_discharge_within_24h_after_landmark"] == 1)).astype(int)

    flow.append(cohort_flow_row("Final analytic cohort", len(cohort)))
    flow.append(cohort_flow_row("De-escalated/stopped", int(cohort["A"].sum())))
    flow.append(cohort_flow_row("Continued broad-spectrum", int((cohort["A"] == 0).sum())))
    _write(pd.DataFrame(flow), paths.tables / "cohort_flow.csv")

    t1u = table1(cohort, CONTINUOUS_VARS, BINARY_VARS)
    cohort_w, ps_fit, ps_diag = fit_stabilized_iptw(cohort, CANDIDATE_PS_VARS)
    t1w = table1(cohort_w, CONTINUOUS_VARS, BINARY_VARS, weight_col="SW_A")
    balance = balance_table(cohort_w, CANDIDATE_PS_VARS)
    weights = weight_summary(cohort_w)
    ps_dist = cohort_w.groupby("A")["ps_den"].describe(percentiles=[.01, .05, .10, .25, .50, .75, .90, .95, .99]).reset_index()
    _write(t1u, paths.tables / "table1_unweighted.csv"); _write(t1w, paths.tables / "table1_weighted.csv")
    _write(balance, paths.diagnostics / "balance_before_after.csv"); _write(weights, paths.diagnostics / "weight_summary.csv"); _write(ps_dist, paths.diagnostics / "propensity_score_distribution.csv")
    _write(ps_diag["standardization_table"], paths.diagnostics / "ps_standardization.csv")
    (paths.diagnostics / "ps_formula.txt").write_text(ps_diag["formula"], encoding="utf-8")

    estimates = _primary_and_secondary_results(cohort_w, cfg, paths)["outcomes"]
    _write(estimates, paths.tables / "primary_secondary_outcomes.csv")
    progressive = _progressive_adjustment(cohort, cfg, paths)
    _write(progressive, paths.tables / "progressive_adjustment.csv")

    # Weight-truncation sensitivity.
    trunc_rows = []
    for low, high in cfg["weighting"]["truncation_percentiles"]:
        lo = cohort_w["SW_A"].quantile(float(low) / 100); hi = cohort_w["SW_A"].quantile(float(high) / 100)
        wcol = f"SW_trunc_{low}_{high}"; tmp = cohort_w.copy(); tmp[wcol] = tmp["SW_A"].clip(lo, hi)
        rt, rc, rd, rr = risks(tmp, "death_by_horizon", wcol)
        trunc_rows.append({"analysis": f"weight truncation {low}-{high}", "low": low, "high": high, "risk_treated": rt, "risk_control": rc, "risk_difference": rd, "risk_ratio": rr})
    _write(pd.DataFrame(trunc_rows), paths.tables / "weight_truncation_sensitivity.csv")

    # Narrowing and stopping contrasts.
    narrow = cohort.loc[cohort["deescalation_type"] != "stopped_all_observed_systemic_antibiotics"].copy(); narrow["A"] = (narrow["deescalation_type"] == "narrowed_or_non_broad_only").astype(int)
    stop = cohort.loc[cohort["deescalation_type"] != "narrowed_or_non_broad_only"].copy(); stop["A"] = (stop["deescalation_type"] == "stopped_all_observed_systemic_antibiotics").astype(int)
    sens_rows = []
    for i, (label, sdf) in enumerate([("Narrowed/non-broad only vs continued broad", narrow), ("Stopped all observed systemic antibiotics vs continued broad", stop)]):
        out, boot = _simple_sensitivity(sdf, label, cfg, seed_offset=500 + i); sens_rows.append(out); _write(boot, paths.diagnostics / f"bootstrap_sensitivity_{i+1}.csv")

    # Near-discharge sensitivity.
    near = cohort.loc[cohort["near_discharge_24h_after_landmark"] == 0].copy()
    out, boot = _simple_sensitivity(near, "Exclude ICU/hospital discharge within 24 h after landmark", cfg, seed_offset=510); sens_rows.append(out); _write(boot, paths.diagnostics / "bootstrap_near_discharge.csv")

    # Strict test-name culture sensitivity. Rebuild from pre-microbiology eligible IDs while retaining final feature rows where possible.
    strict_ids = micro_sets["strict_eligible"]
    strict = cohort.loc[cohort["hadm_id"].astype(int).isin(strict_ids)].copy()
    out, boot = _simple_sensitivity(strict, "Strict test-name culture-only", cfg, seed_offset=520); sens_rows.append(out); _write(boot, paths.diagnostics / "bootstrap_strict_culture.csv")

    # Eventual culture-negative sensitivity.
    eventual_ids = micro_sets["eventual_culture_negative"]
    eventual = cohort.loc[cohort["hadm_id"].astype(int).isin(eventual_ids)].copy()
    out, boot = _simple_sensitivity(eventual, "Eventual culture-negative from specimens collected by 72 h", cfg, seed_offset=530); sens_rows.append(out); _write(boot, paths.diagnostics / "bootstrap_eventual_culture_negative.csv")
    sensitivity = pd.concat(sens_rows, ignore_index=True, sort=False)
    _write(sensitivity, paths.tables / "mortality_sensitivity_analyses.csv")

    # Missing-stop-time deterministic sensitivity reruns exposure eligibility/classification using fully featured cohort.
    stop_rows = []
    for stop_fill in cfg["medications"].get("stop_fill_sensitivity_hours", []):
        alt = prepare_coverage(rx_raw, int(stop_fill)); broad_alt = alt.loc[systemic_antibiotic_mask(alt, broad=True) & alt["coverage_start"].notna()].copy(); any_alt = alt.loc[systemic_antibiotic_mask(alt, broad=False) & alt["coverage_start"].notna()].copy()
        windows = cohort_pre_exposure[["hadm_id", "decision_time"]].copy(); windows["window_start"] = windows["decision_time"] - pd.Timedelta(hours=24); windows["window_end"] = windows["decision_time"]
        pre_alt = overlap_rows(broad_alt, "hadm_id", "coverage_start", "coverage_stop", windows)
        ids = set(pre_alt["hadm_id"]); alt_cohort = cohort.loc[cohort["hadm_id"].isin(ids)].copy(); alt_cohort = classify_treatment(alt_cohort, broad_alt, any_alt)
        if alt_cohort["A"].nunique() == 2:
            w, _, _ = fit_stabilized_iptw(alt_cohort, CANDIDATE_PS_VARS); rt, rc, rd, rr = risks(w, "death_by_horizon", "SW_A")
            stop_rows.append({"stop_fill_hours": stop_fill, "n": len(w), "risk_treated": rt, "risk_control": rc, "risk_difference": rd, "risk_ratio": rr})
    _write(pd.DataFrame(stop_rows), paths.tables / "stop_time_assumption_sensitivity.csv")

    # Audits and figure-ready datasets. No figure rendering occurs here.
    _write(definition_table(), paths.audits / "antibiotic_definition.csv")
    raw_broad = rx.loc[raw_broad_mask(rx)].copy(); excluded = raw_broad.loc[~systemic_antibiotic_mask(raw_broad, broad=True)].copy(); _write(excluded[[c for c in ["hadm_id", "drug", "route", "starttime", "stoptime"] if c in excluded]].head(1000), paths.audits / "excluded_broad_medication_rows_sample.csv")
    _write(progressive, paths.figures / "figure2_progressive_adjustment.csv")
    _write(sensitivity, paths.figures / "figure3_mortality_sensitivities.csv")
    _write(balance, paths.figures / "supplement_love_plot.csv")
    _write(cohort_w[["A", "ps_den", "SW_A"]], paths.figures / "supplement_ps_overlap_weights.csv")
    stewardship = estimates.loc[estimates["analysis"].isin(["hospital-free days", "antibiotic-free days", "normalized systemic antibiotic exposure", "normalized broad-spectrum exposure"])].copy(); _write(stewardship, paths.figures / "supplement_stewardship_outcomes.csv")

    # Patient-level analytic file is local only and gitignored.
    if cfg.get("outputs", {}).get("save_patient_level_cohort", True):
        cohort_w.to_csv(paths.run_dir / "analysis_cohort_weighted.csv", index=False)

    summary = {
        "n": len(cohort_w),
        "deescalated_stopped": int(cohort_w["A"].sum()),
        "continued": int((cohort_w["A"] == 0).sum()),
        "deaths": int(cohort_w["death_by_horizon"].sum()),
        "max_post_smd": float(balance["after"].max()) if len(balance) else None,
        "max_weight": float(cohort_w["SW_A"].max()),
    }
    (paths.run_dir / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_run_manifest(paths.run_dir, cfg, source, extra={"analysis_summary": summary})
    if cfg.get("outputs", {}).get("zip_run_directory", True):
        zip_run(paths.run_dir)
    return paths.run_dir

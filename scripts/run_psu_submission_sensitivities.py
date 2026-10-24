#!/usr/bin/env python3
"""Post hoc PSU submission checks using the frozen cohort and propensity model.

Analyses:
1. Patient-clustered percentile bootstrap of 30-day mortality, resampling unique
   patients and all their eligible encounters with PS refitting in every draw.
2. The calendar-day +30 mortality boundary check, excluding events on exactly the
   30th calendar date, with identical cohort, exposure, and fitted primary weights.

Neither check alters the prespecified/frozen primary estimate or its inference.
All identifiers remain in memory; only aggregate outputs are written.
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import pandas as pd

from run_psu_point_estimates import INJECT as POINT_INJECT
from run_psu_bootstrap_inference import _init_worker, _worker

REPS = 1000
SEED = 20261024
EXPECTED_N = 19841
EXPECTED_PATIENTS = 14569
EXPECTED_EVENTS = 2381
EXPECTED_RD = -0.025598893527824254

EXTRA_INJECT = r'''
    run_submission_sensitivities(df, binary, continuous, sw, effects, con, args.output_dir)
'''


def _weighted_risk(y, a, w, treated):
    take = (a == treated)
    return float(np.sum(y[take] * w[take]) / np.sum(w[take]))


def _measure(y, a, w):
    t = _weighted_risk(y, a, w, 1)
    c = _weighted_risk(y, a, w, 0)
    return {"deescalated_risk": t, "continued_risk": c,
            "risk_difference": t - c, "risk_ratio": t / c if c > 0 else None}


def _write_progress(current, total, phase):
    target = os.environ.get("RUNRELAY_PROGRESS_FILE")
    if not target:
        return
    p = Path(target)
    import time
    payload = {"schema_version": 1, "current": current, "total": total,
               "fraction": current / total, "unit": "bootstrap replicates",
               "phase": phase, "updated_at_epoch": time.time()}
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(json.dumps(payload), encoding="utf-8")
    os.replace(tmp, p)


def _cluster_tasks(group_rows):
    rng = np.random.default_rng(SEED)
    n_patients = len(group_rows)
    for rep in range(REPS):
        sample = rng.integers(0, n_patients, size=n_patients, dtype=np.int32)
        rows = np.concatenate([group_rows[int(j)] for j in sample])
        yield rep, rows


def run_submission_sensitivities(df, binary, continuous, sw, effects, con, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    assert len(df) == EXPECTED_N, "Frozen strict-cohort parity failure"
    assert df["patid"].nunique() == EXPECTED_PATIENTS, "Unique-patient parity failure"
    assert int(df["death_30d"].sum()) == EXPECTED_EVENTS, "Frozen death-event parity failure"

    primary_row = effects.loc[
        (effects["method"] == "stabilized_ate_iptw") &
        (effects["outcome"] == "death_30d")
    ].iloc[0]
    if abs(float(primary_row["difference_A1_minus_A0"]) - EXPECTED_RD) > 1e-4:
        raise RuntimeError("Primary frozen PSU RD parity failure")

    # Original outcome includes death on landmark_date + 30 (calendar-day
    # inclusive). The alternative uses dates +0 through +29 only.
    boundary = con.execute("""
        select patid, encounterid,
               case when death_date is not null
                       and death_date >= cast(anchor_ts + interval 96 hour as date)
                       and death_date < cast(anchor_ts + interval 96 hour as date)
                                         + 30
                    then 1 else 0 end as death_30d_calendar_exclusive
        from cohort
    """).fetchdf()
    df = df.merge(boundary, on=["patid", "encounterid"],
                  validate="one_to_one", how="left")
    assert int(df["death_30d_calendar_exclusive"].isna().sum()) == 0
    a = df["A"].astype(int).to_numpy()
    baseline = df["death_30d"].to_numpy(float)
    exclusive = df["death_30d_calendar_exclusive"].to_numpy(float)
    assert np.all(exclusive <= baseline), "Boundary coding must only remove events"
    baseline_measure = _measure(baseline, a, sw)
    exclusive_measure = _measure(exclusive, a, sw)

    # Preserve the frozen original 30-day label. The alternative is a
    # post hoc endpoint definition sensitivity, not a correction to source data.
    group_rows = [np.asarray(rows, dtype=np.int32)
                  for rows in df.groupby("patid", sort=False).indices.values()]
    if len(group_rows) != EXPECTED_PATIENTS:
        raise RuntimeError("Patient-group formation parity failure")

    _write_progress(0, REPS, "patient-cluster bootstrap")
    jobs = max(1, min((os.cpu_count() or 2) - 1, 8))
    try:
        ctx = get_context("fork")
    except ValueError:
        ctx = get_context("spawn")
    records = []
    errors = {}
    with ProcessPoolExecutor(max_workers=jobs, mp_context=ctx,
                             initializer=_init_worker,
                             initargs=(df, binary, continuous,
                                       [("death_30d", "binary"),
                                        ("death_30d_calendar_exclusive", "binary")])) as pool:
        for finished, (rows, err) in enumerate(
                pool.map(_worker, _cluster_tasks(group_rows), chunksize=1), start=1):
            if err:
                errors[err] = errors.get(err, 0) + 1
            else:
                records.extend(rows)
            if finished % 10 == 0 or finished == REPS:
                _write_progress(finished, REPS, "patient-cluster bootstrap")

    boot = pd.DataFrame(records)
    success_reps = int(boot["rep"].nunique()) if len(boot) else 0
    ci_rows = []
    if not boot.empty:
        for (method, outcome), group in boot.groupby(["method", "outcome"], sort=False):
            for metric in ("difference_A1_minus_A0", "risk_ratio_A1_over_A0"):
                z = pd.to_numeric(group[metric], errors="coerce").replace(
                    [np.inf, -np.inf], np.nan).dropna()
                ci_rows.append({
                    "method": method, "outcome": outcome, "estimand": metric,
                    "lower_95": float(z.quantile(0.025)),
                    "upper_95": float(z.quantile(0.975)),
                    "n_success": int(len(z)), "n_requested": REPS
                })

    pd.DataFrame(ci_rows).to_csv(
        outdir / "psu_patient_cluster_bootstrap_ci.csv", index=False)
    diagnostic = {
        "requested_replicates": REPS, "successful_replicates": success_reps,
        "failed_replicates": REPS - success_reps, "failure_types": errors,
        "seed": SEED, "max_parallel_workers": jobs,
        "resampling_unit": "unique patient, with all eligible encounters",
        "ps_refit_in_each_replicate": True,
        "frozen_primary_analysis_changed": False,
    }
    (outdir / "psu_submission_diagnostics.json").write_text(
        json.dumps(diagnostic, indent=2), encoding="utf-8")

    summary = {
        "type": "post_hoc_submission_robustness",
        "data_safety": "only derived aggregate output; IDs retained locally in memory",
        "cohort_n": len(df),
        "unique_patients": len(group_rows),
        "multiple_encounter_patients": int((df.groupby("patid").size() > 1).sum()),
        "primary_baseline_events": int(np.sum(baseline)),
        "exclusive_calendar_day_30_events": int(np.sum(exclusive)),
        "events_on_calendar_day_30": int(np.sum(baseline - exclusive)),
        "baseline_primary_weighted": baseline_measure,
        "exclusive_boundary_same_primary_weights": exclusive_measure,
        "primary_original_reported_rd": float(primary_row["difference_A1_minus_A0"]),
        "published_encounter_bootstrap_primary_rd_ci": [-0.03612403394167609, -0.015495359974566778],
        "published_encounter_bootstrap_primary_rr_ci": [0.7213435259606157, 0.8757988404579597],
        "cluster_bootstrap_ci": [r for r in ci_rows if
             r["method"] == "stabilized_ate_iptw"],
        "limitations": [
            "Calendar-date boundary check is not exact 30x24-hour follow-up.",
            "New boundary outcome only changes deaths on landmark calendar date +30.",
            "The primary frozen endpoint and encounter-level confidence intervals are unchanged.",
            "Patient-level bootstrap is post hoc and intended to test uncertainty robustness."
        ]
    }
    (outdir / "psu_submission_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({
        "cohort_n": summary["cohort_n"],
        "unique_patients": summary["unique_patients"],
        "events_on_calendar_day_30": summary["events_on_calendar_day_30"],
        "successful_replicates": success_reps
    }))
    if success_reps < 950:
        raise RuntimeError("Too many failed patient-level bootstrap replicates")


def main():
    source_path = Path(__file__).resolve().parent / "audit_psu_ps_balance.py"
    source = source_path.read_text(encoding="utf-8")
    marker = "    max_pre=float(bdf.abs_pre_smd.max()); max_post=float(bdf.abs_post_smd.max()); worst=str(bdf.iloc[0].variable)\n"
    if source.count(marker) != 1:
        raise RuntimeError("Frozen PS source insertion marker missing or ambiguous")
    source = source.replace(
        marker, POINT_INJECT + "\n" + EXTRA_INJECT + "\n" + marker, 1)
    ns = {
        "__name__": "__main__",
        "__file__": str(source_path),
        "run_submission_sensitivities": run_submission_sensitivities
    }
    exec(compile(source, str(source_path), "exec"), ns, ns)


if __name__ == "__main__":
    main()

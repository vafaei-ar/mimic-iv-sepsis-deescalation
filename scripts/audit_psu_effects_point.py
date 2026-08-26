#!/usr/bin/env python3
"""Aggregate-only PSU treatment-effect point estimates after PS/outcome freeze.

This task deliberately reuses the already frozen propensity-score diagnostic code in
process and captures its local analysis frame/weights at function return. It then adds
the separately frozen PSU date-level outcomes and reports only aggregate weighted
point estimates and weighting sensitivities. No bootstrap or outcome-driven model
changes occur here.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import audit_psu_ps_balance as psmod


def wmean(x, w):
    x = np.asarray(x, float); w = np.asarray(w, float)
    m = np.isfinite(x) & np.isfinite(w)
    return float(np.sum(x[m] * w[m]) / np.sum(w[m])) if m.any() and np.sum(w[m]) > 0 else np.nan


def ess(w):
    w = np.asarray(w, float)
    return float((w.sum() ** 2) / np.sum(w ** 2)) if np.sum(w ** 2) > 0 else np.nan


def run_frozen_ps(data_root: Path, work_dir: Path):
    captured = {}

    def tracer(frame, event, arg):
        if event == "return" and frame.f_code.co_name == "main" and frame.f_globals.get("__name__") == psmod.__name__:
            captured.update(frame.f_locals)
        return tracer

    old_argv = sys.argv[:]
    old_trace = sys.gettrace()
    sys.argv = ["audit_psu_ps_balance.py", str(data_root), "--output-dir", str(work_dir)]
    try:
        sys.settrace(tracer)
        psmod.main()
    finally:
        sys.settrace(old_trace)
        sys.argv = old_argv
    required = ["df", "sw", "ps", "y", "con", "broad", "nonbroad", "p_start_expr", "p_end_expr", "pp", "pe"]
    missing = [k for k in required if k not in captured]
    if missing:
        raise RuntimeError(f"Frozen PS capture failed; missing locals: {missing}")
    return captured


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    internal = args.output_dir / "_internal_ps_diagnostic"
    if internal.exists():
        shutil.rmtree(internal)
    internal.mkdir(parents=True)
    z = run_frozen_ps(args.data_root, internal)

    df = z["df"].copy()
    sw = np.asarray(z["sw"], float)
    ps = np.asarray(z["ps"], float)
    y = np.asarray(z["y"], int)
    con = z["con"]
    broad = z["broad"]; nonbroad = z["nonbroad"]
    p_start = z["p_start_expr"]; p_end = z["p_end_expr"]
    pp, pe = z["pp"], z["pe"]
    q, qi = psmod.q, psmod.qi

    if len(df) != 19841 or int(y.sum()) != 5346:
        raise RuntimeError(f"Frozen PS parity failure: n={len(df)}, A1={int(y.sum())}")

    # Frozen outcome construction, matching PSU-FINAL-OUTCOME-FREEZE-0001.
    con.execute("drop table if exists effect_days")
    con.execute("drop table if exists effect_day_flags")
    con.execute("drop table if exists effect_day_triplets")
    con.execute("drop table if exists effect_abx_summary")
    con.execute("drop table if exists effect_outcomes")
    con.execute("""
      create temp table effect_days as
      select c.patid,c.encounterid,
             cast(c.anchor_ts+interval 96 hour as date) analysis_date,
             c.discharge_date,c.death_date,
             gs.day_idx::integer day_idx,
             cast(c.anchor_ts+interval 96 hour as date)+gs.day_idx::integer day_date
      from cohort c cross join generate_series(0,29) as gs(day_idx)
    """)
    any_systemic = f"(({broad}) or ({nonbroad}))"
    con.execute(f"""
      create temp table effect_day_flags as
      select d.patid,d.encounterid,d.day_idx,d.day_date,
        max(case when {any_systemic} and {p_start}<=d.day_date and {p_end}>=d.day_date then 1 else 0 end) any_abx,
        max(case when {broad} and {p_start}<=d.day_date and {p_end}>=d.day_date then 1 else 0 end) broad_abx
      from effect_days d left join p
        on cast(p.{qi(pp)} as varchar)=d.patid and cast(p.{qi(pe)} as varchar)=d.encounterid
      group by 1,2,3,4
    """)
    con.execute("""
      create temp table effect_day_triplets as
      select *,lead(any_abx,1,0) over(partition by patid,encounterid order by day_idx) any_abx_d1,
               lead(any_abx,2,0) over(partition by patid,encounterid order by day_idx) any_abx_d2
      from effect_day_flags
    """)
    con.execute("""
      create temp table effect_abx_summary as
      select patid,encounterid,
        sum(any_abx)::integer antibiotic_days_30d,
        sum(broad_abx)::integer broad_antibiotic_days_30d,
        max(case when day_idx>=7 and day_idx<=27 and any_abx=1 and any_abx_d1=1 and any_abx_d2=1 then 1 else 0 end)::integer late_recurrent_or_persistent_abx_course_30d
      from effect_day_triplets group by 1,2
    """)
    con.execute("""
      create temp table effect_outcomes as
      select c.patid,c.encounterid,
        case when c.death_date is not null and c.death_date>=cast(c.anchor_ts+interval 96 hour as date)
               and c.death_date<=cast(c.anchor_ts+interval 96 hour as date)+30 then 1 else 0 end death_30d,
        case when c.death_date is not null and c.death_date>=cast(c.anchor_ts+interval 96 hour as date)
               and c.death_date<=cast(c.anchor_ts+interval 96 hour as date)+30 then 0.0
             else greatest(0.0,30.0-least(30.0,greatest(0.0,date_diff('day',cast(c.anchor_ts+interval 96 hour as date),c.discharge_date)::double))) end hospital_free_days_30d,
        a.antibiotic_days_30d,
        case when c.death_date is not null and c.death_date>=cast(c.anchor_ts+interval 96 hour as date)
               and c.death_date<=cast(c.anchor_ts+interval 96 hour as date)+30 then 0.0
             else greatest(0.0,30.0-a.antibiotic_days_30d::double) end antibiotic_free_days_30d,
        a.broad_antibiotic_days_30d,
        case when c.death_date is not null and c.death_date>=cast(c.anchor_ts+interval 96 hour as date)
               and c.death_date<=cast(c.anchor_ts+interval 96 hour as date)+30
             then greatest(0.1,least(30.0,date_diff('day',cast(c.anchor_ts+interval 96 hour as date),c.death_date)::double))
             else 30.0 end days_alive_30d,
        a.late_recurrent_or_persistent_abx_course_30d,
        least(30.0,greatest(0.0,date_diff('day',cast(c.anchor_ts+interval 96 hour as date),c.discharge_date)::double)) observable_hospital_days_after_landmark
      from cohort c join effect_abx_summary a using(patid,encounterid)
    """)
    out = con.execute("select * from effect_outcomes").fetchdf()
    out["normalized_antibiotic_exposure_30d"] = (out["antibiotic_days_30d"] / out["days_alive_30d"]).clip(0,1)
    out["normalized_broad_antibiotic_exposure_30d"] = (out["broad_antibiotic_days_30d"] / out["days_alive_30d"]).clip(0,1)
    d = df[["patid","encounterid","A"]].merge(out,on=["patid","encounterid"],how="inner",validate="one_to_one")
    if len(d) != len(df):
        raise RuntimeError(f"Outcome merge parity failure: {len(d)} vs {len(df)}")
    # Align captured PS/weights to df merge order; merge preserves left order.
    d["ps"] = ps; d["sw"] = sw

    # Verify the frozen overall outcome distributions before exposing group effects.
    if int(d.death_30d.sum()) != 2381:
        raise RuntimeError(f"Outcome freeze parity failure: death_30d={int(d.death_30d.sum())}, expected 2381")

    lo1, hi1 = np.quantile(sw, [0.01,0.99])
    lo25, hi25 = np.quantile(sw, [0.025,0.975])
    weight_sets = {
        "primary_stabilized_iptw": sw,
        "overlap": np.where(y==1, 1-ps, ps),
        "trunc_1_99": np.clip(sw, lo1, hi1),
        "trunc_2_5_97_5": np.clip(sw, lo25, hi25),
    }
    binary = ["death_30d","late_recurrent_or_persistent_abx_course_30d"]
    continuous = ["hospital_free_days_30d","antibiotic_free_days_30d","normalized_antibiotic_exposure_30d","normalized_broad_antibiotic_exposure_30d"]
    diagnostic = ["observable_hospital_days_after_landmark"]
    rows = []
    for method,w in weight_sets.items():
        for name in binary:
            x=d[name].to_numpy(float); r1=wmean(x[y==1],w[y==1]); r0=wmean(x[y==0],w[y==0])
            rows.append({"method":method,"outcome":name,"type":"binary","deescalated":r1,"continued":r0,"contrast":"risk_difference","estimate":r1-r0,"risk_ratio":r1/r0 if r0>0 else None})
        for name in continuous+diagnostic:
            x=d[name].to_numpy(float); m1=wmean(x[y==1],w[y==1]); m0=wmean(x[y==0],w[y==0])
            rows.append({"method":method,"outcome":name,"type":"diagnostic" if name in diagnostic else "continuous","deescalated":m1,"continued":m0,"contrast":"mean_difference","estimate":m1-m0,"risk_ratio":None})
    res=pd.DataFrame(rows)
    res.to_csv(args.output_dir/"point_estimates.csv",index=False)

    wrows=[]
    for method,w in weight_sets.items():
        wrows.append({"method":method,"ess_deescalated":ess(w[y==1]),"ess_continued":ess(w[y==0]),"max_weight":float(np.max(w)),"p99_weight":float(np.quantile(w,.99))})
    pd.DataFrame(wrows).to_csv(args.output_dir/"weighting_sensitivity_diagnostics.csv",index=False)

    prim=res[res.method=="primary_stabilized_iptw"].copy()
    summary={
      "privacy_mode":"aggregate_only",
      "strict_cohort_n":19841,
      "deescalated_n":5346,
      "continued_n":14495,
      "ps_specification":"identical to frozen PSU-PS-BALANCE-DIAGNOSTIC-0002; no outcome-driven changes",
      "outcome_specification":"identical to frozen PSU-FINAL-OUTCOME-FREEZE-0001",
      "primary_weighting":"stabilized ATE IPTW",
      "primary_death_rd":float(prim.loc[prim.outcome=="death_30d","estimate"].iloc[0]),
      "primary_death_rr":float(prim.loc[prim.outcome=="death_30d","risk_ratio"].iloc[0]),
      "bootstrap":"not_run_in_this_task",
      "next_step":"If point estimates and diagnostics are coherent, run prespecified bootstrap confidence intervals without changing PS or outcome definitions."
    }
    (args.output_dir/"summary.json").write_text(json.dumps(summary,indent=2))

    shutil.rmtree(internal,ignore_errors=True)


if __name__ == "__main__":
    main()

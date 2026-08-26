#!/usr/bin/env python3
"""Aggregate-only final outcome freeze for the modified PSU external replication.

Builds the established strict 96-hour landmark cohort, defines outcomes from the
96-hour analysis landmark, and exports only overall aggregate distributions and a
frozen outcome dictionary. Exposure groups are deliberately not joined in this task,
so no unadjusted or adjusted treatment effects are computed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb
import pandas as pd

from audit_psu_final_covariate_freeze import (
    BROAD_PATTERN,
    NON_SYSTEMIC_PATTERN,
    find,
    first,
    parse_legacy_codes,
    q,
    qi,
    ts,
)

NON_BROAD_PATTERN = (
    "ceftriaxone|cefazolin|ampicillin|amoxicillin|doxycycline|azithromycin|"
    "metronidazole|clindamycin|cephalexin|ciprofloxacin|levofloxacin|gentamicin|tobramycin"
)
MIN_CELL = 11


def safe(n):
    if n is None:
        return None
    n = int(n)
    return n if n >= MIN_CELL else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    root = args.data_root

    paths = {
        "s": find(root, ["sepsis_encounter"]),
        "p": find(root, ["prescribing"]),
        "d": find(root, ["death"]),
    }
    con = duckdb.connect()
    con.execute("pragma threads=4")
    for k, pth in paths.items():
        con.execute(f"create view {k} as select * from read_parquet({q(str(pth))})")
    cols = {k: set(con.execute(f"describe {k}").fetchdf()["column_name"].astype(str)) for k in paths}

    sp, se = first(cols["s"], ["PATID"]), first(cols["s"], ["ENCOUNTERID", "ENCOUNTER_ID"])
    ad, dc = first(cols["s"], ["ADMIT_DATE"]), first(cols["s"], ["DISCHARGE_DATE"])
    pp, pe = first(cols["p"], ["PATID"]), first(cols["p"], ["ENCOUNTERID", "ENCOUNTER_ID"])
    pod, pot = first(cols["p"], ["RX_ORDER_DATE"]), first(cols["p"], ["RX_ORDER_TIME"])
    pstart, pend = first(cols["p"], ["RX_START_DATE", "START_DATE"]), first(cols["p"], ["RX_END_DATE", "END_DATE"])
    pn, pc, pr = first(cols["p"], ["RAW_RX_MED_NAME"]), first(cols["p"], ["RXNORM_CUI"]), first(cols["p"], ["RX_ROUTE"])
    dp, dd = first(cols["d"], ["PATID"]), first(cols["d"], ["DEATH_DATE"])
    req = [sp, se, ad, dc, pp, pe, pod, pn, pc, dp, dd]
    if any(x is None for x in req):
        raise RuntimeError("Required PSU outcome-freeze fields are missing")

    legacy = root / "PCORnet" / "code" / "config" / "codes_antibiotics.py"
    inc, exc = parse_legacy_codes(legacy)
    incsql = ",".join(q(x) for x in sorted(inc)) or "''"
    excsql = ",".join(q(x) for x in sorted(exc)) or "''"

    pts = ts("p", pod, pot)
    pname = f"lower(coalesce(cast(p.{qi(pn)} as varchar),''))"
    pcode = f"trim(coalesce(cast(p.{qi(pc)} as varchar),''))"
    route = f"upper(trim(coalesce(cast(p.{qi(pr)} as varchar),'')))" if pr else "''"
    bad_route = "('ORAL','RESPIRATORY_TRACT','INHALATION')"
    broad = (
        f"(({pcode} in ({incsql}) or regexp_matches({pname},{q(BROAD_PATTERN)})) "
        f"and {pcode} not in ({excsql}) and not regexp_matches({pname},{q(NON_SYSTEMIC_PATTERN)}) "
        f"and {route} not in {bad_route})"
    )
    any_systemic = (
        f"((({pcode} in ({incsql}) or regexp_matches({pname},{q(BROAD_PATTERN)})) "
        f"or regexp_matches({pname},{q(NON_BROAD_PATTERN)})) "
        f"and {pcode} not in ({excsql}) and not regexp_matches({pname},{q(NON_SYSTEMIC_PATTERN)}) "
        f"and {route} not in {bad_route})"
    )

    con.execute(f"""
      create temp table base as
      select distinct cast({qi(sp)} as varchar) patid,
             cast({qi(se)} as varchar) encounterid,
             try_cast({qi(ad)} as date) admit_date,
             try_cast({qi(dc)} as date) discharge_date
      from s
    """)
    con.execute(f"""
      create temp table anchors as
      select cast(p.{qi(pp)} as varchar) patid,
             cast(p.{qi(pe)} as varchar) encounterid,
             min({pts}) anchor_ts
      from p join base b
        on cast(p.{qi(pp)} as varchar)=b.patid
       and cast(p.{qi(pe)} as varchar)=b.encounterid
      where {broad}
        and {pts}>=cast(b.admit_date as timestamp)
        and {pts}<cast(b.admit_date as timestamp)+interval 24 hour
      group by 1,2
    """)
    con.execute(f"""
      create temp table deaths as
      select cast({qi(dp)} as varchar) patid,
             min(try_cast({qi(dd)} as date)) death_date
      from d group by 1
    """)
    con.execute("""
      create temp table cohort as
      select a.patid,a.encounterid,a.anchor_ts,
             a.anchor_ts+interval 96 hour analysis_time0,
             cast(a.anchor_ts+interval 96 hour as date) analysis_date,
             b.admit_date,b.discharge_date,d.death_date
      from anchors a join base b using(patid,encounterid)
      left join deaths d using(patid)
      where (b.discharge_date is null or b.discharge_date>cast(a.anchor_ts+interval 96 hour as date))
        and (d.death_date is null or d.death_date>cast(a.anchor_ts+interval 96 hour as date))
    """)
    n = int(con.execute("select count(*) from cohort").fetchone()[0])
    if n != 19841:
        raise RuntimeError(f"Strict cohort parity failure: expected 19841, observed {n}")

    # Date-level prescription intervals. Missing end date is conservatively treated as same-day.
    p_start = f"coalesce(try_cast(p.{qi(pstart)} as date),try_cast(p.{qi(pod)} as date))" if pstart else f"try_cast(p.{qi(pod)} as date)"
    p_end = f"coalesce(try_cast(p.{qi(pend)} as date),{p_start})" if pend else p_start

    # One local row per encounter-day over the 30-day horizon. No row-level output is exported.
    con.execute("""
      create temp table days as
      select c.patid,c.encounterid,c.analysis_date,c.discharge_date,c.death_date,
             gs.day_idx::integer day_idx,
             c.analysis_date + gs.day_idx::integer day_date
      from cohort c cross join generate_series(0,29) as gs(day_idx)
    """)
    con.execute(f"""
      create temp table day_flags as
      select d.patid,d.encounterid,d.day_idx,d.day_date,
             max(case when {any_systemic} and {p_start}<=d.day_date and {p_end}>=d.day_date then 1 else 0 end) any_abx,
             max(case when {broad} and {p_start}<=d.day_date and {p_end}>=d.day_date then 1 else 0 end) broad_abx
      from days d left join p
        on cast(p.{qi(pp)} as varchar)=d.patid
       and cast(p.{qi(pe)} as varchar)=d.encounterid
      group by 1,2,3,4
    """)
    con.execute("""
      create temp table day_triplets as
      select *,
             lead(any_abx,1,0) over(partition by patid,encounterid order by day_idx) any_abx_d1,
             lead(any_abx,2,0) over(partition by patid,encounterid order by day_idx) any_abx_d2
      from day_flags
    """)
    con.execute("""
      create temp table abx_summary as
      select patid,encounterid,
             sum(any_abx)::integer antibiotic_days_30d,
             sum(broad_abx)::integer broad_antibiotic_days_30d,
             max(case when day_idx>=7 and day_idx<=27 and any_abx=1 and any_abx_d1=1 and any_abx_d2=1 then 1 else 0 end)::integer late_recurrent_or_persistent_abx_course_30d
      from day_triplets group by 1,2
    """)

    con.execute("""
      create temp table outcomes as
      select c.patid,c.encounterid,
        case when c.death_date is not null
               and c.death_date>=c.analysis_date
               and c.death_date<=c.analysis_date+30 then 1 else 0 end death_30d,
        case when c.death_date is not null
               and c.death_date>=c.analysis_date
               and c.death_date<=c.analysis_date+30 then 0.0
             else greatest(0.0,30.0-least(30.0,greatest(0.0,date_diff('day',c.analysis_date,c.discharge_date)::double))) end hospital_free_days_30d,
        a.antibiotic_days_30d,
        case when c.death_date is not null
               and c.death_date>=c.analysis_date
               and c.death_date<=c.analysis_date+30 then 0.0
             else greatest(0.0,30.0-a.antibiotic_days_30d::double) end antibiotic_free_days_30d,
        a.broad_antibiotic_days_30d,
        case when c.death_date is not null
               and c.death_date>=c.analysis_date
               and c.death_date<=c.analysis_date+30
             then greatest(0.1,least(30.0,date_diff('day',c.analysis_date,c.death_date)::double))
             else 30.0 end days_alive_30d,
        a.late_recurrent_or_persistent_abx_course_30d,
        least(30.0,greatest(0.0,date_diff('day',c.analysis_date,c.discharge_date)::double)) observable_hospital_days_after_landmark
      from cohort c join abx_summary a using(patid,encounterid)
    """)
    con.execute("""
      alter table outcomes add column normalized_antibiotic_exposure_30d double;
      update outcomes set normalized_antibiotic_exposure_30d=least(1.0,greatest(0.0,antibiotic_days_30d/days_alive_30d));
      alter table outcomes add column normalized_broad_antibiotic_exposure_30d double;
      update outcomes set normalized_broad_antibiotic_exposure_30d=least(1.0,greatest(0.0,broad_antibiotic_days_30d/days_alive_30d));
    """)

    dictionary = pd.DataFrame([
        {"outcome":"death_30d","status":"primary_harmonized_date_level","source":"DEATH","definition":"Death date from analysis_time0=first broad order+96h through +30 days, inclusive; cohort already excludes deaths through the 96h landmark.","harmonization_note":"Date-level PSU death timing; closest harmonized mortality outcome."},
        {"outcome":"hospital_free_days_30d","status":"secondary_approximate","source":"sepsis_encounter + DEATH","definition":"30 minus calendar days from 96h landmark date to index discharge, clipped 0-30; set to 0 for death within 30 days.","harmonization_note":"Calendar-date approximation of MIMIC hour-level HFD."},
        {"outcome":"antibiotic_free_days_30d","status":"secondary_modified","source":"PRESCRIBING + DEATH","definition":"30 minus index-encounter calendar days with frozen systemic-antibiotic proxy after the 96h landmark; set to 0 for death within 30 days.","harmonization_note":"PRESCRIBING is date-level and measures ordered therapy, not verified administration."},
        {"outcome":"normalized_antibiotic_exposure_30d","status":"secondary_modified","source":"PRESCRIBING + DEATH","definition":"Calendar days with frozen systemic-antibiotic proxy divided by days alive through 30 days, clipped 0-1.","harmonization_note":"Date-level ordered-treatment burden."},
        {"outcome":"normalized_broad_antibiotic_exposure_30d","status":"secondary_modified","source":"PRESCRIBING + DEATH","definition":"Calendar days with frozen broad-spectrum proxy divided by days alive through 30 days, clipped 0-1.","harmonization_note":"Date-level ordered broad-spectrum burden."},
        {"outcome":"late_recurrent_or_persistent_abx_course_30d","status":"exploratory_modified","source":"PRESCRIBING","definition":"At least 3 consecutive calendar days with systemic-antibiotic proxy beginning on day 7 or later after the 96h landmark.","harmonization_note":"Strongly vulnerable to index-encounter observation/discharge differences."},
        {"outcome":"observable_hospital_days_after_landmark","status":"diagnostic_only","source":"sepsis_encounter","definition":"Calendar days from 96h landmark date to index discharge, clipped 0-30.","harmonization_note":"Used to describe differential outcome observability; not a treatment effect target."},
    ])
    dictionary.to_csv(args.output_dir / "outcome_dictionary.csv", index=False)

    binary = ["death_30d","late_recurrent_or_persistent_abx_course_30d"]
    continuous = [
        "hospital_free_days_30d","antibiotic_free_days_30d","antibiotic_days_30d",
        "broad_antibiotic_days_30d","normalized_antibiotic_exposure_30d",
        "normalized_broad_antibiotic_exposure_30d","observable_hospital_days_after_landmark",
    ]
    rows = []
    for name in binary:
        r = con.execute(f"select count(*) n,sum({name}) events,avg({name}::double) mean from outcomes").fetchone()
        rows.append({"outcome":name,"type":"binary","n":safe(r[0]),"events":safe(r[1]),"mean":r[2],"sd":None,"p25":None,"median":None,"p75":None,"min":0.0,"max":1.0})
    for name in continuous:
        r = con.execute(f"select count({name}),avg({name}),stddev_samp({name}),quantile_cont({name},0.25),quantile_cont({name},0.5),quantile_cont({name},0.75),min({name}),max({name}) from outcomes").fetchone()
        rows.append({"outcome":name,"type":"continuous","n":safe(r[0]),"events":None,"mean":r[1],"sd":r[2],"p25":r[3],"median":r[4],"p75":r[5],"min":r[6],"max":r[7]})
    pd.DataFrame(rows).to_csv(args.output_dir / "outcome_overall_summary.csv", index=False)

    obs = con.execute("""
      select
        count(*) n,
        count(*) filter(where discharge_date is not null) discharge_known,
        count(*) filter(where death_date is not null) any_death_date_known,
        count(*) filter(where death_date is not null and death_date>=analysis_date and death_date<=analysis_date+30) death_30d
      from cohort
    """).fetchone()
    pd.DataFrame([{
        "cohort_n":safe(obs[0]),
        "discharge_known":safe(obs[1]),
        "any_death_date_known":safe(obs[2]),
        "death_30d":safe(obs[3]),
        "prescribing_outcome_clock":"calendar-date intervals; missing RX_END_DATE treated as same-day",
        "analysis_time0":"first systemic broad-spectrum PRESCRIBING proxy + 96 hours",
    }]).to_csv(args.output_dir / "outcome_observability.csv", index=False)

    summary = {
        "privacy_mode":"aggregate_only",
        "strict_cohort_n":safe(n),
        "analysis_time0":"first systemic broad-spectrum PRESCRIBING proxy + 96 hours",
        "primary_outcome":"death_30d",
        "primary_outcome_clock":"date-level",
        "secondary_outcomes":["hospital_free_days_30d","antibiotic_free_days_30d","normalized_antibiotic_exposure_30d","normalized_broad_antibiotic_exposure_30d"],
        "exploratory_outcome":"late_recurrent_or_persistent_abx_course_30d",
        "important_limitation":"PSU medication and discharge outcomes are calendar-date/index-encounter approximations; they are not exact MIMIC hour-level equivalents.",
        "guardrail":"No exposure groups, propensity scores, outcomes-by-treatment summaries, or treatment effects are computed in this task.",
        "next_step":"If distributions are plausible, freeze PSU outcomes and run treatment-effect estimation using the already frozen PS specification plus prespecified weighting sensitivities.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

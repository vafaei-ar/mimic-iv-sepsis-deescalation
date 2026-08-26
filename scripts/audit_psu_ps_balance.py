#!/usr/bin/env python3
"""Aggregate-only PSU propensity-score balance diagnostic.

Builds the frozen strict modified PSU cohort, primary PRESCRIBING exposure, and the
pre-specified frozen PSU covariate set. Fits a propensity model and reports overlap,
weights, missingness and standardized mean differences. No outcomes or treatment
effects are computed or exported.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import statsmodels.api as sm

from audit_psu_final_covariate_freeze import (
    BROAD_PATTERN,
    NON_SYSTEMIC_PATTERN,
    LABS,
    OBS,
    VASO_TERMS,
    VASO_EXCLUDE,
    DX_PATTERNS,
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


def weighted_mean(x, w):
    m = np.isfinite(x) & np.isfinite(w)
    if not m.any() or w[m].sum() <= 0:
        return np.nan
    return np.sum(x[m] * w[m]) / np.sum(w[m])


def weighted_var(x, w):
    m = np.isfinite(x) & np.isfinite(w)
    if not m.any() or w[m].sum() <= 0:
        return np.nan
    mu = weighted_mean(x[m], w[m])
    return np.sum(w[m] * (x[m] - mu) ** 2) / np.sum(w[m])


def smd(x, a, w=None):
    x = np.asarray(x, float); a = np.asarray(a, int)
    if w is None:
        w = np.ones(len(x), float)
    else:
        w = np.asarray(w, float)
    m1, m0 = a == 1, a == 0
    mu1, mu0 = weighted_mean(x[m1], w[m1]), weighted_mean(x[m0], w[m0])
    v1, v0 = weighted_var(x[m1], w[m1]), weighted_var(x[m0], w[m0])
    den = np.sqrt((v1 + v0) / 2.0)
    if not np.isfinite(den) or den == 0:
        return 0.0 if np.isfinite(mu1) and np.isfinite(mu0) and mu1 == mu0 else np.nan
    return (mu1 - mu0) / den


def ess(w):
    w = np.asarray(w, float)
    return float((w.sum() ** 2) / np.sum(w ** 2)) if np.sum(w ** 2) > 0 else np.nan


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
        "dem": find(root, ["sepsis_demographic", "demographic"]),
        "dx": find(root, ["sepsis_diagnosis", "diagnosis"]),
        "lab": find(root, ["lab_reduced", "lab_result_cm"]),
        "obs": find(root, ["obs_clin"]),
        "vit": find(root, ["sepsis_vital", "vital"]),
        "med": find(root, ["med_admin"]),
    }
    con = duckdb.connect(); con.execute("pragma threads=4")
    for k, pth in paths.items():
        con.execute(f"create view {k} as select * from read_parquet({q(str(pth))})")
    cols = {k: set(con.execute(f"describe {k}").fetchdf()["column_name"].astype(str)) for k in paths}

    sp, se = first(cols["s"],["PATID"]), first(cols["s"],["ENCOUNTERID","ENCOUNTER_ID"])
    ad, dc = first(cols["s"],["ADMIT_DATE"]), first(cols["s"],["DISCHARGE_DATE"])
    pp, pe = first(cols["p"],["PATID"]), first(cols["p"],["ENCOUNTERID","ENCOUNTER_ID"])
    pod, pot = first(cols["p"],["RX_ORDER_DATE"]), first(cols["p"],["RX_ORDER_TIME"])
    pstart, pend = first(cols["p"],["RX_START_DATE","START_DATE"]), first(cols["p"],["RX_END_DATE","END_DATE"])
    pn, pc, pr = first(cols["p"],["RAW_RX_MED_NAME"]), first(cols["p"],["RXNORM_CUI"]), first(cols["p"],["RX_ROUTE"])
    dp, dd = first(cols["d"],["PATID"]), first(cols["d"],["DEATH_DATE"])
    legacy = root / "PCORnet" / "code" / "config" / "codes_antibiotics.py"
    inc, exc = parse_legacy_codes(legacy)
    incsql = ",".join(q(x) for x in sorted(inc)) or "''"; excsql = ",".join(q(x) for x in sorted(exc)) or "''"
    pts = ts("p", pod, pot)
    pname = f"lower(coalesce(cast(p.{qi(pn)} as varchar),''))"; pcode = f"trim(coalesce(cast(p.{qi(pc)} as varchar),''))"
    route = f"upper(trim(coalesce(cast(p.{qi(pr)} as varchar),'')))" if pr else "''"
    broad = f"(({pcode} in ({incsql}) or regexp_matches({pname},{q(BROAD_PATTERN)})) and {pcode} not in ({excsql}) and not regexp_matches({pname},{q(NON_SYSTEMIC_PATTERN)}) and {route} not in ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    nonbroad = f"(regexp_matches({pname},{q(NON_BROAD_PATTERN)}) and not regexp_matches({pname},{q(NON_SYSTEMIC_PATTERN)}) and {route} not in ('ORAL','RESPIRATORY_TRACT','INHALATION'))"

    con.execute(f"create temp table base as select distinct cast({qi(sp)} as varchar) patid,cast({qi(se)} as varchar) encounterid,try_cast({qi(ad)} as date) admit_date,try_cast({qi(dc)} as date) discharge_date from s")
    con.execute(f"create temp table anchors as select cast(p.{qi(pp)} as varchar) patid,cast(p.{qi(pe)} as varchar) encounterid,min({pts}) anchor_ts from p join base b on cast(p.{qi(pp)} as varchar)=b.patid and cast(p.{qi(pe)} as varchar)=b.encounterid where {broad} and {pts}>=cast(b.admit_date as timestamp) and {pts}<cast(b.admit_date as timestamp)+interval 24 hour group by 1,2")
    con.execute(f"create temp table deaths as select cast({qi(dp)} as varchar) patid,min(try_cast({qi(dd)} as date)) death_date from d group by 1")
    con.execute("create temp table cohort as select a.*,b.admit_date,b.discharge_date,d.death_date from anchors a join base b using(patid,encounterid) left join deaths d using(patid) where (b.discharge_date is null or b.discharge_date>cast(a.anchor_ts+interval 96 hour as date)) and (d.death_date is null or d.death_date>cast(a.anchor_ts+interval 96 hour as date))")
    n = int(con.execute("select count(*) from cohort").fetchone()[0])

    p_start_expr = f"coalesce(try_cast(p.{qi(pstart)} as date),try_cast(p.{qi(pod)} as date))" if pstart else f"try_cast(p.{qi(pod)} as date)"
    p_end_expr = f"coalesce(try_cast(p.{qi(pend)} as date),{p_start_expr})" if pend else p_start_expr
    con.execute(f"""create temp table exposure as
        select c.patid,c.encounterid,
        max(case when {broad} and {p_start_expr}<=cast(c.anchor_ts+interval 96 hour as date) and {p_end_expr}>=cast(c.anchor_ts+interval 72 hour as date) then 1 else 0 end) broad_72_96,
        max(case when {nonbroad} and {p_start_expr}<=cast(c.anchor_ts+interval 96 hour as date) and {p_end_expr}>=cast(c.anchor_ts+interval 72 hour as date) then 1 else 0 end) nonbroad_72_96
        from cohort c left join p on cast(p.{qi(pp)} as varchar)=c.patid and cast(p.{qi(pe)} as varchar)=c.encounterid group by 1,2""")

    df = con.execute("select c.patid,c.encounterid,c.anchor_ts,c.admit_date,case when e.broad_72_96=1 then 0 else 1 end A from cohort c join exposure e using(patid,encounterid)").fetchdf()

    # demographics
    dcols = cols["dem"]; dpat=first(dcols,["PATID"]); birth=first(dcols,["BIRTH_DATE"]); sex=first(dcols,["SEX"]); race=first(dcols,["RACE"])
    dem = con.execute(f"""select cast(x.{qi(dpat)} as varchar) patid,
        max(date_diff('year',try_cast(x.{qi(birth)} as date),current_date)) age_proxy,
        max(case when upper(trim(cast(x.{qi(sex)} as varchar))) in ('M','MALE','1') then 1 else 0 end) sex_male,
        max(case when upper(trim(cast(x.{qi(race)} as varchar))) in ('05','WHITE') then 1 else 0 end) race_white
        from dem x group by 1""").fetchdf()
    df = df.merge(dem,on="patid",how="left")
    # recompute age at admit, avoiding current-date age proxy
    births = con.execute(f"select cast({qi(dpat)} as varchar) patid,min(try_cast({qi(birth)} as date)) birth_date from dem group by 1").fetchdf()
    df = df.merge(births,on="patid",how="left"); df["age"] = ((pd.to_datetime(df["admit_date"])-pd.to_datetime(df["birth_date"])).dt.days/365.25); df.drop(columns=["age_proxy","birth_date"],inplace=True)

    # diagnosis proxies
    dxcols=cols["dx"]; dxp=first(dxcols,["PATID"]); dxe=first(dxcols,["ENCOUNTERID","ENCOUNTER_ID"]); dxc=first(dxcols,["DX","DIAGNOSIS_CODE"])
    norm=f"upper(regexp_replace(coalesce(cast(x.{qi(dxc)} as varchar),''),'[^A-Z0-9]','','g'))"
    pieces=[]
    for name,pat in DX_PATTERNS.items(): pieces.append(f"max(case when regexp_matches({norm},{q(pat)}) then 1 else 0 end) {name}")
    dxf=con.execute(f"select cast(x.{qi(dxp)} as varchar) patid,cast(x.{qi(dxe)} as varchar) encounterid,{','.join(pieces)} from dx x group by 1,2").fetchdf()
    dxf["comorb"]=(dxf[list(DX_PATTERNS)].max(axis=1)>0).astype(int)
    df=df.merge(dxf[["patid","encounterid","comorb","heart_failure","chronic_kidney"]],on=["patid","encounterid"],how="left")

    # helpers for row-level frozen trajectory features, kept local only
    def merge_sql_feature(sql, names):
        nonlocal df
        z=con.execute(sql).fetchdf(); df=df.merge(z,on=["patid","encounterid"],how="left")

    lcols=cols["lab"]; lp=first(lcols,["PATID"]); le=first(lcols,["ENCOUNTERID","ENCOUNTER_ID"]); lc=first(lcols,["LAB_LOINC"]); lv=first(lcols,["RESULT_NUM"]); ldate=first(lcols,["SPECIMEN_DATE"]); ltime=first(lcols,["SPECIMEN_TIME"])
    lts=ts("x",ldate,ltime,numeric_seconds=True)
    for concept,spec in LABS.items():
        codes=",".join(q(c) for c in spec["codes"]); lo,hi=spec["lo"],spec["hi"]
        val=f"try_cast(x.{qi(lv)} as double)"; filt=f"trim(cast(x.{qi(lc)} as varchar)) in ({codes}) and {val} between {lo} and {hi}"
        if spec["summary"]=="worst=min": early=f"min({val}) filter(where {lts}>=c.anchor_ts and {lts}<c.anchor_ts+interval 24 hour)"; late=f"min({val}) filter(where {lts}>=c.anchor_ts+interval 48 hour and {lts}<c.anchor_ts+interval 72 hour)"
        elif spec["summary"]=="worst=max": early=f"max({val}) filter(where {lts}>=c.anchor_ts and {lts}<c.anchor_ts+interval 24 hour)"; late=f"max({val}) filter(where {lts}>=c.anchor_ts+interval 48 hour and {lts}<c.anchor_ts+interval 72 hour)"
        else:
            # arg_max(value,time) gives the last measurement in the window
            early=f"arg_max({val},{lts}) filter(where {lts}>=c.anchor_ts and {lts}<c.anchor_ts+interval 24 hour)"; late=f"arg_max({val},{lts}) filter(where {lts}>=c.anchor_ts+interval 48 hour and {lts}<c.anchor_ts+interval 72 hour)"
        merge_sql_feature(f"select c.patid,c.encounterid,{early} {concept}_0_24h,{late} {concept}_48_72h from cohort c left join lab x on cast(x.{qi(lp)} as varchar)=c.patid and cast(x.{qi(le)} as varchar)=c.encounterid and {filt} group by 1,2",[concept])
        df[f"{concept}_delta"] = df[f"{concept}_48_72h"] - df[f"{concept}_0_24h"]

    ocols=cols["obs"]; op=first(ocols,["PATID"]); oe=first(ocols,["ENCOUNTERID","ENCOUNTER_ID"]); oc=first(ocols,["OBSCLIN_CODE"]); ov=first(ocols,["OBSCLIN_RESULT_NUM","RESULT_NUM"]); od=first(ocols,["OBSCLIN_START_DATE","OBSCLIN_DATE"]); ot=first(ocols,["OBSCLIN_START_TIME","OBSCLIN_TIME"]); ou=first(ocols,["OBSCLIN_RESULT_UNIT","RESULT_UNIT"])
    ots=ts("x",od,ot,numeric_seconds=True); oval=f"try_cast(x.{qi(ov)} as double)"
    for concept,spec in OBS.items():
        baseval=oval
        if concept=="temperature" and ou:
            baseval=f"case when upper(trim(cast(x.{qi(ou)} as varchar))) in ('F','FAHRENHEIT','[DEGF]') then ({oval}-32)*5/9 else {oval} end"
        filt=f"trim(cast(x.{qi(oc)} as varchar))={q(spec['code'])} and {baseval} between {spec['lo']} and {spec['hi']}"
        fun="min" if concept=="spo2" else "max"
        early=f"{fun}({baseval}) filter(where {ots}>=c.anchor_ts and {ots}<c.anchor_ts+interval 24 hour)"; late=f"{fun}({baseval}) filter(where {ots}>=c.anchor_ts+interval 48 hour and {ots}<c.anchor_ts+interval 72 hour)"
        merge_sql_feature(f"select c.patid,c.encounterid,{early} {concept}_0_24h,{late} {concept}_48_72h from cohort c left join obs x on cast(x.{qi(op)} as varchar)=c.patid and cast(x.{qi(oe)} as varchar)=c.encounterid and {filt} group by 1,2",[concept])
        df[f"{concept}_delta"] = df[f"{concept}_48_72h"] - df[f"{concept}_0_24h"]

    # derived MAP
    vcols=cols["vit"]; vp=first(vcols,["PATID"]); ve=first(vcols,["ENCOUNTERID","ENCOUNTER_ID"]); vd=first(vcols,["MEASURE_DATE","VITAL_DATE"]); vt=first(vcols,["MEASURE_TIME","VITAL_TIME"]); vs=first(vcols,["SYSTOLIC"]); vdia=first(vcols,["DIASTOLIC"])
    vts=ts("x",vd,vt,numeric_seconds=True); sb=f"try_cast(x.{qi(vs)} as double)"; db=f"try_cast(x.{qi(vdia)} as double)"; mv=f"({sb}+2*{db})/3.0"; vf=f"{sb} between 40 and 300 and {db} between 20 and 200 and {mv} between 20 and 200"
    merge_sql_feature(f"select c.patid,c.encounterid,min({mv}) filter(where {vts}>=c.anchor_ts and {vts}<c.anchor_ts+interval 24 hour) map_0_24h,min({mv}) filter(where {vts}>=c.anchor_ts+interval 48 hour and {vts}<c.anchor_ts+interval 72 hour) map_48_72h from cohort c left join vit x on cast(x.{qi(vp)} as varchar)=c.patid and cast(x.{qi(ve)} as varchar)=c.encounterid and {vf} group by 1,2",["map"])
    df["map_delta"]=df["map_48_72h"]-df["map_0_24h"]

    # vasopressor date-span rule
    mcols=cols["med"]; mp=first(mcols,["PATID"]); me=first(mcols,["ENCOUNTERID","ENCOUNTER_ID"]); msd=first(mcols,["MEDADMIN_START_DATE"]); med=first(mcols,["MEDADMIN_STOP_DATE"]); mn=first(mcols,["RAW_MEDADMIN_MED_NAME"])
    mname=f"lower(coalesce(cast(x.{qi(mn)} as varchar),''))"; incl=' or '.join(f"strpos({mname},{q(v)})>0" for v in VASO_TERMS); excl=' or '.join(f"strpos({mname},{q(v)})>0" for v in VASO_EXCLUDE); vaso=f"({incl}) and not ({excl})"
    mstart=f"cast(try_cast(x.{qi(msd)} as date) as timestamp)"; mstop=f"coalesce(cast(try_cast(x.{qi(med)} as date) as timestamp)+interval 23 hour+interval 59 minute+interval 59 second,{mstart}+interval 1 hour)" if med else f"{mstart}+interval 1 hour"
    merge_sql_feature(f"select c.patid,c.encounterid,max(case when {mstart}<c.anchor_ts+interval 24 hour and {mstop}>c.anchor_ts then 1 else 0 end) vasopressor_any_0_24h,max(case when {mstart}<c.anchor_ts+interval 72 hour and {mstop}>c.anchor_ts+interval 48 hour then 1 else 0 end) vasopressor_any_48_72h from cohort c left join med x on cast(x.{qi(mp)} as varchar)=c.patid and cast(x.{qi(me)} as varchar)=c.encounterid and {vaso} group by 1,2",["vaso"])
    df["vasopressor_stopped_before_72h"]=((df["vasopressor_any_0_24h"]==1)&(df["vasopressor_any_48_72h"]==0)).astype(int)

    binary = ["sex_male","race_white","comorb","heart_failure","chronic_kidney","vasopressor_any_0_24h","vasopressor_any_48_72h","vasopressor_stopped_before_72h"]
    continuous = ["age"]
    for concept in list(LABS)+list(OBS)+["map"]:
        continuous += [f"{concept}_0_24h",f"{concept}_48_72h",f"{concept}_delta"]
    features = binary + continuous

    # aggregate missingness before imputation
    miss=[]
    for col in features:
        for av in [0,1]:
            sub=df[df.A==av]; k=int(sub[col].isna().sum())
            miss.append({"variable":col,"A":av,"n":safe(len(sub)),"missing_n":safe(k) if k else 0,"missing_prop":float(k/len(sub)) if len(sub) else None})
    pd.DataFrame(miss).to_csv(args.output_dir/"missingness_by_exposure.csv",index=False)

    X=df[features].copy()
    for col in binary: X[col]=pd.to_numeric(X[col],errors="coerce").fillna(0).clip(0,1)
    for col in continuous:
        X[col]=pd.to_numeric(X[col],errors="coerce")
        medv=float(X[col].median()) if X[col].notna().any() else 0.0
        X[col]=X[col].fillna(medv)
        sd=float(X[col].std(ddof=0)); mu=float(X[col].mean())
        X[col]=(X[col]-mu)/(sd if sd>0 else 1.0); X[col]=X[col].clip(-8,8)
    y=df.A.astype(int).to_numpy()
    Xc=sm.add_constant(X,has_constant="add")
    fit_method="glm"
    try:
        model=sm.GLM(y,Xc,family=sm.families.Binomial()).fit(maxiter=200,disp=0)
        ps=np.asarray(model.predict(Xc),float)
        if not np.all(np.isfinite(ps)): raise ValueError("nonfinite PS")
    except Exception:
        fit_method="regularized_glm_alpha_0.001"
        model=sm.GLM(y,Xc,family=sm.families.Binomial()).fit_regularized(alpha=0.001,L1_wt=0.0,maxiter=200)
        ps=np.asarray(model.predict(Xc),float)
    ps=np.clip(ps,0.001,0.999)
    pa=float(y.mean()); sw=np.where(y==1,pa/ps,(1-pa)/(1-ps))

    bal=[]
    for col in features:
        xv=X[col].to_numpy(float)
        pre=smd(xv,y); post=smd(xv,y,sw)
        bal.append({"variable":col,"pre_smd":pre,"post_smd":post,"abs_pre_smd":abs(pre) if np.isfinite(pre) else None,"abs_post_smd":abs(post) if np.isfinite(post) else None})
    bdf=pd.DataFrame(bal).sort_values("abs_post_smd",ascending=False); bdf.to_csv(args.output_dir/"covariate_balance.csv",index=False)

    psrows=[]
    for av in [0,1]:
        z=ps[y==av]
        psrows.append({"A":av,"n":safe(len(z)),"min":float(np.min(z)),"p01":float(np.quantile(z,.01)),"p05":float(np.quantile(z,.05)),"median":float(np.quantile(z,.5)),"p95":float(np.quantile(z,.95)),"p99":float(np.quantile(z,.99)),"max":float(np.max(z))})
    pd.DataFrame(psrows).to_csv(args.output_dir/"ps_distribution.csv",index=False)

    wrows=[]
    for av in [0,1]:
        z=sw[y==av]
        wrows.append({"A":av,"n":safe(len(z)),"ess":ess(z),"mean_weight":float(np.mean(z)),"p95_weight":float(np.quantile(z,.95)),"p99_weight":float(np.quantile(z,.99)),"max_weight":float(np.max(z))})
    pd.DataFrame(wrows).to_csv(args.output_dir/"weight_diagnostics.csv",index=False)

    max_pre=float(bdf.abs_pre_smd.max()); max_post=float(bdf.abs_post_smd.max()); worst=str(bdf.iloc[0].variable)
    summary={
        "privacy_mode":"aggregate_only",
        "strict_cohort_n":safe(n),
        "deescalated_n":safe(int(y.sum())),
        "continued_n":safe(int((1-y).sum())),
        "exposure":"A=1 no broad PRESCRIBING overlap in day-3 72-96h window; A=0 any broad overlap, using frozen date-level prescribing rule",
        "model":"stabilized ATE IPTW; continuous median-imputed, standardized and clipped +/-8; binary fill 0; PS clipped [0.001,0.999]",
        "fit_method":fit_method,
        "max_abs_pre_smd":max_pre,
        "max_abs_post_smd":max_post,
        "worst_post_balance_variable":worst,
        "treated_ess":ess(sw[y==1]),
        "continued_ess":ess(sw[y==0]),
        "max_weight":float(np.max(sw)),
        "guardrail":"No outcomes or treatment effects computed. Do not tune covariates against outcomes.",
    }
    (args.output_dir/"summary.json").write_text(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()

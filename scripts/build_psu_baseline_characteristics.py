#!/usr/bin/env python3
"""Build aggregate Penn State baseline characteristics for manuscript reporting.

This is a descriptive publication-support task. It reconstructs the frozen strict
Penn State modified-replication cohort and its prespecified propensity covariates,
reproduces the frozen stabilized ATE weights for balance reporting, and exports
aggregate summaries only. No outcomes or treatment-effect estimates are computed.
"""
from __future__ import annotations

import argparse
import json
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
EXPECTED_N = 19841
EXPECTED_DEESC = 5346
EXPECTED_CONT = 14495
EXPECTED_MAX_POST_SMD = 0.022621815
EXPECTED_WORST_POST = "vasopressor_any_0_24h"

CHARACTERISTICS = [
    ("Demographics", "age", "Age, years", "continuous", 1),
    ("Demographics", "sex_male", "Male sex", "binary", 1),
    ("Demographics", "race_white", "White race", "binary", 1),
    ("Comorbidity", "comorb", "Any recorded comorbidity proxy", "binary", 1),
    ("Comorbidity", "heart_failure", "Heart failure", "binary", 1),
    ("Comorbidity", "chronic_kidney", "Chronic kidney disease", "binary", 1),
    ("Early severity", "lactate_0_24h", "Last lactate, 0-24 h, mmol/L", "continuous", 1),
    ("Early severity", "creatinine_0_24h", "Last creatinine, 0-24 h, mg/dL", "continuous", 2),
    ("Early severity", "wbc_0_24h", "Last WBC count, 0-24 h", "continuous", 1),
    ("Early severity", "heart_rate_0_24h", "Maximum heart rate, 0-24 h, beats/min", "continuous", 0),
    ("Early severity", "map_0_24h", "Minimum MAP, 0-24 h, mmHg", "continuous", 0),
    ("Early severity", "vasopressor_any_0_24h", "Vasopressor exposure, 0-24 h", "binary", 1),
    ("Near-decision status", "lactate_48_72h", "Last lactate, 48-72 h, mmol/L", "continuous", 1),
    ("Near-decision status", "creatinine_48_72h", "Last creatinine, 48-72 h, mg/dL", "continuous", 2),
    ("Near-decision status", "wbc_48_72h", "Last WBC count, 48-72 h", "continuous", 1),
    ("Near-decision status", "platelet_48_72h", "Worst platelet count, 48-72 h", "continuous", 0),
    ("Near-decision status", "bilirubin_total_48_72h", "Worst bilirubin, 48-72 h, mg/dL", "continuous", 1),
    ("Near-decision status", "heart_rate_48_72h", "Maximum heart rate, 48-72 h, beats/min", "continuous", 0),
    ("Near-decision status", "resp_rate_48_72h", "Maximum respiratory rate, 48-72 h", "continuous", 0),
    ("Near-decision status", "spo2_48_72h", "Minimum oxygen saturation, 48-72 h, %", "continuous", 0),
    ("Near-decision status", "temperature_48_72h", "Maximum temperature, 48-72 h, degC", "continuous", 1),
    ("Near-decision status", "map_48_72h", "Minimum MAP, 48-72 h, mmHg", "continuous", 0),
    ("Near-decision status", "vasopressor_any_48_72h", "Vasopressor exposure, 48-72 h", "binary", 1),
    ("Clinical trajectory", "vasopressor_stopped_before_72h", "Vasopressor stopped before 72 h", "binary", 1),
]


def weighted_mean(x, w):
    x = np.asarray(x, float)
    w = np.asarray(w, float)
    m = np.isfinite(x) & np.isfinite(w)
    if not m.any() or w[m].sum() <= 0:
        return np.nan
    return float(np.sum(x[m] * w[m]) / np.sum(w[m]))


def weighted_var(x, w):
    x = np.asarray(x, float)
    w = np.asarray(w, float)
    m = np.isfinite(x) & np.isfinite(w)
    if not m.any() or w[m].sum() <= 0:
        return np.nan
    mu = weighted_mean(x[m], w[m])
    return float(np.sum(w[m] * (x[m] - mu) ** 2) / np.sum(w[m]))


def smd(x, a, w=None):
    x = np.asarray(x, float)
    a = np.asarray(a, int)
    w = np.ones(len(x), float) if w is None else np.asarray(w, float)
    m1, m0 = a == 1, a == 0
    mu1, mu0 = weighted_mean(x[m1], w[m1]), weighted_mean(x[m0], w[m0])
    v1, v0 = weighted_var(x[m1], w[m1]), weighted_var(x[m0], w[m0])
    den = np.sqrt((v1 + v0) / 2.0)
    if not np.isfinite(den) or den == 0:
        return 0.0 if np.isfinite(mu1) and np.isfinite(mu0) and mu1 == mu0 else np.nan
    return float((mu1 - mu0) / den)


def _continuous_summary(x: pd.Series) -> dict:
    z = pd.to_numeric(x, errors="coerce")
    obs = z.dropna()
    return {
        "nonmissing_n": int(obs.size),
        "missing_n": int(z.isna().sum()),
        "median": float(obs.median()) if obs.size else np.nan,
        "q1": float(obs.quantile(0.25)) if obs.size else np.nan,
        "q3": float(obs.quantile(0.75)) if obs.size else np.nan,
    }


def _binary_summary(x: pd.Series) -> dict:
    raw = pd.to_numeric(x, errors="coerce")
    z = raw.fillna(0).clip(0, 1)
    return {
        "nonmissing_n": int(raw.notna().sum()),
        "missing_n": int(raw.isna().sum()),
        "positive_n": int(z.sum()),
        "positive_percent": float(100 * z.mean()) if len(z) else np.nan,
    }


def _format_continuous(s: dict, digits: int) -> str:
    if not np.isfinite(s["median"]):
        return "NA"
    fmt = f"{{:.{digits}f}}"
    return f"{fmt.format(s['median'])} [{fmt.format(s['q1'])}, {fmt.format(s['q3'])}]"


def _format_binary(s: dict) -> str:
    n = int(s["positive_n"])
    if 0 < n < MIN_CELL:
        return "<11"
    return f"{n:,} ({s['positive_percent']:.1f}%)"


def reconstruct(data_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, dict]:
    root = data_root
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
    inc, exc = parse_legacy_codes(root / "PCORnet" / "code" / "config" / "codes_antibiotics.py")
    incsql = ",".join(q(x) for x in sorted(inc)) or "''"
    excsql = ",".join(q(x) for x in sorted(exc)) or "''"
    pts = ts("p", pod, pot)
    pname = f"lower(coalesce(cast(p.{qi(pn)} as varchar),''))"
    pcode = f"trim(coalesce(cast(p.{qi(pc)} as varchar),''))"
    route = f"upper(trim(coalesce(cast(p.{qi(pr)} as varchar),'')))" if pr else "''"
    broad = f"(({pcode} in ({incsql}) or regexp_matches({pname},{q(BROAD_PATTERN)})) and {pcode} not in ({excsql}) and not regexp_matches({pname},{q(NON_SYSTEMIC_PATTERN)}) and {route} not in ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    nonbroad = f"(regexp_matches({pname},{q(NON_BROAD_PATTERN)}) and not regexp_matches({pname},{q(NON_SYSTEMIC_PATTERN)}) and {route} not in ('ORAL','RESPIRATORY_TRACT','INHALATION'))"

    con.execute(f"create temp table base as select distinct cast({qi(sp)} as varchar) patid,cast({qi(se)} as varchar) encounterid,try_cast({qi(ad)} as date) admit_date,try_cast({qi(dc)} as date) discharge_date from s")
    con.execute(f"create temp table anchors as select cast(p.{qi(pp)} as varchar) patid,cast(p.{qi(pe)} as varchar) encounterid,min({pts}) anchor_ts from p join base b on cast(p.{qi(pp)} as varchar)=b.patid and cast(p.{qi(pe)} as varchar)=b.encounterid where {broad} and {pts}>=cast(b.admit_date as timestamp) and {pts}<cast(b.admit_date as timestamp)+interval 24 hour group by 1,2")
    con.execute(f"create temp table deaths as select cast({qi(dp)} as varchar) patid,min(try_cast({qi(dd)} as date)) death_date from d group by 1")
    con.execute("create temp table cohort as select a.*,b.admit_date,b.discharge_date,d.death_date from anchors a join base b using(patid,encounterid) left join deaths d using(patid) where (b.discharge_date is null or b.discharge_date>cast(a.anchor_ts+interval 96 hour as date)) and (d.death_date is null or d.death_date>cast(a.anchor_ts+interval 96 hour as date))")
    p_start_expr = f"coalesce(try_cast(p.{qi(pstart)} as date),try_cast(p.{qi(pod)} as date))" if pstart else f"try_cast(p.{qi(pod)} as date)"
    p_end_expr = f"coalesce(try_cast(p.{qi(pend)} as date),{p_start_expr})" if pend else p_start_expr
    con.execute(f"""create temp table exposure as
        select c.patid,c.encounterid,
        max(case when {broad} and {p_start_expr}<=cast(c.anchor_ts+interval 96 hour as date) and {p_end_expr}>=cast(c.anchor_ts+interval 72 hour as date) then 1 else 0 end) broad_72_96,
        max(case when {nonbroad} and {p_start_expr}<=cast(c.anchor_ts+interval 96 hour as date) and {p_end_expr}>=cast(c.anchor_ts+interval 72 hour as date) then 1 else 0 end) nonbroad_72_96
        from cohort c left join p on cast(p.{qi(pp)} as varchar)=c.patid and cast(p.{qi(pe)} as varchar)=c.encounterid group by 1,2""")
    df = con.execute("select c.patid,c.encounterid,c.anchor_ts,c.admit_date,case when e.broad_72_96=1 then 0 else 1 end A from cohort c join exposure e using(patid,encounterid)").fetchdf()

    dcols = cols["dem"]
    dpat = first(dcols, ["PATID"])
    birth = first(dcols, ["BIRTH_DATE"])
    sex = first(dcols, ["SEX"])
    race = first(dcols, ["RACE"])
    dem = con.execute(f"""select cast(x.{qi(dpat)} as varchar) patid,
        max(case when upper(trim(cast(x.{qi(sex)} as varchar))) in ('M','MALE','1') then 1 else 0 end) sex_male,
        max(case when upper(trim(cast(x.{qi(race)} as varchar))) in ('05','WHITE') then 1 else 0 end) race_white
        from dem x group by 1""").fetchdf()
    births = con.execute(f"select cast({qi(dpat)} as varchar) patid,min(try_cast({qi(birth)} as date)) birth_date from dem group by 1").fetchdf()
    df = df.merge(dem, on="patid", how="left").merge(births, on="patid", how="left")
    df["age"] = (pd.to_datetime(df["admit_date"]) - pd.to_datetime(df["birth_date"])).dt.days / 365.25
    df.drop(columns=["birth_date"], inplace=True)

    dxcols = cols["dx"]
    dxp = first(dxcols, ["PATID"])
    dxe = first(dxcols, ["ENCOUNTERID", "ENCOUNTER_ID"])
    dxc = first(dxcols, ["DX", "DIAGNOSIS_CODE"])
    norm = f"upper(regexp_replace(coalesce(cast(x.{qi(dxc)} as varchar),''),'[^A-Z0-9]','','g'))"
    pieces = [f"max(case when regexp_matches({norm},{q(pat)}) then 1 else 0 end) {name}" for name, pat in DX_PATTERNS.items()]
    dxf = con.execute(f"select cast(x.{qi(dxp)} as varchar) patid,cast(x.{qi(dxe)} as varchar) encounterid,{','.join(pieces)} from dx x group by 1,2").fetchdf()
    dxf["comorb"] = (dxf[list(DX_PATTERNS)].max(axis=1) > 0).astype(int)
    df = df.merge(dxf[["patid", "encounterid", "comorb", "heart_failure", "chronic_kidney"]], on=["patid", "encounterid"], how="left")

    def merge_sql_feature(sql):
        nonlocal df
        df = df.merge(con.execute(sql).fetchdf(), on=["patid", "encounterid"], how="left")

    lcols = cols["lab"]
    lp = first(lcols, ["PATID"])
    le = first(lcols, ["ENCOUNTERID", "ENCOUNTER_ID"])
    lc = first(lcols, ["LAB_LOINC"])
    lv = first(lcols, ["RESULT_NUM"])
    ldate = first(lcols, ["SPECIMEN_DATE"])
    ltime = first(lcols, ["SPECIMEN_TIME"])
    lts = ts("x", ldate, ltime, numeric_seconds=True)
    for concept, spec in LABS.items():
        codes = ",".join(q(c) for c in spec["codes"])
        val = f"try_cast(x.{qi(lv)} as double)"
        filt = f"trim(cast(x.{qi(lc)} as varchar)) in ({codes}) and {val} between {spec['lo']} and {spec['hi']}"
        if spec["summary"] == "worst=min":
            early = f"min({val}) filter(where {lts}>=c.anchor_ts and {lts}<c.anchor_ts+interval 24 hour)"
            late = f"min({val}) filter(where {lts}>=c.anchor_ts+interval 48 hour and {lts}<c.anchor_ts+interval 72 hour)"
        elif spec["summary"] == "worst=max":
            early = f"max({val}) filter(where {lts}>=c.anchor_ts and {lts}<c.anchor_ts+interval 24 hour)"
            late = f"max({val}) filter(where {lts}>=c.anchor_ts+interval 48 hour and {lts}<c.anchor_ts+interval 72 hour)"
        else:
            early = f"arg_max({val},{lts}) filter(where {lts}>=c.anchor_ts and {lts}<c.anchor_ts+interval 24 hour)"
            late = f"arg_max({val},{lts}) filter(where {lts}>=c.anchor_ts+interval 48 hour and {lts}<c.anchor_ts+interval 72 hour)"
        merge_sql_feature(f"select c.patid,c.encounterid,{early} {concept}_0_24h,{late} {concept}_48_72h from cohort c left join lab x on cast(x.{qi(lp)} as varchar)=c.patid and cast(x.{qi(le)} as varchar)=c.encounterid and {filt} group by 1,2")
        df[f"{concept}_delta"] = df[f"{concept}_48_72h"] - df[f"{concept}_0_24h"]

    ocols = cols["obs"]
    op = first(ocols, ["PATID"])
    oe = first(ocols, ["ENCOUNTERID", "ENCOUNTER_ID"])
    oc = first(ocols, ["OBSCLIN_CODE"])
    ov = first(ocols, ["OBSCLIN_RESULT_NUM", "RESULT_NUM"])
    od = first(ocols, ["OBSCLIN_START_DATE", "OBSCLIN_DATE"])
    ot = first(ocols, ["OBSCLIN_START_TIME", "OBSCLIN_TIME"])
    ou = first(ocols, ["OBSCLIN_RESULT_UNIT", "RESULT_UNIT"])
    ots = ts("x", od, ot, numeric_seconds=True)
    oval = f"try_cast(x.{qi(ov)} as double)"
    for concept, spec in OBS.items():
        baseval = oval
        if concept == "temperature" and ou:
            baseval = f"case when upper(trim(cast(x.{qi(ou)} as varchar))) in ('F','FAHRENHEIT','[DEGF]') then ({oval}-32)*5/9 else {oval} end"
        filt = f"trim(cast(x.{qi(oc)} as varchar))={q(spec['code'])} and {baseval} between {spec['lo']} and {spec['hi']}"
        fun = "min" if concept == "spo2" else "max"
        early = f"{fun}({baseval}) filter(where {ots}>=c.anchor_ts and {ots}<c.anchor_ts+interval 24 hour)"
        late = f"{fun}({baseval}) filter(where {ots}>=c.anchor_ts+interval 48 hour and {ots}<c.anchor_ts+interval 72 hour)"
        merge_sql_feature(f"select c.patid,c.encounterid,{early} {concept}_0_24h,{late} {concept}_48_72h from cohort c left join obs x on cast(x.{qi(op)} as varchar)=c.patid and cast(x.{qi(oe)} as varchar)=c.encounterid and {filt} group by 1,2")
        df[f"{concept}_delta"] = df[f"{concept}_48_72h"] - df[f"{concept}_0_24h"]

    vcols = cols["vit"]
    vp = first(vcols, ["PATID"])
    ve = first(vcols, ["ENCOUNTERID", "ENCOUNTER_ID"])
    vd = first(vcols, ["MEASURE_DATE", "VITAL_DATE"])
    vt = first(vcols, ["MEASURE_TIME", "VITAL_TIME"])
    vs = first(vcols, ["SYSTOLIC"])
    vdia = first(vcols, ["DIASTOLIC"])
    vts = ts("x", vd, vt, numeric_seconds=True)
    sb = f"try_cast(x.{qi(vs)} as double)"
    db = f"try_cast(x.{qi(vdia)} as double)"
    mv = f"({sb}+2*{db})/3.0"
    vf = f"{sb} between 40 and 300 and {db} between 20 and 200 and {mv} between 20 and 200"
    merge_sql_feature(f"select c.patid,c.encounterid,min({mv}) filter(where {vts}>=c.anchor_ts and {vts}<c.anchor_ts+interval 24 hour) map_0_24h,min({mv}) filter(where {vts}>=c.anchor_ts+interval 48 hour and {vts}<c.anchor_ts+interval 72 hour) map_48_72h from cohort c left join vit x on cast(x.{qi(vp)} as varchar)=c.patid and cast(x.{qi(ve)} as varchar)=c.encounterid and {vf} group by 1,2")
    df["map_delta"] = df["map_48_72h"] - df["map_0_24h"]

    mcols = cols["med"]
    mp = first(mcols, ["PATID"])
    me = first(mcols, ["ENCOUNTERID", "ENCOUNTER_ID"])
    msd = first(mcols, ["MEDADMIN_START_DATE"])
    med = first(mcols, ["MEDADMIN_STOP_DATE"])
    mn = first(mcols, ["RAW_MEDADMIN_MED_NAME"])
    mname = f"lower(coalesce(cast(x.{qi(mn)} as varchar),''))"
    incl = " or ".join(f"strpos({mname},{q(v)})>0" for v in VASO_TERMS)
    excl = " or ".join(f"strpos({mname},{q(v)})>0" for v in VASO_EXCLUDE)
    vaso = f"({incl}) and not ({excl})"
    mstart = f"cast(try_cast(x.{qi(msd)} as date) as timestamp)"
    mstop = f"coalesce(cast(try_cast(x.{qi(med)} as date) as timestamp)+interval 23 hour+interval 59 minute+interval 59 second,{mstart}+interval 1 hour)" if med else f"{mstart}+interval 1 hour"
    merge_sql_feature(f"select c.patid,c.encounterid,max(case when {mstart}<c.anchor_ts+interval 24 hour and {mstop}>c.anchor_ts then 1 else 0 end) vasopressor_any_0_24h,max(case when {mstart}<c.anchor_ts+interval 72 hour and {mstop}>c.anchor_ts+interval 48 hour then 1 else 0 end) vasopressor_any_48_72h from cohort c left join med x on cast(x.{qi(mp)} as varchar)=c.patid and cast(x.{qi(me)} as varchar)=c.encounterid and {vaso} group by 1,2")
    df["vasopressor_stopped_before_72h"] = ((df["vasopressor_any_0_24h"] == 1) & (df["vasopressor_any_48_72h"] == 0)).astype(int)

    binary = ["sex_male", "race_white", "comorb", "heart_failure", "chronic_kidney", "vasopressor_any_0_24h", "vasopressor_any_48_72h", "vasopressor_stopped_before_72h"]
    continuous = ["age"]
    for concept in list(LABS) + list(OBS) + ["map"]:
        continuous += [f"{concept}_0_24h", f"{concept}_48_72h", f"{concept}_delta"]
    features = binary + continuous
    X = df[features].copy()
    for col in binary:
        X[col] = pd.to_numeric(X[col], errors="coerce").fillna(0).clip(0, 1)
    for col in continuous:
        X[col] = pd.to_numeric(X[col], errors="coerce")
        medv = float(X[col].median()) if X[col].notna().any() else 0.0
        X[col] = X[col].fillna(medv)
        sd = float(X[col].std(ddof=0))
        mu = float(X[col].mean())
        X[col] = ((X[col] - mu) / (sd if sd > 0 else 1.0)).clip(-8, 8)
    y = df["A"].astype(int).to_numpy()
    Xc = sm.add_constant(X, has_constant="add")
    try:
        model = sm.GLM(y, Xc, family=sm.families.Binomial()).fit(maxiter=200, disp=0)
        ps = np.asarray(model.predict(Xc), float)
        if not np.all(np.isfinite(ps)):
            raise ValueError("nonfinite PS")
        fit_method = "glm"
    except Exception:
        model = sm.GLM(y, Xc, family=sm.families.Binomial()).fit_regularized(alpha=0.001, L1_wt=0.0, maxiter=200)
        ps = np.asarray(model.predict(Xc), float)
        fit_method = "regularized_glm_alpha_0.001"
    ps = np.clip(ps, 0.001, 0.999)
    pa = float(y.mean())
    sw = np.where(y == 1, pa / ps, (1 - pa) / (1 - ps))
    bal = []
    for col in features:
        xv = X[col].to_numpy(float)
        pre = smd(xv, y)
        post = smd(xv, y, sw)
        bal.append({
            "variable": col,
            "abs_pre_smd": abs(pre) if np.isfinite(pre) else np.nan,
            "abs_post_smd": abs(post) if np.isfinite(post) else np.nan,
        })
    bdf = pd.DataFrame(bal).sort_values("abs_post_smd", ascending=False)
    meta = {
        "fit_method": fit_method,
        "max_post": float(bdf["abs_post_smd"].max()),
        "worst_post": str(bdf.iloc[0]["variable"]),
    }
    return df, bdf, sw, meta


def build_table(df: pd.DataFrame, balance: pd.DataFrame, meta: dict) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    n1 = int((df["A"] == 1).sum())
    n0 = int((df["A"] == 0).sum())
    if (len(df), n1, n0) != (EXPECTED_N, EXPECTED_DEESC, EXPECTED_CONT):
        raise RuntimeError(f"Frozen PSU cohort parity failed: n={len(df)}, A=1 {n1}, A=0 {n0}")
    if abs(meta["max_post"] - EXPECTED_MAX_POST_SMD) > 5e-4 or meta["worst_post"] != EXPECTED_WORST_POST:
        raise RuntimeError(f"Frozen PSU balance parity failed: max post-SMD {meta['max_post']:.6f}, worst {meta['worst_post']}")
    b = balance.set_index("variable")
    formatted = []
    detailed = []
    for section, var, label, kind, digits in CHARACTERISTICS:
        if var not in df.columns:
            raise RuntimeError(f"Missing frozen PSU descriptive variable: {var}")
        s1 = _binary_summary(df.loc[df.A == 1, var]) if kind == "binary" else _continuous_summary(df.loc[df.A == 1, var])
        s0 = _binary_summary(df.loc[df.A == 0, var]) if kind == "binary" else _continuous_summary(df.loc[df.A == 0, var])
        d1 = _format_binary(s1) if kind == "binary" else _format_continuous(s1, digits)
        d0 = _format_binary(s0) if kind == "binary" else _format_continuous(s0, digits)
        pre = float(b.loc[var, "abs_pre_smd"]) if var in b.index else np.nan
        post = float(b.loc[var, "abs_post_smd"]) if var in b.index else np.nan
        formatted.append({
            "section": section,
            "characteristic": label,
            f"De-escalation or stopping (n={EXPECTED_DEESC:,})": d1,
            f"Continued broad-spectrum (n={EXPECTED_CONT:,})": d0,
            "Absolute SMD before weighting": f"{pre:.3f}",
            "Absolute SMD after weighting": f"{post:.3f}",
        })
        detailed.append({
            "section": section,
            "variable": var,
            "characteristic": label,
            "type": kind,
            "deescalated_or_stopped_display": d1,
            "continued_broad_display": d0,
            "deescalated_missing_n": s1["missing_n"] if s1["missing_n"] == 0 or s1["missing_n"] >= MIN_CELL else None,
            "continued_missing_n": s0["missing_n"] if s0["missing_n"] == 0 or s0["missing_n"] >= MIN_CELL else None,
            "absolute_smd_before_weighting": pre,
            "absolute_smd_after_weighting": post,
        })
    metadata = {
        "purpose": "Aggregate Penn State baseline characteristics for manuscript/ESM reporting.",
        "cohort_n": EXPECTED_N,
        "deescalation_or_stopping_n": EXPECTED_DEESC,
        "continued_broad_n": EXPECTED_CONT,
        "summary_convention": "Continuous variables: median [IQR] from observed values; categorical variables: n (%) with missing binary values coded as zero to match the frozen PSU PS convention.",
        "smd_convention": "Absolute SMDs before and after the frozen PSU stabilized ATE IPTW; continuous variables use the frozen median-imputation/standardization rule for balance calculations.",
        "max_post_weighting_absolute_smd": meta["max_post"],
        "worst_post_weighting_variable": meta["worst_post"],
        "fit_method": meta["fit_method"],
        "modified_replication_note": "PSU lacks validated exact ICU-type, urine-output, GCS, FiO2, and MIMIC-equivalent microbiology-intensity terms; these are intentionally not shown as if harmonized.",
        "outcomes_analyzed": False,
        "row_level_artifacts_exported": False,
    }
    return pd.DataFrame(detailed), pd.DataFrame(formatted), metadata


def write_markdown(formatted: pd.DataFrame, metadata: dict, path: Path) -> None:
    g1 = f"De-escalation or stopping (n={EXPECTED_DEESC:,})"
    g0 = f"Continued broad-spectrum (n={EXPECTED_CONT:,})"
    lines = [
        "# Penn State baseline characteristics",
        "",
        "Continuous variables are median [IQR]; categorical variables are n (%). No hypothesis-test p-values are shown.",
        "",
        f"| Characteristic | {g1} | {g0} | Absolute SMD before | Absolute SMD after |",
        "|---|---:|---:|---:|---:|",
    ]
    last = None
    for _, row in formatted.iterrows():
        if row["section"] != last:
            lines.append(f"| **{row['section']}** |  |  |  |  |")
            last = row["section"]
        lines.append(f"| {row['characteristic']} | {row[g1]} | {row[g0]} | {row['Absolute SMD before weighting']} | {row['Absolute SMD after weighting']} |")
    lines += [
        "",
        "Notes: Penn State is a modified external replication. Variables unavailable or not faithfully harmonizable in the frozen PSU model are not substituted with post hoc proxies. Binary missing values are coded as zero in the frozen PSU propensity implementation; continuous descriptive summaries use observed values, while SMD calculations use the frozen median-imputation and standardization convention.",
        "",
        f"Primary PSU balance check: maximum post-weighting absolute SMD {metadata['max_post_weighting_absolute_smd']:.3f} ({metadata['worst_post_weighting_variable']}).",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    df, balance, _, fitmeta = reconstruct(args.data_root)
    detailed, formatted, metadata = build_table(df, balance, fitmeta)
    detailed.to_csv(args.output_dir / "psu_baseline_characteristics_detailed.csv", index=False)
    formatted.to_csv(args.output_dir / "psu_baseline_characteristics_formatted.csv", index=False)
    write_markdown(formatted, metadata, args.output_dir / "psu_baseline_characteristics.md")
    (args.output_dir / "psu_baseline_characteristics_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

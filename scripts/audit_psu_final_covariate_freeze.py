#!/usr/bin/env python3
"""Final aggregate-only PSU covariate freeze before propensity-score diagnostics.

This audit freezes the modified PSU external-replication covariate definitions using:
- the strict 96-hour hospital-clock cohort;
- exact core-lab LOINCs with SPECIMEN_TIME decoded as seconds since midnight;
- exact OBS_CLIN LOINCs for HR/RR/SpO2/temperature;
- MAP derived from paired VITAL systolic/diastolic values when available;
- conservative MED_ADMIN vasopressor agents with interval-overlap logic; and
- code-based chronic comorbidity proxies aligned to the MIMIC concepts.

Only aggregate counts and definition tables are exported. No identifiers, patient rows,
raw result values, propensity scores, or treatment effects are exported.
"""
from __future__ import annotations

import argparse
import ast
import json
import re
from pathlib import Path

import duckdb
import pandas as pd

MIN_CELL = 11
BROAD_PATTERN = (
    "vancomycin|vancocin|piperacillin|tazobactam|zosyn|cefepime|ceftazidime|fortaz|"
    "meropenem|merrem|imipenem|primaxin|aztreonam|azactam|linezolid|zyvox|"
    "daptomycin|cubicin|ceftolozane|zerbaxa|avibactam|avycaz"
)
NON_SYSTEMIC_PATTERN = "cayston|inhal|nebul|tablet|capsule|oral solution|oral suspension|by mouth|\\bpo\\b"

LABS = {
    "lactate": {"codes": ["2524-7", "32693-4", "19239-3"], "lo": 0.1, "hi": 30.0, "summary": "last"},
    "creatinine": {"codes": ["2160-0", "38483-4"], "lo": 0.1, "hi": 30.0, "summary": "last"},
    "wbc": {"codes": ["6690-2"], "lo": 0.1, "hi": 200.0, "summary": "last"},
    "platelet": {"codes": ["777-3"], "lo": 1.0, "hi": 2000.0, "summary": "worst=min"},
    "bilirubin_total": {"codes": ["1975-2"], "lo": 0.01, "hi": 80.0, "summary": "worst=max"},
}
OBS = {
    "heart_rate": {"code": "8867-4", "lo": 20.0, "hi": 250.0, "summary": "max"},
    "resp_rate": {"code": "9279-1", "lo": 1.0, "hi": 90.0, "summary": "max"},
    "spo2": {"code": "59408-5", "lo": 20.0, "hi": 100.0, "summary": "min"},
    "temperature": {"code": "8310-5", "lo": 25.0, "hi": 45.0, "summary": "max; convert F to C"},
}
VASO_TERMS = ["norepinephrine", "levophed", "phenylephrine", "vasopressin", "epinephrine", "dopamine"]
VASO_EXCLUDE = ["racepinephrine", "nasal", "ophthalm", "lidocaine", "topical"]

# Modified external-replication chronic-condition proxies. Codes are normalized by
# upper-casing and removing punctuation before matching.
DX_PATTERNS = {
    "heart_failure": r"^(I50|428)",
    "chronic_kidney": r"^(N18|585)",
    "diabetes": r"^(E1[0-4]|250)",
    "copd": r"^(J44|491|492|496)",
    "chronic_liver": r"^(K70|K72|K74|571)",
    "malignancy": r"^(C[0-9]{2}|14[0-9]|15[0-9]|16[0-9]|17[0-9]|18[0-9]|19[0-9]|20[0-8])",
}


def q(v: str) -> str:
    return "'" + str(v).replace("'", "''") + "'"


def qi(v: str) -> str:
    return '"' + str(v).replace('"', '""') + '"'


def safe(n):
    if n is None:
        return None
    n = int(n)
    return n if n >= MIN_CELL else None


def find(root: Path, stems: list[str], required: bool = True) -> Path | None:
    for stem in stems:
        candidates = sorted(root.glob(f"PCORnet/parquet/**/{stem}*.parquet"))
        if not candidates:
            candidates = sorted(root.glob(f"PCORnet/parquet/**/{stem.replace('_', '')}*.parquet"))
        if candidates:
            exact = [p for p in candidates if p.stem.lower() == stem.lower()]
            return exact[0] if exact else max(candidates, key=lambda p: p.stat().st_size)
    if required:
        raise FileNotFoundError(stems)
    return None


def first(cols: set[str], names: list[str]) -> str | None:
    lut = {c.upper(): c for c in cols}
    for name in names:
        if name.upper() in lut:
            return lut[name.upper()]
    return None


def parse_legacy_codes(path: Path) -> tuple[set[str], set[str]]:
    tree = ast.parse(path.read_text(errors="ignore"))
    include, exclude = set(), set()
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        names = [t.id.lower() for t in node.targets if isinstance(t, ast.Name)]
        try:
            value = ast.literal_eval(node.value)
        except Exception:
            continue
        for name in names:
            if "exclude" in name and isinstance(value, (set, list, tuple)):
                exclude |= {str(x).strip() for x in value if re.fullmatch(r"\d{4,9}", str(x).strip())}
            elif isinstance(value, (list, tuple)):
                for item in value:
                    if isinstance(item, (list, tuple)) and item and re.fullmatch(r"\d{4,9}", str(item[0]).strip()):
                        include.add(str(item[0]).strip())
    return include - exclude, exclude


def ts(alias: str, date_col: str | None, time_col: str | None, numeric_seconds: bool = False) -> str | None:
    if not date_col:
        return None
    d = f"try_cast({alias}.{qi(date_col)} as date)"
    if not time_col:
        return f"cast({d} as timestamp)"
    raw = f"trim(cast({alias}.{qi(time_col)} as varchar))"
    num = f"try_cast({raw} as double)"
    if numeric_seconds:
        return (
            f"case when {d} is null then null "
            f"when {num} between 0 and 86399 then cast({d} as timestamp)+({num}*interval 1 second) "
            f"when strpos({raw},':')>0 then try_cast(cast({alias}.{qi(date_col)} as varchar)||' '||{raw} as timestamp) "
            f"else null end"
        )
    return (
        f"case when {d} is null then null "
        f"when strpos({raw},':')>0 then try_cast(cast({alias}.{qi(date_col)} as varchar)||' '||{raw} as timestamp) "
        f"when {num} between 0 and 86399 then cast({d} as timestamp)+({num}*interval 1 second) "
        f"else try_cast(cast({alias}.{qi(date_col)} as varchar)||' '||{raw} as timestamp) end"
    )


def main() -> None:
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
        "dem": find(root, ["sepsis_demographic", "demographic"], required=False),
        "dx": find(root, ["sepsis_diagnosis", "diagnosis"], required=False),
        "lab": find(root, ["lab_reduced", "lab_result_cm"], required=False),
        "obs": find(root, ["obs_clin"], required=False),
        "vit": find(root, ["sepsis_vital", "vital"], required=False),
        "med": find(root, ["med_admin"], required=False),
    }
    con = duckdb.connect()
    con.execute("pragma threads=4")
    for key, path in paths.items():
        if path is not None:
            con.execute(f"create view {key} as select * from read_parquet({q(str(path))})")
    cols = {k: set(con.execute(f"describe {k}").fetchdf()["column_name"].astype(str)) for k, p in paths.items() if p is not None}

    # Strict modified PSU cohort, identical to the established feasibility audits.
    sp, se = first(cols["s"], ["PATID"]), first(cols["s"], ["ENCOUNTERID", "ENCOUNTER_ID"])
    ad, dc = first(cols["s"], ["ADMIT_DATE"]), first(cols["s"], ["DISCHARGE_DATE"])
    pp, pe = first(cols["p"], ["PATID"]), first(cols["p"], ["ENCOUNTERID", "ENCOUNTER_ID"])
    pod, pot = first(cols["p"], ["RX_ORDER_DATE"]), first(cols["p"], ["RX_ORDER_TIME"])
    pn, pc, pr = first(cols["p"], ["RAW_RX_MED_NAME"]), first(cols["p"], ["RXNORM_CUI"]), first(cols["p"], ["RX_ROUTE"])
    dp, dd = first(cols["d"], ["PATID"]), first(cols["d"], ["DEATH_DATE"])
    if any(x is None for x in [sp, se, ad, dc, pp, pe, pod, pn, pc, dp, dd]):
        raise RuntimeError("Required strict-cohort fields missing")
    legacy = root / "PCORnet" / "code" / "config" / "codes_antibiotics.py"
    inc, exc = parse_legacy_codes(legacy)
    incsql = ",".join(q(x) for x in sorted(inc)) or "''"
    excsql = ",".join(q(x) for x in sorted(exc)) or "''"
    pts = ts("p", pod, pot)
    name = f"lower(coalesce(cast(p.{qi(pn)} as varchar),''))"
    code = f"trim(coalesce(cast(p.{qi(pc)} as varchar),''))"
    route = f"upper(trim(coalesce(cast(p.{qi(pr)} as varchar),'')))" if pr else "''"
    broad = (
        f"(({code} in ({incsql}) or regexp_matches({name},{q(BROAD_PATTERN)})) and {code} not in ({excsql}) "
        f"and not regexp_matches({name},{q(NON_SYSTEMIC_PATTERN)}) and {route} not in ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    )
    con.execute(f"create temp table base as select distinct cast({qi(sp)} as varchar) patid,cast({qi(se)} as varchar) encounterid,try_cast({qi(ad)} as date) admit_date,try_cast({qi(dc)} as date) discharge_date from s")
    con.execute(f"create temp table anchors as select cast(p.{qi(pp)} as varchar) patid,cast(p.{qi(pe)} as varchar) encounterid,min({pts}) anchor_ts from p join base b on cast(p.{qi(pp)} as varchar)=b.patid and cast(p.{qi(pe)} as varchar)=b.encounterid where {broad} and {pts}>=cast(b.admit_date as timestamp) and {pts}<cast(b.admit_date as timestamp)+interval 24 hour group by 1,2")
    con.execute(f"create temp table deaths as select cast({qi(dp)} as varchar) patid,min(try_cast({qi(dd)} as date)) death_date from d group by 1")
    con.execute("create temp table cohort as select a.*,b.admit_date,b.discharge_date,d.death_date from anchors a join base b using(patid,encounterid) left join deaths d using(patid) where (b.discharge_date is null or b.discharge_date>cast(a.anchor_ts+interval 96 hour as date)) and (d.death_date is null or d.death_date>cast(a.anchor_ts+interval 96 hour as date))")
    n = int(con.execute("select count(*) from cohort").fetchone()[0])

    dictionary: list[dict] = []
    coverage: list[dict] = []
    dx_rows: list[dict] = []
    vaso_rows: list[dict] = []

    def add_dict(variable, domain, status, source, definition, harmonization_note=""):
        dictionary.append({"variable": variable, "domain": domain, "status": status, "source": source, "definition": definition, "harmonization_note": harmonization_note})

    def add_cov(variable, window, count, status="retained"):
        coverage.append({"variable": variable, "window": window, "count": safe(count), "cohort_n": safe(n), "proportion": (count / n if n else None), "status": status})

    # Baseline demographics.
    if paths["dem"] is not None:
        c = cols["dem"]
        pat = first(c, ["PATID"]); birth = first(c, ["BIRTH_DATE"]); sex = first(c, ["SEX"]); race = first(c, ["RACE"])
        if pat and birth:
            cnt = int(con.execute(f"select count(distinct c.encounterid) from cohort c join dem x on cast(x.{qi(pat)} as varchar)=c.patid where try_cast(x.{qi(birth)} as date) is not null").fetchone()[0])
            add_cov("age", "baseline", cnt)
        if pat and sex:
            cnt = int(con.execute(f"select count(distinct c.encounterid) from cohort c join dem x on cast(x.{qi(pat)} as varchar)=c.patid where trim(cast(x.{qi(sex)} as varchar))<>''").fetchone()[0])
            add_cov("sex_male", "baseline", cnt)
        if pat and race:
            cnt = int(con.execute(f"select count(distinct c.encounterid) from cohort c join dem x on cast(x.{qi(pat)} as varchar)=c.patid where trim(cast(x.{qi(race)} as varchar))<>''").fetchone()[0])
            add_cov("race_white", "baseline", cnt)
    add_dict("age", "baseline", "retain", "DEMOGRAPHIC", "Age at hospital admission from BIRTH_DATE and ADMIT_DATE.")
    add_dict("sex_male", "baseline", "retain", "DEMOGRAPHIC", "Binary male indicator from PCORnet SEX.")
    add_dict("race_white", "baseline", "retain_modified", "DEMOGRAPHIC", "Binary White indicator using PCORnet RACE code 05 or literal WHITE.", "Modified coding implementation; same clinical construct as MIMIC.")

    # Diagnosis-based chronic-condition proxies.
    if paths["dx"] is not None:
        c = cols["dx"]
        pat, enc = first(c, ["PATID"]), first(c, ["ENCOUNTERID", "ENCOUNTER_ID"])
        dx = first(c, ["DX", "DIAGNOSIS_CODE"])
        if pat and enc and dx:
            norm = f"upper(regexp_replace(coalesce(cast(x.{qi(dx)} as varchar),''),'[^A-Z0-9]','','g'))"
            concept_counts = {}
            for concept, pattern in DX_PATTERNS.items():
                cnt = int(con.execute(f"select count(distinct c.encounterid) from cohort c join dx x on cast(x.{qi(pat)} as varchar)=c.patid and cast(x.{qi(enc)} as varchar)=c.encounterid where regexp_matches({norm},{q(pattern)})").fetchone()[0])
                concept_counts[concept] = cnt
                dx_rows.append({"construct": concept, "normalized_code_regex": pattern, "encounters": safe(cnt), "cohort_n": safe(n), "proportion": cnt / n if n else None})
            anypat = "|".join(f"(?:{p})" for p in DX_PATTERNS.values())
            anycnt = int(con.execute(f"select count(distinct c.encounterid) from cohort c join dx x on cast(x.{qi(pat)} as varchar)=c.patid and cast(x.{qi(enc)} as varchar)=c.encounterid where regexp_matches({norm},{q(anypat)})").fetchone()[0])
            dx_rows.append({"construct": "comorb", "normalized_code_regex": anypat, "encounters": safe(anycnt), "cohort_n": safe(n), "proportion": anycnt / n if n else None})
            add_cov("comorb", "baseline", anycnt)
            add_cov("heart_failure", "baseline", concept_counts.get("heart_failure", 0))
            add_cov("chronic_kidney", "baseline", concept_counts.get("chronic_kidney", 0))
    add_dict("comorb", "baseline", "retain_modified", "DIAGNOSIS", "Any diabetes, chronic kidney disease, heart failure, COPD, chronic liver disease/cirrhosis, or malignancy by frozen ICD-9/10 prefix sets.", "Code-based chronic proxy replaces MIMIC diagnosis-title matching.")
    add_dict("heart_failure", "baseline", "retain_modified", "DIAGNOSIS", "ICD-10 I50* or ICD-9 428* on the index encounter.")
    add_dict("chronic_kidney", "baseline", "retain_modified", "DIAGNOSIS", "ICD-10 N18* or ICD-9 585* on the index encounter.")

    # Core laboratory trajectories using corrected specimen clock.
    if paths["lab"] is not None:
        c = cols["lab"]
        pat, enc = first(c, ["PATID"]), first(c, ["ENCOUNTERID", "ENCOUNTER_ID"])
        loinc, val = first(c, ["LAB_LOINC", "LOINC"]), first(c, ["RESULT_NUM"])
        sd, st = first(c, ["SPECIMEN_DATE"]), first(c, ["SPECIMEN_TIME"])
        spec_ts = ts("l", sd, st, numeric_seconds=True)
        if all(x is not None for x in [pat, enc, loinc, val, spec_ts]):
            for concept, cfg in LABS.items():
                codes = ",".join(q(x) for x in cfg["codes"])
                valid = f"cast(l.{qi(loinc)} as varchar) in ({codes}) and try_cast(l.{qi(val)} as double) between {cfg['lo']} and {cfg['hi']}"
                for win, a, b in [("0_24h", 0, 24), ("48_72h", 48, 72), ("pre72", 0, 72)]:
                    cnt = int(con.execute(f"select count(distinct c.encounterid) from cohort c join lab l on cast(l.{qi(pat)} as varchar)=c.patid and cast(l.{qi(enc)} as varchar)=c.encounterid where {valid} and {spec_ts}>=c.anchor_ts+interval {a} hour and {spec_ts}<c.anchor_ts+interval {b} hour").fetchone()[0])
                    add_cov(concept, win, cnt)
                add_dict(concept, "laboratory", "retain", "LAB_RESULT_CM", f"Exact LOINC {cfg['codes']}; RESULT_NUM in [{cfg['lo']}, {cfg['hi']}]; specimen timestamp = SPECIMEN_DATE + SPECIMEN_TIME seconds since midnight; early 0-24h and late 48-72h summaries use {cfg['summary']}.")

    # Exact OBS_CLIN physiologic signals.
    if paths["obs"] is not None:
        c = cols["obs"]
        pat, enc = first(c, ["PATID"]), first(c, ["ENCOUNTERID", "ENCOUNTER_ID"])
        codecol = first(c, ["OBSCLIN_CODE"]); val = first(c, ["OBSCLIN_RESULT_NUM"]); unit = first(c, ["OBSCLIN_RESULT_UNIT", "RAW_OBSCLIN_UNIT"])
        od, ot = first(c, ["OBSCLIN_START_DATE"]), first(c, ["OBSCLIN_START_TIME"])
        ots = ts("o", od, ot)
        if all(x is not None for x in [pat, enc, codecol, val, ots]):
            for concept, cfg in OBS.items():
                rawv = f"try_cast(o.{qi(val)} as double)"
                if concept == "temperature" and unit:
                    u = f"upper(trim(coalesce(cast(o.{qi(unit)} as varchar),'')))"
                    vexpr = f"case when {u} in ('F','FAH','FAHRENHEIT') then ({rawv}-32)*5.0/9.0 else {rawv} end"
                else:
                    vexpr = rawv
                valid = f"cast(o.{qi(codecol)} as varchar)={q(cfg['code'])} and ({vexpr}) between {cfg['lo']} and {cfg['hi']}"
                for win, a, b in [("0_24h", 0, 24), ("48_72h", 48, 72), ("pre72", 0, 72)]:
                    cnt = int(con.execute(f"select count(distinct c.encounterid) from cohort c join obs o on cast(o.{qi(pat)} as varchar)=c.patid and cast(o.{qi(enc)} as varchar)=c.encounterid where {valid} and {ots}>=c.anchor_ts+interval {a} hour and {ots}<c.anchor_ts+interval {b} hour").fetchone()[0])
                    add_cov(concept, win, cnt)
                add_dict(concept, "physiology", "retain", "OBS_CLIN", f"Exact OBSCLIN_CODE {cfg['code']}; plausible range [{cfg['lo']}, {cfg['hi']}]; {cfg['summary']} within each window.")

    # MAP derived from paired systolic/diastolic values in VITAL; retain only if measurable.
    map_status = "unavailable"
    if paths["vit"] is not None:
        c = cols["vit"]
        pat, enc = first(c, ["PATID"]), first(c, ["ENCOUNTERID", "ENCOUNTER_ID"])
        sbp, dbp = first(c, ["SYSTOLIC"]), first(c, ["DIASTOLIC"])
        vd, vt = first(c, ["MEASURE_DATE", "VITAL_DATE"]), first(c, ["MEASURE_TIME", "VITAL_TIME"])
        vts = ts("v", vd, vt)
        if all(x is not None for x in [pat, enc, sbp, dbp, vts]):
            mapexpr = f"(try_cast(v.{qi(sbp)} as double)+2*try_cast(v.{qi(dbp)} as double))/3.0"
            valid = f"try_cast(v.{qi(sbp)} as double) between 40 and 300 and try_cast(v.{qi(dbp)} as double) between 20 and 200 and ({mapexpr}) between 20 and 200"
            counts = []
            for win, a, b in [("0_24h", 0, 24), ("48_72h", 48, 72), ("pre72", 0, 72)]:
                cnt = int(con.execute(f"select count(distinct c.encounterid) from cohort c join vit v on cast(v.{qi(pat)} as varchar)=c.patid and cast(v.{qi(enc)} as varchar)=c.encounterid where {valid} and {vts}>=c.anchor_ts+interval {a} hour and {vts}<c.anchor_ts+interval {b} hour").fetchone()[0])
                counts.append(cnt); add_cov("map", win, cnt)
            map_status = "retain_derived" if counts and counts[-1] >= MIN_CELL else "unavailable"
    add_dict("map", "physiology", map_status, "VITAL", "Derived MAP=(SBP+2*DBP)/3 from paired same-row systolic/diastolic values; plausible SBP 40-300, DBP 20-200, MAP 20-200; minimum within window.", "Derived rather than direct MAP observation.")

    # Vasopressor overlap. Missing/invalid stop -> one-hour interval, matching MIMIC repair logic.
    if paths["med"] is not None:
        c = cols["med"]
        pat, enc = first(c, ["PATID"]), first(c, ["ENCOUNTERID", "ENCOUNTER_ID"])
        namecol = first(c, ["RAW_MEDADMIN_MED_NAME", "MEDADMIN_MED_NAME"])
        sd, st = first(c, ["MEDADMIN_START_DATE", "START_DATE"]), first(c, ["MEDADMIN_START_TIME", "START_TIME"])
        ed, et = first(c, ["MEDADMIN_STOP_DATE", "STOP_DATE"]), first(c, ["MEDADMIN_STOP_TIME", "STOP_TIME"])
        start = ts("m", sd, st); stopraw = ts("m", ed, et) if ed else None
        if all(x is not None for x in [pat, enc, namecol, start]):
            nm = f"lower(coalesce(cast(m.{qi(namecol)} as varchar),''))"
            inc = " or ".join(f"strpos({nm},{q(x)})>0" for x in VASO_TERMS)
            excp = " or ".join(f"strpos({nm},{q(x)})>0" for x in VASO_EXCLUDE)
            stop = f"case when {stopraw} is null or {stopraw}<={start} then {start}+interval 1 hour else {stopraw} end" if stopraw else f"{start}+interval 1 hour"
            basejoin = f"cast(m.{qi(pat)} as varchar)=c.patid and cast(m.{qi(enc)} as varchar)=c.encounterid and ({inc}) and not ({excp}) and {start} is not null"
            vals = {}
            for var, a, b in [("vasopressor_any_0_24h", 0, 24), ("vasopressor_any_48_72h", 48, 72), ("vasopressor_any_pre72", 0, 72)]:
                cnt = int(con.execute(f"select count(distinct c.encounterid) from cohort c join med m on {basejoin} where {start}<c.anchor_ts+interval {b} hour and {stop}>c.anchor_ts+interval {a} hour").fetchone()[0])
                vals[var] = cnt; add_cov(var, f"{a}_{b}h", cnt)
            stopped = int(con.execute(f"select count(*) from (select c.encounterid,max(case when {start}<c.anchor_ts+interval 24 hour and {stop}>c.anchor_ts then 1 else 0 end) early,max(case when {start}<c.anchor_ts+interval 72 hour and {stop}>c.anchor_ts+interval 48 hour then 1 else 0 end) late from cohort c left join med m on {basejoin} group by 1) z where early=1 and late=0").fetchone()[0])
            add_cov("vasopressor_stopped_before_72h", "trajectory", stopped)
            vaso_rows.extend([
                {"construct": "vasopressor_any_0_24h", "encounters": safe(vals.get("vasopressor_any_0_24h")), "cohort_n": safe(n), "proportion": vals.get("vasopressor_any_0_24h", 0) / n if n else None},
                {"construct": "vasopressor_any_48_72h", "encounters": safe(vals.get("vasopressor_any_48_72h")), "cohort_n": safe(n), "proportion": vals.get("vasopressor_any_48_72h", 0) / n if n else None},
                {"construct": "vasopressor_stopped_before_72h", "encounters": safe(stopped), "cohort_n": safe(n), "proportion": stopped / n if n else None},
            ])
    add_dict("vasopressor_any_0_24h", "severity", "retain", "MED_ADMIN", "Validated norepinephrine/levophed, phenylephrine, vasopressin, epinephrine, or dopamine interval overlaps [anchor,anchor+24h); exclude racepinephrine/nasal/ophthalmic/lidocaine/topical false positives.")
    add_dict("vasopressor_any_48_72h", "severity", "retain", "MED_ADMIN", "Same validated vasopressor interval-overlap rule for [anchor+48h,anchor+72h).")
    add_dict("vasopressor_stopped_before_72h", "trajectory", "retain", "MED_ADMIN", "Early vasopressor overlap=1 and late overlap=0.")

    # Prespecified harmonization exclusions/modified terms.
    for variable, domain, status, source, definition, note in [
        ("gcs_total_48_72h", "physiology", "exclude_primary", "OBS_CLIN", "No validated PSU GCS mapping in current extract.", "Also excluded as a direct term from corrected MIMIC primary PS."),
        ("fio2_48_72h", "physiology", "exclude_primary", "OBS_CLIN", "Sparse PSU FiO2 signal; not retained.", "Also excluded as a direct term from corrected MIMIC primary PS."),
        ("vent_proc", "severity", "exclude_primary", "PROCEDURES", "No validated procedures source in current PSU extract.", "Unavailable rather than opportunistically substituted."),
        ("urine_output_ml_48_72h", "severity", "exclude_primary", "none", "No harmonized urine-output source frozen for PSU.", "Unavailable rather than opportunistically substituted."),
        ("hours_admit_to_icu", "baseline", "exclude_primary", "none", "Exact ICU admission clock unavailable; PSU anchor is hospital-clock modified replication.", "Structural site difference."),
        ("icu_type_indicators", "baseline", "exclude_primary", "none", "No validated ICU care-unit field.", "Structural site difference."),
        ("microbiology_intensity_terms", "diagnostic_intensity", "exclude_primary", "LAB_RESULT_CM", "Parent clinical-culture result positivity/availability cannot be faithfully linked in current PSU extract.", "Do not substitute keyword-derived microbiology terms."),
        ("broad_abx_hours_pre72", "antibiotic_intensity", "exclude_primary", "PRESCRIBING", "PRESCRIBING duration is date-level; exact exposure hours are not harmonizable.", "Retain agent/class intensity terms later, not exact hours."),
        ("bmi_pre72", "baseline", "exclude_primary", "VITAL", "Height/weight exist but a harmonized pre-72h BMI construction was not validated in the freeze sequence.", "Avoid adding a site-specific term after seeing balance."),
    ]:
        add_dict(variable, domain, status, source, definition, note)

    pd.DataFrame(dictionary).to_csv(args.output_dir / "final_covariate_dictionary.csv", index=False)
    pd.DataFrame(coverage).to_csv(args.output_dir / "final_covariate_coverage.csv", index=False)
    pd.DataFrame(dx_rows).to_csv(args.output_dir / "diagnosis_definition_summary.csv", index=False)
    pd.DataFrame(vaso_rows).to_csv(args.output_dir / "vasopressor_window_overlap.csv", index=False)

    retained = [r["variable"] for r in dictionary if r["status"].startswith("retain")]
    excluded = [r["variable"] for r in dictionary if r["status"].startswith("exclude") or r["status"] == "unavailable"]
    summary = {
        "privacy_mode": "aggregate_only",
        "strict_cohort_n": safe(n),
        "replication_type": "modified external replication using hospital clock",
        "lab_clock": "SPECIMEN_DATE plus SPECIMEN_TIME decoded as seconds since midnight",
        "lab_result_availability_warning": "SPECIMEN and RESULT clocks are identical in the audited extract; do not interpret RESULT time as independent result-availability time",
        "retained_covariate_constructs": retained,
        "excluded_or_unavailable_constructs": excluded,
        "map_status": map_status,
        "next_step": "build PSU analysis dataset and run propensity-score balance diagnostics without treatment-effect estimation",
        "guardrail": "No propensity score or treatment effects fit in this task.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

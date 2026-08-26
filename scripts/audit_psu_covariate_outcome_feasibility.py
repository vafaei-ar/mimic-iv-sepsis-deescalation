#!/usr/bin/env python3
"""Aggregate-only PSU covariate and outcome feasibility audit.

Builds the frozen strict modified PSU 96-hour landmark cohort, then quantifies whether
key MIMIC covariate domains and candidate outcomes are reproducible from available
PCORnet/local sources. This is a feasibility audit only: no patient rows, identifiers,
free-text values, propensity scores, or treatment effects are exported.
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
LAB_PATTERNS = {
    "lactate": r"lactate|lactic acid",
    "creatinine": r"creatinine",
    "wbc": r"white blood|\\bwbc\\b",
    "platelet": r"platelet",
    "bilirubin": r"bilirubin",
}
OBS_PATTERNS = {
    "heart_rate": r"heart rate|pulse rate",
    "resp_rate": r"respiratory rate|resp rate",
    "spo2": r"spo2|oxygen saturation|o2 saturation",
    "temperature": r"temperature|temp",
    "gcs": r"glasgow|\\bgcs\\b",
    "fio2": r"fio2|fraction.*inspired.*oxygen",
}
VASO_PATTERN = r"norepinephrine|levophed|phenylephrine|vasopressin|epinephrine|dopamine"
VENT_PATTERN = r"mechanical ventilation|respiratory ventilation|intubation|endotracheal"


def q(v: str) -> str:
    return "'" + v.replace("'", "''") + "'"


def qi(v: str) -> str:
    return '"' + v.replace('"', '""') + '"'


def safe(n: int | float | None):
    if n is None:
        return None
    n = int(n)
    return n if n >= MIN_CELL else None


def find_parquet(root: Path, stems: list[str], required: bool = True) -> Path | None:
    for stem in stems:
        c = sorted(root.glob(f"PCORnet/parquet/**/{stem}*.parquet"))
        if not c:
            c = sorted(root.glob(f"PCORnet/parquet/**/{stem.replace('_', '')}*.parquet"))
        if c:
            exact = [p for p in c if p.stem.lower() == stem.lower()]
            return exact[0] if exact else max(c, key=lambda p: p.stat().st_size)
    if required:
        raise FileNotFoundError(f"No parquet found for {stems}")
    return None


def first_present(cols: set[str], names: list[str]) -> str | None:
    lut = {c.upper(): c for c in cols}
    for n in names:
        if n.upper() in lut:
            return lut[n.upper()]
    return None


def parse_legacy_codes(path: Path) -> tuple[set[str], set[str]]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    tree = ast.parse(text)
    include: set[str] = set()
    exclude: set[str] = set()
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        targets = [t.id for t in node.targets if isinstance(t, ast.Name)]
        if not targets:
            continue
        try:
            value = ast.literal_eval(node.value)
        except Exception:
            continue
        for target in targets:
            low = target.lower()
            if "exclude" in low and isinstance(value, (set, list, tuple)):
                for item in value:
                    tok = str(item).strip()
                    if re.fullmatch(r"\d{4,9}", tok):
                        exclude.add(tok)
            elif isinstance(value, (list, tuple)):
                for item in value:
                    if isinstance(item, (list, tuple)) and item:
                        tok = str(item[0]).strip()
                        if re.fullmatch(r"\d{4,9}", tok):
                            include.add(tok)
    return include - exclude, exclude


def describe(con: duckdb.DuckDBPyConnection, view: str) -> set[str]:
    return set(con.execute(f"DESCRIBE {view}").fetchdf()["column_name"].astype(str))


def timestamp_expr(alias: str, date_col: str | None, time_col: str | None) -> str | None:
    if not date_col:
        return None
    if time_col:
        return (
            f"try_cast(cast({alias}.{qi(date_col)} AS VARCHAR) || ' ' || "
            f"coalesce(cast({alias}.{qi(time_col)} AS VARCHAR),'00:00:00') AS TIMESTAMP)"
        )
    return f"try_cast({alias}.{qi(date_col)} AS TIMESTAMP)"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    root = args.data_root

    paths = {
        "sepsis": find_parquet(root, ["sepsis_encounter"]),
        "prescribing": find_parquet(root, ["prescribing"]),
        "death": find_parquet(root, ["death"]),
        "demographic": find_parquet(root, ["sepsis_demographic", "demographic"], required=False),
        "diagnosis": find_parquet(root, ["sepsis_diagnosis", "diagnosis"], required=False),
        "vital": find_parquet(root, ["sepsis_vital", "vital"], required=False),
        "lab": find_parquet(root, ["lab_reduced", "lab_result_cm"], required=False),
        "med_admin": find_parquet(root, ["med_admin"], required=False),
        "obs_clin": find_parquet(root, ["obs_clin"], required=False),
        "procedures": find_parquet(root, ["procedures", "procedure"], required=False),
    }
    legacy = root / "PCORnet" / "code" / "config" / "codes_antibiotics.py"
    if not legacy.exists():
        raise FileNotFoundError(legacy)
    include_codes, exclude_codes = parse_legacy_codes(legacy)
    include_sql = ",".join(q(x) for x in sorted(include_codes)) or "''"
    exclude_sql = ",".join(q(x) for x in sorted(exclude_codes)) or "''"

    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    for name, path in paths.items():
        if path is not None:
            con.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet({q(str(path))})")

    s_cols, p_cols, d_cols = describe(con, "sepsis"), describe(con, "prescribing"), describe(con, "death")
    s_patid = first_present(s_cols, ["PATID"]); s_enc = first_present(s_cols, ["ENCOUNTERID", "ENCOUNTER_ID"])
    admit = first_present(s_cols, ["ADMIT_DATE", "ADMITDATE"]); discharge = first_present(s_cols, ["DISCHARGE_DATE", "DISCHARGEDATE"])
    p_patid = first_present(p_cols, ["PATID"]); p_enc = first_present(p_cols, ["ENCOUNTERID", "ENCOUNTER_ID"])
    p_order_date = first_present(p_cols, ["RX_ORDER_DATE", "ORDER_DATE"]); p_order_time = first_present(p_cols, ["RX_ORDER_TIME", "ORDER_TIME"])
    p_name = first_present(p_cols, ["RAW_RX_MED_NAME", "RX_MED_NAME", "MEDICATION_NAME"]); p_code = first_present(p_cols, ["RXNORM_CUI", "RAW_RXNORM_CUI", "RXNORM"])
    p_route = first_present(p_cols, ["RX_ROUTE", "RAW_RX_ROUTE", "ROUTE"])
    d_patid = first_present(d_cols, ["PATID"]); death_date = first_present(d_cols, ["DEATH_DATE", "DEATHDATE"])
    required = [s_patid,s_enc,admit,discharge,p_patid,p_enc,p_order_date,p_name,p_code,d_patid,death_date]
    if any(x is None for x in required):
        raise RuntimeError("Required core fields missing")

    p_order_ts = timestamp_expr("p", p_order_date, p_order_time)
    name_expr = f"lower(coalesce(cast(p.{qi(p_name)} AS VARCHAR),''))"
    code_expr = f"trim(coalesce(cast(p.{qi(p_code)} AS VARCHAR),''))"
    route_expr = f"upper(trim(coalesce(cast(p.{qi(p_route)} AS VARCHAR),'')))" if p_route else "''"
    broad = (
        f"(({code_expr} IN ({include_sql}) OR regexp_matches({name_expr},{q(BROAD_PATTERN)})) "
        f"AND {code_expr} NOT IN ({exclude_sql}) AND NOT regexp_matches({name_expr},{q(NON_SYSTEMIC_PATTERN)}) "
        f"AND {route_expr} NOT IN ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    )

    con.execute(f"""
      CREATE TEMP TABLE base AS
      SELECT DISTINCT cast({qi(s_patid)} AS VARCHAR) patid,
             cast({qi(s_enc)} AS VARCHAR) encounterid,
             try_cast({qi(admit)} AS DATE) admit_date,
             try_cast({qi(discharge)} AS DATE) discharge_date
      FROM sepsis WHERE {qi(s_patid)} IS NOT NULL AND {qi(s_enc)} IS NOT NULL
    """)
    con.execute(f"""
      CREATE TEMP TABLE anchor_orders AS
      SELECT cast(p.{qi(p_patid)} AS VARCHAR) patid,
             cast(p.{qi(p_enc)} AS VARCHAR) encounterid,
             {p_order_ts} order_ts
      FROM prescribing p JOIN base b
        ON cast(p.{qi(p_patid)} AS VARCHAR)=b.patid
       AND cast(p.{qi(p_enc)} AS VARCHAR)=b.encounterid
      WHERE {broad} AND {p_order_ts} IS NOT NULL AND b.admit_date IS NOT NULL
        AND {p_order_ts} >= cast(b.admit_date AS TIMESTAMP)
        AND {p_order_ts} < cast(b.admit_date AS TIMESTAMP)+INTERVAL 24 HOUR
    """)
    con.execute("CREATE TEMP TABLE anchors AS SELECT patid,encounterid,min(order_ts) anchor_ts FROM anchor_orders GROUP BY 1,2")
    con.execute(f"CREATE TEMP TABLE deaths AS SELECT cast({qi(d_patid)} AS VARCHAR) patid,min(try_cast({qi(death_date)} AS DATE)) death_date FROM death GROUP BY 1")
    con.execute("""
      CREATE TEMP TABLE cohort AS
      SELECT a.patid,a.encounterid,a.anchor_ts,b.admit_date,b.discharge_date,d.death_date
      FROM anchors a JOIN base b USING(patid,encounterid) LEFT JOIN deaths d USING(patid)
      WHERE (b.discharge_date IS NULL OR b.discharge_date > cast(a.anchor_ts+INTERVAL 96 HOUR AS DATE))
        AND (d.death_date IS NULL OR d.death_date > cast(a.anchor_ts+INTERVAL 96 HOUR AS DATE))
    """)
    cohort_n = int(con.execute("SELECT count(*) FROM cohort").fetchone()[0])

    source_rows = []
    for source in ["demographic","diagnosis","vital","lab","med_admin","obs_clin","procedures","death"]:
        if source not in paths or paths[source] is None:
            source_rows.append({"source": source, "present": False, "cohort_encounters_with_any_row": None})
            continue
        cols = describe(con, source)
        pat = first_present(cols, ["PATID"])
        enc = first_present(cols, ["ENCOUNTERID", "ENCOUNTER_ID"])
        if not pat:
            source_rows.append({"source": source, "present": True, "cohort_encounters_with_any_row": None})
            continue
        if enc:
            sql = f"SELECT count(DISTINCT c.encounterid) FROM cohort c JOIN {source} x ON cast(x.{qi(pat)} AS VARCHAR)=c.patid AND cast(x.{qi(enc)} AS VARCHAR)=c.encounterid"
        else:
            sql = f"SELECT count(DISTINCT c.encounterid) FROM cohort c JOIN {source} x ON cast(x.{qi(pat)} AS VARCHAR)=c.patid"
        n = int(con.execute(sql).fetchone()[0])
        source_rows.append({"source": source, "present": True, "cohort_encounters_with_any_row": safe(n)})
    pd.DataFrame(source_rows).to_csv(args.output_dir / "source_cohort_coverage.csv", index=False)

    cov_rows: list[dict] = []
    def cov(domain: str, construct: str, status: str, source: str, count=None, note=""):
        cov_rows.append({"domain": domain, "construct": construct, "status": status, "source": source, "cohort_encounters_with_signal": safe(count) if count is not None else None, "cohort_n": safe(cohort_n), "note": note})

    # Demographics and diagnoses: availability gate only; coding details are deferred.
    if paths["demographic"] is not None:
        dem_cols = describe(con, "demographic")
        cov("baseline","demographics","available","DEMOGRAPHIC",None,"Patient-level demographic source present; exact age/sex/race derivation requires field-level mapping in the build step.")
    else:
        cov("baseline","demographics","not_available","DEMOGRAPHIC")
    if paths["diagnosis"] is not None:
        cov("baseline","comorbidity_diagnoses","available","DIAGNOSIS",None,"Encounter diagnosis source present; chronic-condition code sets still need harmonized definitions.")
    else:
        cov("baseline","comorbidity_diagnoses","not_available","DIAGNOSIS")

    # Lab trajectory coverage by name/LOINC metadata; no raw names exported.
    if paths["lab"] is not None:
        lc = describe(con, "lab")
        l_pat = first_present(lc,["PATID"]); l_enc = first_present(lc,["ENCOUNTERID","ENCOUNTER_ID"])
        l_name = first_present(lc,["RAW_LAB_NAME","LAB_NAME"]); l_num = first_present(lc,["RESULT_NUM","LAB_RESULT_NUM","RESULT_NUMERIC"])
        l_date = first_present(lc,["SPECIMEN_DATE","RESULT_DATE","LAB_ORDER_DATE"]); l_time = first_present(lc,["SPECIMEN_TIME","RESULT_TIME"])
        l_ts = timestamp_expr("l", l_date, l_time)
        if l_pat and l_enc and l_name and l_num and l_ts:
            name_l = f"lower(coalesce(cast(l.{qi(l_name)} AS VARCHAR),''))"
            for concept, pattern in LAB_PATTERNS.items():
                for win, lo, hi in [("0_24h",0,24),("48_72h",48,72),("pre72",0,72)]:
                    n = int(con.execute(f"""
                      SELECT count(DISTINCT c.encounterid)
                      FROM cohort c JOIN lab l
                        ON cast(l.{qi(l_pat)} AS VARCHAR)=c.patid AND cast(l.{qi(l_enc)} AS VARCHAR)=c.encounterid
                      WHERE regexp_matches({name_l},{q(pattern)})
                        AND try_cast(l.{qi(l_num)} AS DOUBLE) IS NOT NULL
                        AND {l_ts} >= c.anchor_ts + INTERVAL {lo} HOUR
                        AND {l_ts} <= c.anchor_ts + INTERVAL {hi} HOUR
                    """).fetchone()[0])
                    cov("trajectory",f"{concept}_{win}","measurable","LAB_RESULT_CM",n,"Concept identified from local lab-name metadata; mapping should be frozen before final analysis.")
        else:
            for concept in LAB_PATTERNS:
                cov("trajectory",concept,"field_mapping_incomplete","LAB_RESULT_CM")

    # OBS_CLIN candidate ICU-style vital/neurologic variables.
    if paths["obs_clin"] is not None:
        oc = describe(con,"obs_clin")
        o_pat = first_present(oc,["PATID"]); o_enc = first_present(oc,["ENCOUNTERID","ENCOUNTER_ID"])
        o_name = first_present(oc,["RAW_OBSCLIN_NAME","OBSCLIN_NAME","RAW_OBS_NAME","OBS_NAME"])
        o_date = first_present(oc,["OBSCLIN_START_DATE","OBS_DATE","MEASURE_DATE"]); o_time = first_present(oc,["OBSCLIN_START_TIME","OBS_TIME","MEASURE_TIME"])
        o_ts = timestamp_expr("o",o_date,o_time)
        if o_pat and o_enc and o_name and o_ts:
            o_name_expr = f"lower(coalesce(cast(o.{qi(o_name)} AS VARCHAR),''))"
            for concept, pattern in OBS_PATTERNS.items():
                for win, lo, hi in [("0_24h",0,24),("48_72h",48,72)]:
                    n = int(con.execute(f"""
                      SELECT count(DISTINCT c.encounterid)
                      FROM cohort c JOIN obs_clin o
                        ON cast(o.{qi(o_pat)} AS VARCHAR)=c.patid AND cast(o.{qi(o_enc)} AS VARCHAR)=c.encounterid
                      WHERE regexp_matches({o_name_expr},{q(pattern)})
                        AND {o_ts} >= c.anchor_ts + INTERVAL {lo} HOUR
                        AND {o_ts} <= c.anchor_ts + INTERVAL {hi} HOUR
                    """).fetchone()[0])
                    cov("trajectory",f"{concept}_{win}","candidate_signal","OBS_CLIN",n,"Name-based feasibility signal only; code/value semantics require a separate frozen map before modeling.")
        else:
            for concept in OBS_PATTERNS:
                cov("trajectory",concept,"field_mapping_incomplete","OBS_CLIN")
    else:
        for concept in OBS_PATTERNS:
            cov("trajectory",concept,"not_available","OBS_CLIN")

    # VITAL is expected to cover BP/anthropometrics rather than all MIMIC ICU vitals.
    if paths["vital"] is not None:
        vc = describe(con,"vital")
        v_pat = first_present(vc,["PATID"]); v_enc = first_present(vc,["ENCOUNTERID","ENCOUNTER_ID"])
        v_date = first_present(vc,["MEASURE_DATE","VITAL_DATE"]); v_time = first_present(vc,["MEASURE_TIME","VITAL_TIME"])
        v_ts = timestamp_expr("v",v_date,v_time)
        if v_pat and v_enc and v_ts:
            n0 = int(con.execute(f"SELECT count(DISTINCT c.encounterid) FROM cohort c JOIN vital v ON cast(v.{qi(v_pat)} AS VARCHAR)=c.patid AND cast(v.{qi(v_enc)} AS VARCHAR)=c.encounterid WHERE {v_ts}>=c.anchor_ts AND {v_ts}<=c.anchor_ts+INTERVAL 72 HOUR").fetchone()[0])
            cov("trajectory","vital_table_pre72_any","available","VITAL",n0,"VITAL has timestamped measurements, but concept coverage is limited compared with MIMIC chartevents.")

    # Vasopressor trajectory from MED_ADMIN.
    if paths["med_admin"] is not None:
        mc = describe(con,"med_admin")
        m_pat = first_present(mc,["PATID"]); m_enc = first_present(mc,["ENCOUNTERID","ENCOUNTER_ID"])
        m_name = first_present(mc,["RAW_MEDADMIN_MED_NAME","MEDADMIN_MED_NAME","MEDICATION_NAME"])
        m_date = first_present(mc,["MEDADMIN_START_DATE","START_DATE"]); m_time = first_present(mc,["MEDADMIN_START_TIME","START_TIME"])
        m_ts = timestamp_expr("m",m_date,m_time)
        if m_pat and m_enc and m_name and m_ts:
            mname = f"lower(coalesce(cast(m.{qi(m_name)} AS VARCHAR),''))"
            for win,lo,hi in [("0_24h",0,24),("48_72h",48,72)]:
                n = int(con.execute(f"SELECT count(DISTINCT c.encounterid) FROM cohort c JOIN med_admin m ON cast(m.{qi(m_pat)} AS VARCHAR)=c.patid AND cast(m.{qi(m_enc)} AS VARCHAR)=c.encounterid WHERE regexp_matches({mname},{q(VASO_PATTERN)}) AND {m_ts}>=c.anchor_ts+INTERVAL {lo} HOUR AND {m_ts}<=c.anchor_ts+INTERVAL {hi} HOUR").fetchone()[0])
                cov("trajectory",f"vasopressor_any_{win}","measurable","MED_ADMIN",n,"Administration timestamps available; interval logic can be harmonized in the final build.")

    # Mechanical ventilation/source-control feasibility from procedures.
    if paths["procedures"] is not None:
        pc = describe(con,"procedures")
        pr_pat = first_present(pc,["PATID"]); pr_enc = first_present(pc,["ENCOUNTERID","ENCOUNTER_ID"])
        pr_name = first_present(pc,["RAW_PX_NAME","PX_NAME","PROCEDURE_NAME"])
        if pr_pat and pr_enc and pr_name:
            pname = f"lower(coalesce(cast(r.{qi(pr_name)} AS VARCHAR),''))"
            n = int(con.execute(f"SELECT count(DISTINCT c.encounterid) FROM cohort c JOIN procedures r ON cast(r.{qi(pr_pat)} AS VARCHAR)=c.patid AND cast(r.{qi(pr_enc)} AS VARCHAR)=c.encounterid WHERE regexp_matches({pname},{q(VENT_PATTERN)})").fetchone()[0])
            cov("baseline","mechanical_ventilation_procedure","candidate_signal","PROCEDURES",n,"Name-based feasibility only; exact code set/timing must be frozen before analysis.")
        else:
            cov("baseline","mechanical_ventilation_procedure","field_mapping_incomplete","PROCEDURES")

    pd.DataFrame(cov_rows).to_csv(args.output_dir / "covariate_trajectory_coverage.csv", index=False)

    out_rows = []
    def out(name,status,source,note,count=None):
        out_rows.append({"outcome":name,"status":status,"source":source,"cohort_encounters_with_observable_signal":safe(count) if count is not None else None,"cohort_n":safe(cohort_n),"note":note})
    death_any = int(con.execute("SELECT count(*) FROM cohort WHERE death_date IS NOT NULL").fetchone()[0])
    death30 = int(con.execute("SELECT count(*) FROM cohort WHERE death_date IS NOT NULL AND death_date <= cast(anchor_ts+INTERVAL 30 DAY AS DATE)").fetchone()[0])
    discharge_known = int(con.execute("SELECT count(*) FROM cohort WHERE discharge_date IS NOT NULL").fetchone()[0])
    out("30_day_mortality","feasible_date_level","DEATH", "Death is date-level, not exact time; suitable for 30-day mortality with explicit date-level semantics.", death30)
    out("all_observed_death","feasible_date_level","DEATH","Out-of-system completeness still depends on local death-source capture.",death_any)
    out("hospital_length_of_stay","approximate_date_level","sepsis_encounter","Admission and discharge are calendar dates; exact hours are unavailable.",discharge_known)
    out("hospital_free_days_30","approximate_date_level","sepsis_encounter + DEATH","Can be approximated with date-level discharge and death-to-zero, but not exactly harmonized to MIMIC hour-level timing.",discharge_known)
    out("antibiotic_free_days","candidate","PRESCRIBING / MED_ADMIN","Medication data are available, but interval semantics differ: PRESCRIBING is partly date-level and MED_ADMIN is timestamped.",cohort_n)
    out("readmission_30d","not_validated","current sepsis encounter subset","Current sepsis_encounter extract should not be assumed to contain all subsequent hospitalizations; do not use until a complete encounter source is verified.")
    out("late_recurrent_persistent_antibiotic_use","exploratory_candidate","PRESCRIBING / MED_ADMIN","Technically derivable, but remains discharge/observation-time sensitive and should stay exploratory.",cohort_n)
    pd.DataFrame(out_rows).to_csv(args.output_dir / "outcome_feasibility.csv", index=False)

    fields = []
    for name,path in paths.items():
        if path is None:
            fields.append({"source":name,"present":False,"columns":None})
        else:
            cols = sorted(describe(con,name))
            fields.append({"source":name,"present":True,"columns":";".join(cols)})
    pd.DataFrame(fields).to_csv(args.output_dir / "source_field_inventory.csv", index=False)

    summary = {
        "privacy_mode":"aggregate_only_no_ids_no_patient_rows_no_free_text_export",
        "minimum_reported_cell":MIN_CELL,
        "strict_modified_psu_cohort":safe(cohort_n),
        "clock":"hospital admission calendar date at midnight; exact ICU entry unavailable",
        "exposure":"frozen systemic broad-spectrum PRESCRIBING proxy; strict 96h landmark",
        "purpose":"covariate and outcome feasibility before any propensity-score fitting or treatment-effect estimation",
        "microbiology":"day-3 culture positivity remains not faithfully recoverable and is not used here",
        "guardrail":"Do not fit PS/effect models until covariate mappings and outcome definitions are frozen from this audit."
    }
    (args.output_dir/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")


if __name__ == "__main__":
    main()

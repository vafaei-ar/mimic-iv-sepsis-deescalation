#!/usr/bin/env python3
"""Aggregate-only feasibility audit of PSU 72-96h antibiotic exposure classification.

Uses the frozen modified PSU hospital-clock cohort and systemic broad-spectrum medication
proxy. Primary classification is PRESCRIBING-based. MED_ADMIN is summarized as a
sensitivity source. Exports aggregate counts only, never identifiers or patient rows.
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
NON_BROAD_PATTERN = (
    "ceftriaxone|cefazolin|ampicillin|amoxicillin|doxycycline|azithromycin|"
    "metronidazole|clindamycin|cephalexin|ciprofloxacin|levofloxacin|gentamicin|tobramycin"
)
NON_SYSTEMIC_PATTERN = "cayston|inhal|nebul|tablet|capsule|oral solution|oral suspension|by mouth|\\bpo\\b"


def q(v: str) -> str:
    return "'" + v.replace("'", "''") + "'"


def qident(v: str) -> str:
    return '"' + v.replace('"', '""') + '"'


def find_parquet(root: Path, stem: str) -> Path:
    c = sorted(root.glob(f"PCORnet/parquet/**/{stem}*.parquet"))
    if not c:
        c = sorted(root.glob(f"PCORnet/parquet/**/{stem.replace('_', '')}*.parquet"))
    if not c:
        raise FileNotFoundError(f"No parquet found for {stem}")
    exact = [p for p in c if p.stem.lower() == stem.lower()]
    return exact[0] if exact else max(c, key=lambda p: p.stat().st_size)


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


def safe(n: int) -> int | None:
    return int(n) if int(n) >= MIN_CELL else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    root = args.data_root
    sepsis = find_parquet(root, "sepsis_encounter")
    prescribing = find_parquet(root, "prescribing")
    med_admin = find_parquet(root, "med_admin")
    death = find_parquet(root, "death")
    legacy = root / "PCORnet" / "code" / "config" / "codes_antibiotics.py"
    include_codes, exclude_codes = parse_legacy_codes(legacy)
    include_sql = ",".join(q(x) for x in sorted(include_codes)) or "''"
    exclude_sql = ",".join(q(x) for x in sorted(exclude_codes)) or "''"

    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    for name, path in [("s", sepsis), ("p", prescribing), ("m", med_admin), ("d", death)]:
        con.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet({q(str(path))})")

    cols = {name: set(con.execute(f"DESCRIBE {name}").fetchdf()["column_name"].astype(str)) for name in ["s", "p", "m", "d"]}
    s_patid = first_present(cols["s"], ["PATID"]); s_enc = first_present(cols["s"], ["ENCOUNTERID", "ENCOUNTER_ID"])
    admit = first_present(cols["s"], ["ADMIT_DATE", "ADMITDATE"]); discharge = first_present(cols["s"], ["DISCHARGE_DATE", "DISCHARGEDATE"])
    p_patid = first_present(cols["p"], ["PATID"]); p_enc = first_present(cols["p"], ["ENCOUNTERID", "ENCOUNTER_ID"])
    p_order_date = first_present(cols["p"], ["RX_ORDER_DATE", "ORDER_DATE"]); p_order_time = first_present(cols["p"], ["RX_ORDER_TIME", "ORDER_TIME"])
    p_start = first_present(cols["p"], ["RX_START_DATE", "START_DATE"]); p_end = first_present(cols["p"], ["RX_END_DATE", "END_DATE"])
    p_name = first_present(cols["p"], ["RAW_RX_MED_NAME", "RX_MED_NAME", "MEDICATION_NAME"]); p_code = first_present(cols["p"], ["RXNORM_CUI", "RAW_RXNORM_CUI", "RXNORM"])
    p_route = first_present(cols["p"], ["RX_ROUTE", "RAW_RX_ROUTE", "ROUTE"])
    m_patid = first_present(cols["m"], ["PATID"]); m_enc = first_present(cols["m"], ["ENCOUNTERID", "ENCOUNTER_ID"])
    m_start_date = first_present(cols["m"], ["MEDADMIN_START_DATE", "START_DATE"]); m_start_time = first_present(cols["m"], ["MEDADMIN_START_TIME", "START_TIME"])
    m_stop_date = first_present(cols["m"], ["MEDADMIN_STOP_DATE", "STOP_DATE"]); m_stop_time = first_present(cols["m"], ["MEDADMIN_STOP_TIME", "STOP_TIME"])
    m_name = first_present(cols["m"], ["RAW_MEDADMIN_MED_NAME", "MEDADMIN_MED_NAME", "MEDICATION_NAME"]); m_code = first_present(cols["m"], ["MEDADMIN_CODE", "RXNORM_CUI"])
    m_route = first_present(cols["m"], ["MEDADMIN_ROUTE", "ROUTE"])
    d_patid = first_present(cols["d"], ["PATID"]); death_date = first_present(cols["d"], ["DEATH_DATE", "DEATHDATE"])

    required = [s_patid,s_enc,admit,discharge,p_patid,p_enc,p_order_date,p_name,p_code,m_patid,m_enc,m_start_date,m_name,m_code,d_patid,death_date]
    if any(x is None for x in required):
        raise RuntimeError("Required source columns are missing")

    p_order_ts = f"try_cast(cast(p.{qident(p_order_date)} AS VARCHAR) || ' ' || coalesce(cast(p.{qident(p_order_time)} AS VARCHAR),'00:00:00') AS TIMESTAMP)" if p_order_time else f"try_cast(p.{qident(p_order_date)} AS TIMESTAMP)"
    p_name_expr = f"lower(coalesce(cast(p.{qident(p_name)} AS VARCHAR),''))"; p_code_expr = f"trim(coalesce(cast(p.{qident(p_code)} AS VARCHAR),''))"
    p_route_expr = f"upper(trim(coalesce(cast(p.{qident(p_route)} AS VARCHAR),'')))" if p_route else "''"
    p_broad = f"(({p_code_expr} IN ({include_sql}) OR regexp_matches({p_name_expr},{q(BROAD_PATTERN)})) AND {p_code_expr} NOT IN ({exclude_sql}) AND NOT regexp_matches({p_name_expr},{q(NON_SYSTEMIC_PATTERN)}) AND {p_route_expr} NOT IN ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    p_nonbroad = f"(regexp_matches({p_name_expr},{q(NON_BROAD_PATTERN)}) AND NOT regexp_matches({p_name_expr},{q(NON_SYSTEMIC_PATTERN)}) AND {p_route_expr} NOT IN ('ORAL','RESPIRATORY_TRACT','INHALATION'))"

    con.execute(f"CREATE TEMP TABLE base AS SELECT DISTINCT cast({qident(s_patid)} AS VARCHAR) patid, cast({qident(s_enc)} AS VARCHAR) encounterid, try_cast({qident(admit)} AS DATE) admit_date, try_cast({qident(discharge)} AS DATE) discharge_date FROM s WHERE {qident(s_patid)} IS NOT NULL AND {qident(s_enc)} IS NOT NULL")
    con.execute(f"CREATE TEMP TABLE anchor_orders AS SELECT cast(p.{qident(p_patid)} AS VARCHAR) patid, cast(p.{qident(p_enc)} AS VARCHAR) encounterid, {p_order_ts} order_ts FROM p JOIN base b ON cast(p.{qident(p_patid)} AS VARCHAR)=b.patid AND cast(p.{qident(p_enc)} AS VARCHAR)=b.encounterid WHERE {p_broad} AND {p_order_ts} >= cast(b.admit_date AS TIMESTAMP) AND {p_order_ts} < cast(b.admit_date AS TIMESTAMP)+INTERVAL 24 HOUR")
    con.execute("CREATE TEMP TABLE anchors AS SELECT patid,encounterid,min(order_ts) anchor_ts FROM anchor_orders GROUP BY 1,2")
    con.execute(f"CREATE TEMP TABLE deaths AS SELECT cast({qident(d_patid)} AS VARCHAR) patid,min(try_cast({qident(death_date)} AS DATE)) death_date FROM d GROUP BY 1")
    con.execute("CREATE TEMP TABLE cohort AS SELECT a.*,b.discharge_date,d.death_date FROM anchors a JOIN base b USING(patid,encounterid) LEFT JOIN deaths d USING(patid) WHERE (b.discharge_date IS NULL OR b.discharge_date > cast(a.anchor_ts+INTERVAL 96 HOUR AS DATE)) AND (d.death_date IS NULL OR d.death_date > cast(a.anchor_ts+INTERVAL 96 HOUR AS DATE))")

    # PRESCRIBING intervals are date-level. Treat missing end date as same-day start for this feasibility audit.
    p_start_expr = f"coalesce(try_cast(p.{qident(p_start)} AS DATE), try_cast(p.{qident(p_order_date)} AS DATE))" if p_start else f"try_cast(p.{qident(p_order_date)} AS DATE)"
    p_end_expr = f"coalesce(try_cast(p.{qident(p_end)} AS DATE), {p_start_expr})" if p_end else p_start_expr
    con.execute(f"""
      CREATE TEMP TABLE p_flags AS
      SELECT c.patid,c.encounterid,
        max(CASE WHEN {p_broad} AND {p_start_expr} <= cast(c.anchor_ts+INTERVAL 96 HOUR AS DATE) AND {p_end_expr} >= cast(c.anchor_ts+INTERVAL 72 HOUR AS DATE) THEN 1 ELSE 0 END) broad_72_96,
        max(CASE WHEN {p_nonbroad} AND {p_start_expr} <= cast(c.anchor_ts+INTERVAL 96 HOUR AS DATE) AND {p_end_expr} >= cast(c.anchor_ts+INTERVAL 72 HOUR AS DATE) THEN 1 ELSE 0 END) nonbroad_72_96
      FROM cohort c LEFT JOIN p ON cast(p.{qident(p_patid)} AS VARCHAR)=c.patid AND cast(p.{qident(p_enc)} AS VARCHAR)=c.encounterid
      GROUP BY 1,2
    """)

    pclass = con.execute("""SELECT CASE WHEN broad_72_96=1 THEN 'continued_broad' WHEN nonbroad_72_96=1 THEN 'deescalated_narrowed_nonbroad' ELSE 'deescalated_stopped_all' END exposure_group, count(*)::BIGINT count FROM p_flags GROUP BY 1 ORDER BY count DESC""").fetchdf()
    pclass["count"] = pclass["count"].map(lambda x: safe(int(x)))
    pclass.to_csv(args.output_dir / "prescribing_exposure_counts.csv", index=False)

    # MED_ADMIN sensitivity uses actual start/stop timestamps where available; if stop absent, use start timestamp.
    m_name_expr = f"lower(coalesce(cast(m.{qident(m_name)} AS VARCHAR),''))"; m_code_expr = f"trim(coalesce(cast(m.{qident(m_code)} AS VARCHAR),''))"
    m_route_expr = f"upper(trim(coalesce(cast(m.{qident(m_route)} AS VARCHAR),'')))" if m_route else "''"
    m_broad = f"(({m_code_expr} IN ({include_sql}) OR regexp_matches({m_name_expr},{q(BROAD_PATTERN)})) AND {m_code_expr} NOT IN ({exclude_sql}) AND NOT regexp_matches({m_name_expr},{q(NON_SYSTEMIC_PATTERN)}) AND {m_route_expr} NOT IN ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    m_nonbroad = f"(regexp_matches({m_name_expr},{q(NON_BROAD_PATTERN)}) AND NOT regexp_matches({m_name_expr},{q(NON_SYSTEMIC_PATTERN)}) AND {m_route_expr} NOT IN ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    m_start_ts = f"try_cast(cast(m.{qident(m_start_date)} AS VARCHAR) || ' ' || coalesce(cast(m.{qident(m_start_time)} AS VARCHAR),'00:00:00') AS TIMESTAMP)" if m_start_time else f"try_cast(m.{qident(m_start_date)} AS TIMESTAMP)"
    if m_stop_date:
        m_stop_ts = f"coalesce(try_cast(cast(m.{qident(m_stop_date)} AS VARCHAR) || ' ' || coalesce(cast(m.{qident(m_stop_time)} AS VARCHAR),'23:59:59') AS TIMESTAMP), {m_start_ts})" if m_stop_time else f"coalesce(try_cast(m.{qident(m_stop_date)} AS TIMESTAMP),{m_start_ts})"
    else:
        m_stop_ts = m_start_ts
    con.execute(f"""
      CREATE TEMP TABLE m_flags AS
      SELECT c.patid,c.encounterid,
        max(CASE WHEN {m_broad} AND {m_start_ts} < c.anchor_ts+INTERVAL 96 HOUR AND {m_stop_ts} >= c.anchor_ts+INTERVAL 72 HOUR THEN 1 ELSE 0 END) broad_72_96,
        max(CASE WHEN {m_nonbroad} AND {m_start_ts} < c.anchor_ts+INTERVAL 96 HOUR AND {m_stop_ts} >= c.anchor_ts+INTERVAL 72 HOUR THEN 1 ELSE 0 END) nonbroad_72_96
      FROM cohort c LEFT JOIN m ON cast(m.{qident(m_patid)} AS VARCHAR)=c.patid AND cast(m.{qident(m_enc)} AS VARCHAR)=c.encounterid
      GROUP BY 1,2
    """)
    mclass = con.execute("""SELECT CASE WHEN broad_72_96=1 THEN 'continued_broad' WHEN nonbroad_72_96=1 THEN 'deescalated_narrowed_nonbroad' ELSE 'deescalated_stopped_all' END exposure_group, count(*)::BIGINT count FROM m_flags GROUP BY 1 ORDER BY count DESC""").fetchdf()
    mclass["count"] = mclass["count"].map(lambda x: safe(int(x)))
    mclass.to_csv(args.output_dir / "med_admin_exposure_counts.csv", index=False)

    reclass = con.execute("""
      SELECT
        CASE WHEN p.broad_72_96=1 THEN 'continued_broad' WHEN p.nonbroad_72_96=1 THEN 'deescalated_narrowed_nonbroad' ELSE 'deescalated_stopped_all' END prescribing_group,
        CASE WHEN m.broad_72_96=1 THEN 'continued_broad' WHEN m.nonbroad_72_96=1 THEN 'deescalated_narrowed_nonbroad' ELSE 'deescalated_stopped_all' END med_admin_group,
        count(*)::BIGINT count
      FROM p_flags p JOIN m_flags m USING(patid,encounterid)
      GROUP BY 1,2 ORDER BY 1,2
    """).fetchdf()
    reclass["count"] = reclass["count"].map(lambda x: safe(int(x)))
    reclass.to_csv(args.output_dir / "source_reclassification_matrix.csv", index=False)

    agree = int(con.execute("""SELECT count(*) FROM p_flags p JOIN m_flags m USING(patid,encounterid) WHERE (CASE WHEN p.broad_72_96=1 THEN 0 WHEN p.nonbroad_72_96=1 THEN 1 ELSE 2 END)=(CASE WHEN m.broad_72_96=1 THEN 0 WHEN m.nonbroad_72_96=1 THEN 1 ELSE 2 END)""").fetchone()[0])
    total = int(con.execute("SELECT count(*) FROM cohort").fetchone()[0])
    pd.DataFrame([{"strict_96h_cohort": safe(total), "source_agreement_count": safe(agree), "source_agreement_proportion": agree/total if total else None}]).to_csv(args.output_dir / "source_agreement_summary.csv", index=False)

    summary = {
      "privacy_mode": "aggregate_only_no_ids_no_patient_rows_no_free_text_export",
      "minimum_reported_cell": MIN_CELL,
      "cohort": "strict modified PSU 96h landmark cohort from hospital-clock feasibility audit",
      "primary_source": "PRESCRIBING",
      "sensitivity_source": "MED_ADMIN",
      "exposure_window": "anchor +72h through +96h",
      "classification": "any broad overlap => continued broad; otherwise any non-broad overlap => narrowed/non-broad; otherwise stopped all",
      "prescribing_interval_caveat": "RX_START_DATE/RX_END_DATE are date-level; missing end date is treated as same-day start for feasibility only",
      "microbiology": "not applied because linked day-3 culture positivity is not faithfully recoverable",
      "guardrail": "Feasibility/reclassification only; do not estimate treatment effects from this audit."
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

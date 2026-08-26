#!/usr/bin/env python3
"""Aggregate-only feasibility audit for a modified PSU external replication cohort.

The PSU extract cannot faithfully reproduce the MIMIC ICU clock or day-3 positive-culture
availability phenotype. This audit therefore evaluates a clearly labeled hospital-clock
modified replication using the frozen systemic broad-spectrum antibiotic proxy in
PRESCRIBING. It exports aggregate counts only and never patient identifiers or rows.
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
NON_SYSTEMIC_PATTERN = (
    "cayston|inhal|nebul|tablet|capsule|oral solution|oral suspension|by mouth|\\bpo\\b"
)


def q(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def find_parquet(root: Path, stem: str) -> Path:
    candidates = sorted(root.glob(f"PCORnet/parquet/**/{stem}*.parquet"))
    if not candidates:
        candidates = sorted(root.glob(f"PCORnet/parquet/**/{stem.replace('_', '')}*.parquet"))
    if not candidates:
        raise FileNotFoundError(f"No parquet found for {stem}")
    exact = [p for p in candidates if p.stem.lower() == stem.lower()]
    return exact[0] if exact else max(candidates, key=lambda p: p.stat().st_size)


def first_present(cols: set[str], names: list[str]) -> str | None:
    lut = {c.upper(): c for c in cols}
    for name in names:
        if name.upper() in lut:
            return lut[name.upper()]
    return None


def qident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


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
                    token = str(item).strip()
                    if re.fullmatch(r"\d{4,9}", token):
                        exclude.add(token)
            elif isinstance(value, (list, tuple)):
                for item in value:
                    if isinstance(item, (list, tuple)) and item:
                        token = str(item[0]).strip()
                        if re.fullmatch(r"\d{4,9}", token):
                            include.add(token)
    return include - exclude, exclude


def safe(v: int) -> int | None:
    return int(v) if int(v) >= MIN_CELL else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    root = args.data_root
    sepsis = find_parquet(root, "sepsis_encounter")
    prescribing = find_parquet(root, "prescribing")
    death = find_parquet(root, "death")
    legacy = root / "PCORnet" / "code" / "config" / "codes_antibiotics.py"
    if not legacy.exists():
        raise FileNotFoundError(f"Missing legacy antibiotic map: {legacy}")

    include_codes, exclude_codes = parse_legacy_codes(legacy)
    include_sql = ",".join(q(x) for x in sorted(include_codes)) or "''"
    exclude_sql = ",".join(q(x) for x in sorted(exclude_codes)) or "''"

    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    con.execute(f"CREATE VIEW s AS SELECT * FROM read_parquet({q(str(sepsis))})")
    con.execute(f"CREATE VIEW p AS SELECT * FROM read_parquet({q(str(prescribing))})")
    con.execute(f"CREATE VIEW d AS SELECT * FROM read_parquet({q(str(death))})")

    s_cols = set(con.execute("DESCRIBE s").fetchdf()["column_name"].astype(str))
    p_cols = set(con.execute("DESCRIBE p").fetchdf()["column_name"].astype(str))
    d_cols = set(con.execute("DESCRIBE d").fetchdf()["column_name"].astype(str))

    s_patid = first_present(s_cols, ["PATID"])
    s_enc = first_present(s_cols, ["ENCOUNTERID", "ENCOUNTER_ID"])
    admit = first_present(s_cols, ["ADMIT_DATE", "ADMITDATE"])
    discharge = first_present(s_cols, ["DISCHARGE_DATE", "DISCHARGEDATE"])
    p_patid = first_present(p_cols, ["PATID"])
    p_enc = first_present(p_cols, ["ENCOUNTERID", "ENCOUNTER_ID"])
    order_date = first_present(p_cols, ["RX_ORDER_DATE", "ORDER_DATE"])
    order_time = first_present(p_cols, ["RX_ORDER_TIME", "ORDER_TIME"])
    med_name = first_present(p_cols, ["RAW_RX_MED_NAME", "RX_MED_NAME", "MEDICATION_NAME"])
    rxnorm = first_present(p_cols, ["RXNORM_CUI", "RAW_RXNORM_CUI", "RXNORM"])
    route = first_present(p_cols, ["RX_ROUTE", "RAW_RX_ROUTE", "ROUTE"])
    d_patid = first_present(d_cols, ["PATID"])
    death_date = first_present(d_cols, ["DEATH_DATE", "DEATHDATE"])

    required = {
        "sepsis_patid": s_patid,
        "sepsis_encounter": s_enc,
        "admit_date": admit,
        "discharge_date": discharge,
        "prescribing_patid": p_patid,
        "prescribing_encounter": p_enc,
        "rx_order_date": order_date,
        "med_name": med_name,
        "rxnorm": rxnorm,
        "death_patid": d_patid,
        "death_date": death_date,
    }
    missing = [k for k, v in required.items() if v is None]
    if missing:
        raise RuntimeError("Required columns missing: " + ", ".join(missing))

    if order_time:
        order_ts = (
            f"try_cast(cast(p.{qident(order_date)} AS VARCHAR) || ' ' || "
            f"coalesce(cast(p.{qident(order_time)} AS VARCHAR),'00:00:00') AS TIMESTAMP)"
        )
    else:
        order_ts = f"try_cast(p.{qident(order_date)} AS TIMESTAMP)"

    name_expr = f"lower(coalesce(cast(p.{qident(med_name)} AS VARCHAR),''))"
    code_expr = f"trim(coalesce(cast(p.{qident(rxnorm)} AS VARCHAR),''))"
    route_expr = f"upper(trim(coalesce(cast(p.{qident(route)} AS VARCHAR),'')))" if route else "''"
    systemic = (
        f"(({code_expr} IN ({include_sql}) OR regexp_matches({name_expr}, {q(BROAD_PATTERN)})) "
        f"AND {code_expr} NOT IN ({exclude_sql}) "
        f"AND NOT regexp_matches({name_expr}, {q(NON_SYSTEMIC_PATTERN)}) "
        f"AND {route_expr} NOT IN ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    )

    con.execute(f"""
        CREATE TEMP TABLE base AS
        SELECT DISTINCT
          cast(s.{qident(s_patid)} AS VARCHAR) AS patid,
          cast(s.{qident(s_enc)} AS VARCHAR) AS encounterid,
          try_cast(s.{qident(admit)} AS DATE) AS admit_date,
          try_cast(s.{qident(discharge)} AS DATE) AS discharge_date
        FROM s
        WHERE s.{qident(s_patid)} IS NOT NULL AND s.{qident(s_enc)} IS NOT NULL
    """)

    con.execute(f"""
        CREATE TEMP TABLE qualifying_orders AS
        SELECT
          cast(p.{qident(p_patid)} AS VARCHAR) AS patid,
          cast(p.{qident(p_enc)} AS VARCHAR) AS encounterid,
          {order_ts} AS order_ts
        FROM p
        INNER JOIN base b
          ON cast(p.{qident(p_patid)} AS VARCHAR)=b.patid
         AND cast(p.{qident(p_enc)} AS VARCHAR)=b.encounterid
        WHERE {systemic}
          AND {order_ts} IS NOT NULL
          AND b.admit_date IS NOT NULL
          AND {order_ts} >= cast(b.admit_date AS TIMESTAMP)
          AND {order_ts} < cast(b.admit_date AS TIMESTAMP) + INTERVAL 24 HOUR
    """)

    con.execute("""
        CREATE TEMP TABLE anchors AS
        SELECT patid, encounterid, min(order_ts) AS anchor_ts
        FROM qualifying_orders
        GROUP BY 1,2
    """)

    con.execute(f"""
        CREATE TEMP TABLE death_by_patient AS
        SELECT cast({qident(d_patid)} AS VARCHAR) AS patid,
               min(try_cast({qident(death_date)} AS DATE)) AS death_date
        FROM d
        WHERE {qident(d_patid)} IS NOT NULL
        GROUP BY 1
    """)

    counts = []
    def add(step: str, sql: str) -> None:
        n = int(con.execute(sql).fetchone()[0])
        counts.append({"step": step, "count": safe(n)})

    add("sepsis_encounters_total", "SELECT count(*) FROM base")
    add("sepsis_encounters_with_admit_date", "SELECT count(*) FROM base WHERE admit_date IS NOT NULL")
    add("sepsis_encounters_with_discharge_date", "SELECT count(*) FROM base WHERE discharge_date IS NOT NULL")
    add("encounters_with_systemic_broad_order_in_first_24h", "SELECT count(*) FROM anchors")
    add(
        "anchors_with_discharge_after_day3_date",
        """SELECT count(*) FROM anchors a JOIN base b USING(patid,encounterid)
           WHERE b.discharge_date IS NULL OR b.discharge_date > cast(a.anchor_ts + INTERVAL 72 HOUR AS DATE)""",
    )
    add(
        "strict_96h_landmark_discharge_only",
        """SELECT count(*) FROM anchors a JOIN base b USING(patid,encounterid)
           WHERE b.discharge_date IS NULL OR b.discharge_date > cast(a.anchor_ts + INTERVAL 96 HOUR AS DATE)""",
    )
    add(
        "lenient_96h_landmark_discharge_only",
        """SELECT count(*) FROM anchors a JOIN base b USING(patid,encounterid)
           WHERE b.discharge_date IS NULL OR b.discharge_date >= cast(a.anchor_ts + INTERVAL 96 HOUR AS DATE)""",
    )
    add(
        "strict_96h_landmark_discharge_and_death",
        """SELECT count(*) FROM anchors a JOIN base b USING(patid,encounterid)
           LEFT JOIN death_by_patient d USING(patid)
           WHERE (b.discharge_date IS NULL OR b.discharge_date > cast(a.anchor_ts + INTERVAL 96 HOUR AS DATE))
             AND (d.death_date IS NULL OR d.death_date > cast(a.anchor_ts + INTERVAL 96 HOUR AS DATE))""",
    )
    add(
        "lenient_96h_landmark_discharge_and_death",
        """SELECT count(*) FROM anchors a JOIN base b USING(patid,encounterid)
           LEFT JOIN death_by_patient d USING(patid)
           WHERE (b.discharge_date IS NULL OR b.discharge_date >= cast(a.anchor_ts + INTERVAL 96 HOUR AS DATE))
             AND (d.death_date IS NULL OR d.death_date >= cast(a.anchor_ts + INTERVAL 96 HOUR AS DATE))""",
    )

    pd.DataFrame(counts).to_csv(args.output_dir / "cohort_step_counts.csv", index=False)

    timing = con.execute(f"""
        SELECT
          count(*)::BIGINT AS prescribing_rows,
          count(*) FILTER (WHERE {qident(order_date)} IS NOT NULL)::BIGINT AS order_date_present,
          {('count(*) FILTER (WHERE ' + qident(order_time) + ' IS NOT NULL)::BIGINT') if order_time else '0::BIGINT'} AS order_time_present,
          count(*) FILTER (WHERE {qident(med_name)} IS NOT NULL)::BIGINT AS med_name_present,
          count(*) FILTER (WHERE {qident(rxnorm)} IS NOT NULL)::BIGINT AS rxnorm_present
        FROM p
    """).fetchdf()
    timing.to_csv(args.output_dir / "prescribing_time_completeness.csv", index=False)

    landmark_bounds = con.execute("""
        SELECT
          count(*)::BIGINT AS anchored_encounters,
          count(*) FILTER (WHERE b.discharge_date = cast(a.anchor_ts + INTERVAL 96 HOUR AS DATE))::BIGINT AS discharge_on_landmark_calendar_date,
          count(*) FILTER (WHERE d.death_date = cast(a.anchor_ts + INTERVAL 96 HOUR AS DATE))::BIGINT AS death_on_landmark_calendar_date
        FROM anchors a
        JOIN base b USING(patid,encounterid)
        LEFT JOIN death_by_patient d USING(patid)
    """).fetchdf()
    landmark_bounds.to_csv(args.output_dir / "landmark_date_ambiguity.csv", index=False)

    field_map = pd.DataFrame([{"construct": k, "column": v or "<NONE>"} for k, v in {
        **required,
        "rx_order_time": order_time,
        "route": route,
    }.items()])
    field_map.to_csv(args.output_dir / "field_map.csv", index=False)

    summary = {
        "privacy_mode": "aggregate_only_no_ids_no_patient_rows_no_free_text_export",
        "minimum_reported_cell": MIN_CELL,
        "clock": "modified external replication: hospital admission calendar date at midnight; exact ICU entry unavailable",
        "anchor": "first systemic broad-spectrum PRESCRIBING order within 24h after hospital admission date",
        "antibiotic_phenotype": "reviewed legacy RxNorm inclusion plus frozen broad-spectrum text fallback; explicit oral/inhaled exclusions; route field not required because predominantly UN",
        "microbiology": "not applied in this feasibility audit because linked day-3 culture positivity is not faithfully recoverable from the current PSU PCORnet extract",
        "landmark": "anchor +96h; strict and lenient bounds reported because discharge/death are date-level in available sources",
        "guardrail": "Feasibility only. Do not estimate treatment effects or call this an exact MIMIC replication from this audit.",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

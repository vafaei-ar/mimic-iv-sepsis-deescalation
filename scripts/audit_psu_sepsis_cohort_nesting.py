#!/usr/bin/env python3
"""Confirm that the frozen strict PSU analytic cohort is nested in sepsis_encounter.

This is a provenance audit only. It reconstructs the frozen strict cohort from the
local restricted PSU/PCORnet sources and reports aggregate counts showing whether any
analytic encounters fall outside the upstream sepsis_encounter source. It does not
validate how the upstream adapted sepsis definition itself was implemented and never
exports patient- or encounter-level rows.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb

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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    root = args.data_root
    sepsis_path = find(root, ["sepsis_encounter"])
    prescribing_path = find(root, ["prescribing"])
    death_path = find(root, ["death"])

    con = duckdb.connect()
    con.execute("pragma threads=4")
    con.execute(f"create view s as select * from read_parquet({q(str(sepsis_path))})")
    con.execute(f"create view p as select * from read_parquet({q(str(prescribing_path))})")
    con.execute(f"create view d as select * from read_parquet({q(str(death_path))})")

    sc = set(con.execute("describe s").fetchdf()["column_name"].astype(str))
    pc = set(con.execute("describe p").fetchdf()["column_name"].astype(str))
    dc = set(con.execute("describe d").fetchdf()["column_name"].astype(str))

    sp, se = first(sc, ["PATID"]), first(sc, ["ENCOUNTERID", "ENCOUNTER_ID"])
    ad, dis = first(sc, ["ADMIT_DATE"]), first(sc, ["DISCHARGE_DATE"])
    pp, pe = first(pc, ["PATID"]), first(pc, ["ENCOUNTERID", "ENCOUNTER_ID"])
    pod, pot = first(pc, ["RX_ORDER_DATE"]), first(pc, ["RX_ORDER_TIME"])
    pn, pcode, proute = first(pc, ["RAW_RX_MED_NAME"]), first(pc, ["RXNORM_CUI"]), first(pc, ["RX_ROUTE"])
    dp, dd = first(dc, ["PATID"]), first(dc, ["DEATH_DATE"])

    required = {"sepsis PATID": sp, "sepsis ENCOUNTERID": se, "ADMIT_DATE": ad, "DISCHARGE_DATE": dis,
                "prescribing PATID": pp, "prescribing ENCOUNTERID": pe, "RX_ORDER_DATE": pod,
                "RAW_RX_MED_NAME": pn, "RXNORM_CUI": pcode, "death PATID": dp, "DEATH_DATE": dd}
    missing = [name for name, value in required.items() if value is None]
    if missing:
        raise RuntimeError(f"Required PSU fields unavailable: {missing}")

    legacy = root / "PCORnet" / "code" / "config" / "codes_antibiotics.py"
    include_codes, exclude_codes = parse_legacy_codes(legacy)
    incsql = ",".join(q(x) for x in sorted(include_codes)) or "''"
    excsql = ",".join(q(x) for x in sorted(exclude_codes)) or "''"

    pts = ts("p", pod, pot)
    pname = f"lower(coalesce(cast(p.{qi(pn)} as varchar),''))"
    pcode_expr = f"trim(coalesce(cast(p.{qi(pcode)} as varchar),''))"
    route = f"upper(trim(coalesce(cast(p.{qi(proute)} as varchar),'')))" if proute else "''"
    broad = (
        f"(({pcode_expr} in ({incsql}) or regexp_matches({pname},{q(BROAD_PATTERN)})) "
        f"and {pcode_expr} not in ({excsql}) and not regexp_matches({pname},{q(NON_SYSTEMIC_PATTERN)}) "
        f"and {route} not in ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    )

    con.execute(
        f"create temp table sepsis_ids as "
        f"select distinct cast({qi(sp)} as varchar) patid, cast({qi(se)} as varchar) encounterid, "
        f"try_cast({qi(ad)} as date) admit_date, try_cast({qi(dis)} as date) discharge_date from s"
    )
    con.execute(
        f"create temp table anchors as "
        f"select cast(p.{qi(pp)} as varchar) patid, cast(p.{qi(pe)} as varchar) encounterid, min({pts}) anchor_ts "
        f"from p join sepsis_ids s on cast(p.{qi(pp)} as varchar)=s.patid and cast(p.{qi(pe)} as varchar)=s.encounterid "
        f"where {broad} and {pts}>=cast(s.admit_date as timestamp) "
        f"and {pts}<cast(s.admit_date as timestamp)+interval 24 hour group by 1,2"
    )
    con.execute(
        f"create temp table deaths as select cast({qi(dp)} as varchar) patid, "
        f"min(try_cast({qi(dd)} as date)) death_date from d group by 1"
    )
    con.execute(
        "create temp table strict_cohort as "
        "select a.patid,a.encounterid,a.anchor_ts,s.admit_date,s.discharge_date,d.death_date "
        "from anchors a join sepsis_ids s using(patid,encounterid) left join deaths d using(patid) "
        "where (s.discharge_date is null or s.discharge_date>cast(a.anchor_ts+interval 96 hour as date)) "
        "and (d.death_date is null or d.death_date>cast(a.anchor_ts+interval 96 hour as date))"
    )

    source_n = int(con.execute("select count(*) from sepsis_ids").fetchone()[0])
    strict_n = int(con.execute("select count(*) from strict_cohort").fetchone()[0])
    unmatched_n = int(
        con.execute(
            "select count(*) from strict_cohort c left join sepsis_ids s using(patid,encounterid) "
            "where s.patid is null"
        ).fetchone()[0]
    )
    distinct_strict_n = int(con.execute("select count(*) from (select distinct patid,encounterid from strict_cohort)").fetchone()[0])

    result = {
        "audit": "PSU frozen strict cohort nesting in upstream sepsis_encounter",
        "upstream_sepsis_encounters": source_n,
        "strict_analytic_cohort_rows": strict_n,
        "strict_analytic_distinct_encounters": distinct_strict_n,
        "analytic_encounters_not_in_sepsis_encounter": unmatched_n,
        "all_analytic_encounters_nested": unmatched_n == 0 and strict_n == distinct_strict_n,
        "expected_frozen_strict_cohort_n": 19841,
        "matches_frozen_strict_cohort_n": strict_n == 19841,
        "interpretation": (
            "This confirms provenance/nesting only. It does not independently validate the upstream implementation "
            "of the adapted sepsis definition (suspected/confirmed infection plus absolute modified SOFA >=2 without GCS)."
        ),
        "data_safety": "Aggregate counts only; no patient- or encounter-level data are exported.",
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))

    if not result["all_analytic_encounters_nested"] or not result["matches_frozen_strict_cohort_n"]:
        raise SystemExit("PSU sepsis cohort nesting audit failed; inspect output JSON")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Aggregate-only audit of FACILITYID values in PSU/PCORnet parquet sources.

The PSU source extract used by the publication pipeline is stored as parquet under
PCORnet/parquet. This audit scans only parquet schemas, selects encounter/admission-like
files containing FACILITYID, and writes aggregate counts only. It never emits patient-
or encounter-level rows.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb


def q(v: str) -> str:
    return "'" + v.replace("'", "''") + "'"


def qi(v: str) -> str:
    return '"' + v.replace('"', '""') + '"'


def parquet_candidates(root: Path) -> list[Path]:
    base = root / "PCORnet" / "parquet"
    if not base.exists():
        base = root
    return sorted(p.resolve() for p in base.glob("**/*.parquet") if p.is_file())


def columns(con: duckdb.DuckDBPyConnection, path: Path) -> list[str]:
    return list(
        con.execute(f"DESCRIBE SELECT * FROM read_parquet({q(str(path))})")
        .fetchdf()["column_name"]
        .astype(str)
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    con = duckdb.connect()
    con.execute("PRAGMA threads=4")

    all_parquet = parquet_candidates(args.data_root)
    schema_matches: list[tuple[Path, list[str]]] = []
    for path in all_parquet:
        try:
            cols = columns(con, path)
        except Exception:
            continue
        lut = {c.upper(): c for c in cols}
        if "FACILITYID" not in lut:
            continue
        name = path.stem.lower().replace("_", "")
        if any(token in name for token in ("admission", "encounter", "sepsis")):
            schema_matches.append((path, cols))

    # If the extract uses a nonstandard filename, fall back to any parquet table that
    # contains FACILITYID. The output records the source filename so interpretation is explicit.
    if not schema_matches:
        for path in all_parquet:
            try:
                cols = columns(con, path)
            except Exception:
                continue
            if "FACILITYID" in {c.upper() for c in cols}:
                schema_matches.append((path, cols))

    if not schema_matches:
        raise RuntimeError("No parquet source containing FACILITYID was found under the PSU data root")

    summaries = []
    for path, cols in schema_matches:
        lut = {c.upper(): c for c in cols}
        fac = lut["FACILITYID"]
        pat = lut.get("PATID")
        enc = lut.get("ENCOUNTERID") or lut.get("ENCOUNTER_ID")

        select_parts = [
            f"coalesce(nullif(trim(cast({qi(fac)} AS VARCHAR)), ''), '<MISSING>') AS facility_id",
            "count(*)::BIGINT AS rows",
        ]
        if pat:
            select_parts.append(f"count(DISTINCT cast({qi(pat)} AS VARCHAR))::BIGINT AS unique_patients")
        if enc:
            select_parts.append(f"count(DISTINCT cast({qi(enc)} AS VARCHAR))::BIGINT AS unique_encounters")

        df = con.execute(
            "SELECT " + ", ".join(select_parts)
            + f" FROM read_parquet({q(str(path))}) GROUP BY 1 ORDER BY rows DESC, facility_id"
        ).fetchdf()

        facilities = []
        for row in df.to_dict(orient="records"):
            item = {"facility_id": str(row["facility_id"]), "rows": int(row["rows"])}
            if "unique_patients" in row:
                item["unique_patients"] = int(row["unique_patients"])
            if "unique_encounters" in row:
                item["unique_encounters"] = int(row["unique_encounters"])
            facilities.append(item)

        summaries.append(
            {
                "source_file": str(path.relative_to(args.data_root)),
                "n_facility_ids": len(facilities),
                "facilities": facilities,
            }
        )

    result = {
        "audit": "PSU FACILITYID aggregate audit",
        "matching_sources": len(summaries),
        "sources": summaries,
        "note": "Aggregate counts only; no patient- or encounter-level data are included.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

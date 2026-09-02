#!/usr/bin/env python3
"""Aggregate-only audit of FACILITYID values in the PSU ADMISSION table.

This script reads the local restricted Penn State/PCORnet source but writes only
sanitized aggregate facility counts. It never emits patient- or encounter-level rows.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def _candidate_files(root: Path) -> list[Path]:
    out: list[Path] = []
    for pattern in ("**/*ADMISSION*.csv", "**/*admission*.csv", "**/*ADMISSION*.tsv", "**/*admission*.tsv", "**/*ADMISSION*.txt", "**/*admission*.txt"):
        out.extend(root.glob(pattern))
    return sorted({p.resolve() for p in out if p.is_file()})


def _read_header(path: Path) -> tuple[str, list[str]]:
    for sep in ("|", "\t", ","):
        try:
            cols = list(pd.read_csv(path, sep=sep, nrows=0, dtype=str).columns)
        except Exception:
            continue
        if len(cols) > 1:
            return sep, cols
    raise RuntimeError(f"Could not determine delimiter for {path}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    candidates = _candidate_files(args.data_root)
    matches: list[tuple[Path, str, list[str]]] = []
    for path in candidates:
        try:
            sep, cols = _read_header(path)
        except Exception:
            continue
        lookup = {c.upper(): c for c in cols}
        if "FACILITYID" in lookup:
            matches.append((path, sep, cols))

    if not matches:
        raise RuntimeError("No ADMISSION-like source file containing FACILITYID was found")

    summaries = []
    for path, sep, cols in matches:
        lookup = {c.upper(): c for c in cols}
        usecols = [lookup["FACILITYID"]]
        if "PATID" in lookup:
            usecols.append(lookup["PATID"])
        df = pd.read_csv(path, sep=sep, usecols=usecols, dtype=str, low_memory=False)
        facility_col = lookup["FACILITYID"]
        df[facility_col] = df[facility_col].fillna("<MISSING>").astype(str).str.strip()
        rows = []
        for facility_id, g in df.groupby(facility_col, dropna=False):
            item = {
                "facility_id": str(facility_id),
                "admission_rows": int(len(g)),
            }
            if "PATID" in lookup:
                item["unique_patients"] = int(g[lookup["PATID"]].nunique(dropna=True))
            rows.append(item)
        summaries.append(
            {
                "source_file": str(path.relative_to(args.data_root)),
                "n_rows": int(len(df)),
                "n_facility_ids": int(df[facility_col].nunique(dropna=False)),
                "facilities": sorted(rows, key=lambda x: (-x["admission_rows"], x["facility_id"])),
            }
        )

    result = {
        "audit": "PSU ADMISSION FACILITYID aggregate audit",
        "data_root": str(args.data_root),
        "matching_admission_files": len(summaries),
        "sources": summaries,
        "note": "Aggregate counts only; no patient- or encounter-level data are included.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Aggregate-only PSU covariate mapping audit before PS modeling.

Produces candidate code/name/unit maps for labs, OBS_CLIN physiologic variables,
vasopressors, and baseline demographic/diagnosis fields. No patient-level rows,
identifiers, free-text result values, propensity scores, or treatment effects are exported.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb
import pandas as pd

MIN_CELL = 11
LAB_PATTERNS = {
    "lactate": r"lactate|lactic acid",
    "creatinine": r"creatinine",
    "wbc": r"white blood|\\bwbc\\b",
    "platelet": r"platelet",
    "bilirubin": r"bilirubin",
}
OBS_PATTERNS = {
    "heart_rate": r"^heart rate$|pulse rate",
    "resp_rate": r"respiratory rate|resp rate",
    "spo2": r"spo2|oxygen saturation|o2 saturation",
    "map": r"mean arterial|arterial pressure mean|blood pressure mean|\\bmap\\b",
    "temperature": r"temperature|^temp$",
    "gcs": r"glasgow|\\bgcs\\b",
    "fio2": r"fio2|fraction.*inspired.*oxygen",
}
VASO_PATTERN = r"norepinephrine|levophed|phenylephrine|vasopressin|epinephrine|dopamine"


def q(v: str) -> str:
    return "'" + v.replace("'", "''") + "'"


def find_parquet(root: Path, stems: list[str], required: bool = True) -> Path | None:
    for stem in stems:
        candidates = sorted(root.glob(f"PCORnet/parquet/**/{stem}*.parquet"))
        if not candidates:
            candidates = sorted(root.glob(f"PCORnet/parquet/**/{stem.replace('_', '')}*.parquet"))
        if candidates:
            exact = [p for p in candidates if p.stem.lower() == stem.lower()]
            return exact[0] if exact else max(candidates, key=lambda p: p.stat().st_size)
    if required:
        raise FileNotFoundError(f"No parquet found for {stems}")
    return None


def suppress(df: pd.DataFrame, count_cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in count_cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
            out.loc[out[col] < MIN_CELL, col] = pd.NA
    return out


def describe(con: duckdb.DuckDBPyConnection, view: str) -> set[str]:
    return set(con.execute(f"DESCRIBE {view}").fetchdf()["column_name"].astype(str))


def first_present(cols: set[str], names: list[str]) -> str | None:
    lut = {c.upper(): c for c in cols}
    for n in names:
        if n.upper() in lut:
            return lut[n.upper()]
    return None


def qi(v: str) -> str:
    return '"' + v.replace('"', '""') + '"'


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    root = args.data_root

    paths = {
        "sepsis": find_parquet(root, ["sepsis_encounter"]),
        "lab": find_parquet(root, ["lab_reduced", "lab_result_cm"]),
        "obs": find_parquet(root, ["obs_clin"]),
        "med": find_parquet(root, ["med_admin"]),
        "dem": find_parquet(root, ["sepsis_demographic", "demographic"]),
        "dx": find_parquet(root, ["sepsis_diagnosis", "diagnosis"]),
    }

    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    for name, path in paths.items():
        con.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet({q(str(path))})")

    # Restrict maps to the known sepsis encounter universe to avoid unrelated-system metadata.
    sc = describe(con, "sepsis")
    s_pat = first_present(sc, ["PATID"])
    s_enc = first_present(sc, ["ENCOUNTERID", "ENCOUNTER_ID"])
    if not s_pat or not s_enc:
        raise RuntimeError("sepsis encounter identifiers unavailable")
    con.execute(
        f"CREATE TEMP TABLE sepsis_ids AS SELECT DISTINCT cast({qi(s_pat)} AS VARCHAR) patid, "
        f"cast({qi(s_enc)} AS VARCHAR) encounterid FROM sepsis"
    )

    # Lab candidate map: concept, LOINC, raw lab name, unit, aggregate rows and encounters.
    lc = describe(con, "lab")
    l_pat = first_present(lc, ["PATID"])
    l_enc = first_present(lc, ["ENCOUNTERID", "ENCOUNTER_ID"])
    l_name = first_present(lc, ["RAW_LAB_NAME", "LAB_NAME"])
    l_loinc = first_present(lc, ["LAB_LOINC", "LOINC"])
    l_unit = first_present(lc, ["RESULT_UNIT", "RAW_UNIT", "LAB_RESULT_UNIT"])
    l_num = first_present(lc, ["RESULT_NUM", "LAB_RESULT_NUM", "RESULT_NUMERIC"])
    if not all([l_pat, l_enc, l_name, l_loinc, l_unit, l_num]):
        raise RuntimeError("Required lab mapping fields unavailable")
    lab_frames = []
    for concept, pattern in LAB_PATTERNS.items():
        df = con.execute(f"""
            SELECT {q(concept)} AS concept,
                   cast(l.{qi(l_loinc)} AS VARCHAR) AS code,
                   cast(l.{qi(l_name)} AS VARCHAR) AS raw_name,
                   cast(l.{qi(l_unit)} AS VARCHAR) AS unit,
                   count(*)::BIGINT AS rows,
                   count(DISTINCT cast(l.{qi(l_enc)} AS VARCHAR))::BIGINT AS encounters,
                   count(*) FILTER (WHERE try_cast(l.{qi(l_num)} AS DOUBLE) IS NOT NULL)::BIGINT AS numeric_rows
            FROM lab l
            JOIN sepsis_ids s ON cast(l.{qi(l_pat)} AS VARCHAR)=s.patid
                              AND cast(l.{qi(l_enc)} AS VARCHAR)=s.encounterid
            WHERE regexp_matches(lower(coalesce(cast(l.{qi(l_name)} AS VARCHAR),'')), {q(pattern)})
            GROUP BY 1,2,3,4
            HAVING count(*) >= {MIN_CELL}
            ORDER BY encounters DESC, rows DESC
            LIMIT 100
        """).fetchdf()
        lab_frames.append(df)
    lab_map = pd.concat(lab_frames, ignore_index=True) if lab_frames else pd.DataFrame()
    suppress(lab_map, ["rows", "encounters", "numeric_rows"]).to_csv(args.output_dir / "lab_candidate_map.csv", index=False)

    # OBS_CLIN map: code/name/unit plus numeric coverage. This is where HR/RR/SpO2/MAP/temp are frozen.
    oc = describe(con, "obs")
    o_pat = first_present(oc, ["PATID"])
    o_enc = first_present(oc, ["ENCOUNTERID", "ENCOUNTER_ID"])
    o_name = first_present(oc, ["RAW_OBSCLIN_NAME", "OBSCLIN_NAME", "RAW_OBS_NAME", "OBS_NAME"])
    o_code = first_present(oc, ["OBSCLIN_CODE", "RAW_OBSCLIN_CODE", "OBS_CODE"])
    o_unit = first_present(oc, ["OBSCLIN_RESULT_UNIT", "RAW_OBSCLIN_UNIT", "OBS_RESULT_UNIT"])
    o_num = first_present(oc, ["OBSCLIN_RESULT_NUM", "OBS_RESULT_NUM", "RESULT_NUM"])
    if not all([o_pat, o_enc, o_name, o_code, o_unit, o_num]):
        raise RuntimeError("Required OBS_CLIN mapping fields unavailable")
    obs_frames = []
    for concept, pattern in OBS_PATTERNS.items():
        df = con.execute(f"""
            SELECT {q(concept)} AS concept,
                   cast(o.{qi(o_code)} AS VARCHAR) AS code,
                   cast(o.{qi(o_name)} AS VARCHAR) AS raw_name,
                   cast(o.{qi(o_unit)} AS VARCHAR) AS unit,
                   count(*)::BIGINT AS rows,
                   count(DISTINCT cast(o.{qi(o_enc)} AS VARCHAR))::BIGINT AS encounters,
                   count(*) FILTER (WHERE try_cast(o.{qi(o_num)} AS DOUBLE) IS NOT NULL)::BIGINT AS numeric_rows
            FROM obs o
            JOIN sepsis_ids s ON cast(o.{qi(o_pat)} AS VARCHAR)=s.patid
                              AND cast(o.{qi(o_enc)} AS VARCHAR)=s.encounterid
            WHERE regexp_matches(lower(coalesce(cast(o.{qi(o_name)} AS VARCHAR),'')), {q(pattern)})
            GROUP BY 1,2,3,4
            HAVING count(*) >= {MIN_CELL}
            ORDER BY encounters DESC, rows DESC
            LIMIT 100
        """).fetchdf()
        obs_frames.append(df)
    obs_map = pd.concat(obs_frames, ignore_index=True) if obs_frames else pd.DataFrame()
    suppress(obs_map, ["rows", "encounters", "numeric_rows"]).to_csv(args.output_dir / "obsclin_candidate_map.csv", index=False)

    # Vasopressor candidate map from MED_ADMIN.
    mc = describe(con, "med")
    m_pat = first_present(mc, ["PATID"])
    m_enc = first_present(mc, ["ENCOUNTERID", "ENCOUNTER_ID"])
    m_name = first_present(mc, ["RAW_MEDADMIN_MED_NAME", "MEDADMIN_MED_NAME", "MEDICATION_NAME"])
    m_code = first_present(mc, ["MEDADMIN_CODE", "RAW_MEDADMIN_CODE", "RXNORM_CUI"])
    m_route = first_present(mc, ["MEDADMIN_ROUTE", "RAW_MEDADMIN_ROUTE", "ROUTE"])
    if not all([m_pat, m_enc, m_name, m_code, m_route]):
        raise RuntimeError("Required MED_ADMIN mapping fields unavailable")
    vaso = con.execute(f"""
        SELECT cast(m.{qi(m_code)} AS VARCHAR) AS code,
               cast(m.{qi(m_name)} AS VARCHAR) AS raw_name,
               cast(m.{qi(m_route)} AS VARCHAR) AS route,
               count(*)::BIGINT AS rows,
               count(DISTINCT cast(m.{qi(m_enc)} AS VARCHAR))::BIGINT AS encounters
        FROM med m
        JOIN sepsis_ids s ON cast(m.{qi(m_pat)} AS VARCHAR)=s.patid
                          AND cast(m.{qi(m_enc)} AS VARCHAR)=s.encounterid
        WHERE regexp_matches(lower(coalesce(cast(m.{qi(m_name)} AS VARCHAR),'')), {q(VASO_PATTERN)})
        GROUP BY 1,2,3
        HAVING count(*) >= {MIN_CELL}
        ORDER BY encounters DESC, rows DESC
        LIMIT 150
    """).fetchdf()
    suppress(vaso, ["rows", "encounters"]).to_csv(args.output_dir / "vasopressor_candidate_map.csv", index=False)

    # Baseline field completeness only. No identifiers or patient rows are exported.
    dc = describe(con, "dem")
    d_pat = first_present(dc, ["PATID"])
    baseline_rows = []
    for construct, candidates in {
        "birth_date": ["BIRTH_DATE"],
        "sex": ["SEX", "RAW_SEX"],
        "race": ["RACE", "RAW_RACE"],
        "hispanic": ["HISPANIC", "RAW_HISPANIC"],
        "ruca_code": ["RUCA_CODE"],
    }.items():
        col = first_present(dc, candidates)
        if col and d_pat:
            n = int(con.execute(f"SELECT count(DISTINCT cast({qi(d_pat)} AS VARCHAR)) FROM dem WHERE {qi(col)} IS NOT NULL").fetchone()[0])
            baseline_rows.append({"construct": construct, "source": "DEMOGRAPHIC", "field": col, "nonmissing_patients": n if n >= MIN_CELL else None})
        else:
            baseline_rows.append({"construct": construct, "source": "DEMOGRAPHIC", "field": None, "nonmissing_patients": None})

    xc = describe(con, "dx")
    for construct, candidates in {
        "diagnosis_code": ["DX"],
        "raw_diagnosis_code": ["RAW_DX"],
        "principal_diagnosis_flag": ["PDX", "RAW_PDX"],
        "diagnosis_date": ["DX_DATE"],
    }.items():
        col = first_present(xc, candidates)
        baseline_rows.append({"construct": construct, "source": "DIAGNOSIS", "field": col, "nonmissing_patients": None})
    pd.DataFrame(baseline_rows).to_csv(args.output_dir / "baseline_field_map.csv", index=False)

    # Prespecified harmonization recommendations, before any treatment-effect model is run.
    recommendations = pd.DataFrame([
        {"construct": "age", "status": "retain", "rule": "derive from BIRTH_DATE and hospital admission date"},
        {"construct": "sex", "status": "retain", "rule": "use PCORnet SEX; report missing/other descriptively"},
        {"construct": "race", "status": "retain_modified", "rule": "harmonize only categories supported by PSU coding; do not force unavailable MIMIC categories"},
        {"construct": "comorbidity", "status": "retain_modified", "rule": "freeze ICD code sets before modeling; use diagnosis history available in PSU"},
        {"construct": "lactate_creatinine_wbc_platelet_bilirubin", "status": "retain_if_map_validated", "rule": "use frozen code/name/unit map from this audit plus MIMIC plausibility ranges"},
        {"construct": "heart_rate_resp_rate_spo2_map_temperature", "status": "retain_if_map_validated", "rule": "use exact OBS_CLIN code/unit map from this audit and prespecified plausibility ranges"},
        {"construct": "vasopressors", "status": "retain", "rule": "use MED_ADMIN timestamps with frozen vasopressor agent map"},
        {"construct": "gcs", "status": "exclude_from_primary", "rule": "do not substitute opportunistically; retain only if a clearly valid code/value map emerges"},
        {"construct": "fio2", "status": "exclude_from_primary", "rule": "coverage is too sparse for primary harmonized model; may remain descriptive/sensitivity only"},
        {"construct": "mechanical_ventilation_procedure", "status": "unavailable_current_extract", "rule": "do not fabricate a proxy unless a validated source is identified prospectively"},
        {"construct": "ruca_code", "status": "descriptive_psu_specific", "rule": "do not silently add to harmonized primary PS model"},
    ])
    recommendations.to_csv(args.output_dir / "harmonization_recommendations.csv", index=False)

    summary = {
        "privacy_mode": "aggregate_only_no_ids_no_patient_rows_no_result_free_text",
        "minimum_reported_cell": MIN_CELL,
        "purpose": "freeze exact PSU covariate mappings before propensity-score balance diagnostics",
        "sources": ["LAB_RESULT_CM", "OBS_CLIN", "MED_ADMIN", "DEMOGRAPHIC", "DIAGNOSIS"],
        "guardrail": "No propensity score or treatment effect is estimated by this audit. Freeze mappings before outcome-aware modeling.",
        "primary_harmonization": "retain only covariates with defensible PSU code/unit semantics; explicitly exclude non-harmonizable constructs rather than replacing them post hoc",
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

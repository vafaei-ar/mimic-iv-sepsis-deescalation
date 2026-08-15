from __future__ import annotations

import gzip
import zipfile
from pathlib import Path
from typing import Callable, Iterable

import pandas as pd


def _resolve_member_name(zip_path: Path, rel_path: str) -> str:
    with zipfile.ZipFile(zip_path) as zf:
        matches = [name for name in zf.namelist() if name.endswith(rel_path)]
    if not matches:
        raise FileNotFoundError(rel_path)
    if len(matches) > 1:
        matches = sorted(matches, key=len)
    return matches[0]


def table_exists(source: str | Path, rel_path: str) -> bool:
    source = Path(source)
    if source.is_dir():
        return (source / rel_path).exists()
    try:
        _resolve_member_name(source, rel_path)
        return True
    except Exception:
        return False


def read_csv(
    source: str | Path,
    rel_path: str,
    usecols: list[str] | None = None,
    parse_dates: list[str] | None = None,
    nrows: int | None = None,
) -> pd.DataFrame:
    source = Path(source)
    kwargs = dict(usecols=usecols, parse_dates=parse_dates, low_memory=False, nrows=nrows)
    if source.is_dir():
        return pd.read_csv(source / rel_path, compression="gzip", **kwargs)
    member = _resolve_member_name(source, rel_path)
    with zipfile.ZipFile(source) as zf, zf.open(member) as raw, gzip.GzipFile(fileobj=raw) as f:
        return pd.read_csv(f, **kwargs)


def get_columns(source: str | Path, rel_path: str) -> list[str]:
    return list(read_csv(source, rel_path, nrows=0).columns)


def read_csv_filtered(
    source: str | Path,
    rel_path: str,
    usecols: list[str] | None = None,
    parse_dates: list[str] | None = None,
    filter_func: Callable[[pd.DataFrame], pd.Series] | None = None,
    chunksize: int = 300_000,
) -> pd.DataFrame:
    source = Path(source)
    if source.is_dir():
        iterator = pd.read_csv(
            source / rel_path,
            usecols=usecols,
            parse_dates=parse_dates,
            compression="gzip",
            chunksize=chunksize,
            low_memory=False,
        )
        closer = None
    else:
        zf = zipfile.ZipFile(source)
        raw = zf.open(_resolve_member_name(source, rel_path))
        gz = gzip.GzipFile(fileobj=raw)
        iterator = pd.read_csv(
            gz,
            usecols=usecols,
            parse_dates=parse_dates,
            chunksize=chunksize,
            low_memory=False,
        )
        closer = (gz, raw, zf)

    pieces: list[pd.DataFrame] = []
    try:
        for chunk in iterator:
            sub = chunk if filter_func is None else chunk.loc[filter_func(chunk)]
            if len(sub):
                pieces.append(sub.copy())
    finally:
        if closer:
            for obj in closer:
                obj.close()

    if not pieces:
        return pd.DataFrame(columns=usecols)
    return pd.concat(pieces, ignore_index=True)


def master_admissions(source: str | Path) -> pd.DataFrame:
    adm = read_csv(
        source,
        "hosp/admissions.csv.gz",
        usecols=["subject_id", "hadm_id", "admittime", "dischtime", "deathtime", "hospital_expire_flag", "race"],
        parse_dates=["admittime", "dischtime", "deathtime"],
    )
    pat = read_csv(
        source,
        "hosp/patients.csv.gz",
        usecols=["subject_id", "gender", "anchor_age", "dod"],
        parse_dates=["dod"],
    )
    icu = read_csv(
        source,
        "icu/icustays.csv.gz",
        usecols=["subject_id", "hadm_id", "stay_id", "first_careunit", "last_careunit", "intime", "outtime"],
        parse_dates=["intime", "outtime"],
    )
    icu = icu.sort_values(["hadm_id", "intime"]).drop_duplicates("hadm_id", keep="first")
    out = adm.merge(pat, on="subject_id", how="left").merge(icu, on=["subject_id", "hadm_id"], how="left")
    out["dod"] = pd.to_datetime(out["dod"], errors="coerce")
    return out


def diagnoses(source: str | Path, hadm_ids: Iterable[int] | None = None) -> pd.DataFrame:
    usecols = ["subject_id", "hadm_id", "icd_code", "icd_version"]
    if hadm_ids is None:
        dx = read_csv(source, "hosp/diagnoses_icd.csv.gz", usecols=usecols)
    else:
        ids = set(int(x) for x in hadm_ids)
        dx = read_csv_filtered(
            source,
            "hosp/diagnoses_icd.csv.gz",
            usecols=usecols,
            filter_func=lambda c: c["hadm_id"].isin(ids),
        )
    dictionary = read_csv(
        source,
        "hosp/d_icd_diagnoses.csv.gz",
        usecols=["icd_code", "icd_version", "long_title"],
    )
    dx = dx.merge(dictionary, on=["icd_code", "icd_version"], how="left")
    dx["long_title"] = dx["long_title"].fillna("")
    dx["long_title_lower"] = dx["long_title"].str.lower()
    dx["icd_code_clean"] = dx["icd_code"].astype(str).str.replace(".", "", regex=False).str.upper()
    return dx


def procedures(source: str | Path, hadm_ids: Iterable[int] | None = None) -> pd.DataFrame:
    usecols = ["subject_id", "hadm_id", "chartdate", "icd_code", "icd_version"]
    parse_dates = ["chartdate"]
    if hadm_ids is None:
        px = read_csv(source, "hosp/procedures_icd.csv.gz", usecols=usecols, parse_dates=parse_dates)
    else:
        ids = set(int(x) for x in hadm_ids)
        px = read_csv_filtered(
            source,
            "hosp/procedures_icd.csv.gz",
            usecols=usecols,
            parse_dates=parse_dates,
            filter_func=lambda c: c["hadm_id"].isin(ids),
        )
    dictionary = read_csv(
        source,
        "hosp/d_icd_procedures.csv.gz",
        usecols=["icd_code", "icd_version", "long_title"],
    )
    px = px.merge(dictionary, on=["icd_code", "icd_version"], how="left")
    px["long_title"] = px["long_title"].fillna("")
    px["long_title_lower"] = px["long_title"].str.lower()
    return px


def prescriptions(source: str | Path, hadm_ids: Iterable[int] | None = None) -> pd.DataFrame:
    cols = ["subject_id", "hadm_id", "starttime", "stoptime", "drug", "route", "dose_val_rx", "dose_unit_rx"]
    if hadm_ids is None:
        rx = read_csv(source, "hosp/prescriptions.csv.gz", usecols=cols, parse_dates=["starttime", "stoptime"])
    else:
        ids = set(int(x) for x in hadm_ids)
        rx = read_csv_filtered(
            source,
            "hosp/prescriptions.csv.gz",
            usecols=cols,
            parse_dates=["starttime", "stoptime"],
            filter_func=lambda c: c["hadm_id"].isin(ids),
        )
    rx["drug"] = rx["drug"].fillna("")
    rx["drug_lower"] = rx["drug"].str.lower()
    rx["route"] = rx["route"].fillna("")
    rx["route_lower"] = rx["route"].str.lower().str.strip()
    rx["dose_val_num"] = pd.to_numeric(rx["dose_val_rx"], errors="coerce")
    return rx


def microbiology(source: str | Path, hadm_ids: Iterable[int]) -> pd.DataFrame:
    available = set(get_columns(source, "hosp/microbiologyevents.csv.gz"))
    wanted = [
        "subject_id", "hadm_id", "charttime", "chartdate", "storetime", "storedate",
        "org_name", "spec_type_desc", "test_name",
    ]
    cols = [c for c in wanted if c in available]
    required = {"subject_id", "hadm_id", "org_name", "spec_type_desc", "test_name"}
    missing = sorted(required - set(cols))
    if missing:
        raise ValueError(f"Missing required microbiologyevents columns: {missing}")
    ids = set(int(x) for x in hadm_ids)
    parse_dates = [c for c in ["charttime", "chartdate", "storetime", "storedate"] if c in cols]
    d = read_csv_filtered(
        source,
        "hosp/microbiologyevents.csv.gz",
        usecols=cols,
        parse_dates=parse_dates,
        filter_func=lambda c: c["hadm_id"].isin(ids),
    )
    for c in ["charttime", "chartdate", "storetime", "storedate"]:
        if c not in d:
            d[c] = pd.NaT
        d[c] = pd.to_datetime(d[c], errors="coerce")
    for c in ["org_name", "spec_type_desc", "test_name"]:
        d[c] = d[c].fillna("")
        d[f"{c}_lower"] = d[c].str.lower()
    d["specimen_time"] = d["charttime"].fillna(d["chartdate"])
    d["result_available_time"] = d["storetime"].fillna(d["storedate"])
    return d


def inputevents(source: str | Path, stay_ids: Iterable[int]) -> pd.DataFrame:
    ids = set(int(x) for x in stay_ids)
    cols = [
        "subject_id", "hadm_id", "stay_id", "starttime", "endtime", "itemid",
        "amount", "amountuom", "ordercategoryname", "statusdescription",
    ]
    d_items = read_csv(source, "icu/d_items.csv.gz", usecols=["itemid", "label"])
    d = read_csv_filtered(
        source,
        "icu/inputevents.csv.gz",
        usecols=cols,
        parse_dates=["starttime", "endtime"],
        filter_func=lambda c: c["stay_id"].isin(ids),
    )
    d = d.merge(d_items, on="itemid", how="left")
    d["label"] = d["label"].fillna("")
    d["label_lower"] = d["label"].str.lower()
    d["amount_num"] = pd.to_numeric(d["amount"], errors="coerce")
    return d

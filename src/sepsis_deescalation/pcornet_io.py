from __future__ import annotations

from pathlib import Path
from typing import Any

import duckdb
import pandas as pd


def read_table(path: str | Path, columns: list[str] | None = None) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path, columns=columns)
    if suffix in {".csv", ".gz"}:
        return pd.read_csv(path, usecols=columns, low_memory=False)
    raise ValueError(f"Unsupported PCORnet source format: {path}")


def columns_for(path: str | Path) -> list[str]:
    path = Path(path)
    if path.suffix.lower() == ".parquet":
        con = duckdb.connect()
        try:
            rows = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{path.as_posix()}')").fetchall()
            return [r[0] for r in rows]
        finally:
            con.close()
    return list(pd.read_csv(path, nrows=0).columns)


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d.columns = [str(c).strip().lower() for c in d.columns]
    return d


def combine_date_time(date: pd.Series, time: pd.Series | None = None) -> pd.Series:
    """Combine PCORnet date and time columns into timestamps.

    PCORnet time fields may arrive as numeric HHMM values or text such as HH:MM / HH:MM:SS.
    Missing or unparseable times fall back to midnight on the supplied date.
    """
    base = pd.to_datetime(date, errors="coerce").dt.normalize()
    if time is None:
        return base

    if pd.api.types.is_numeric_dtype(time):
        vals = pd.to_numeric(time, errors="coerce")
        hours = (vals // 100).where(vals.notna())
        minutes = (vals % 100).where(vals.notna())
        valid = hours.between(0, 23) & minutes.between(0, 59)
        offset = pd.to_timedelta(hours.where(valid, 0).fillna(0), unit="h") + pd.to_timedelta(
            minutes.where(valid, 0).fillna(0), unit="m"
        )
        return base + offset

    text = time.astype("string").str.strip()
    hhmm = text.str.extract(r"^(?P<hour>\d{1,2}):(?P<minute>\d{2})(?::(?P<second>\d{2}(?:\.\d+)?))?$")
    hours = pd.to_numeric(hhmm["hour"], errors="coerce")
    minutes = pd.to_numeric(hhmm["minute"], errors="coerce")
    seconds = pd.to_numeric(hhmm["second"], errors="coerce").fillna(0)
    valid = hours.between(0, 23) & minutes.between(0, 59) & seconds.between(0, 59.999999)
    offset = (
        pd.to_timedelta(hours.where(valid, 0).fillna(0), unit="h")
        + pd.to_timedelta(minutes.where(valid, 0).fillna(0), unit="m")
        + pd.to_timedelta(seconds.where(valid, 0).fillna(0), unit="s")
    )
    return base + offset


def resolve_path(cfg: dict[str, Any], key: str) -> Path:
    raw = cfg.get("paths", {}).get(key)
    if not raw:
        raise KeyError(f"Missing paths.{key} in PCORnet configuration")
    return Path(raw).expanduser()

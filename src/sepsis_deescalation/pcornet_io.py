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
    base = pd.to_datetime(date, errors="coerce")
    if time is None:
        return base
    t = time.copy()
    if pd.api.types.is_numeric_dtype(t):
        vals = pd.to_numeric(t, errors="coerce")
        hours = (vals // 100).fillna(0).astype(int)
        minutes = (vals % 100).fillna(0).astype(int)
        return base + pd.to_timedelta(hours, unit="h") + pd.to_timedelta(minutes, unit="m")
    text = t.astype(str).str.strip()
    parsed = pd.to_timedelta(text.where(text.str.match(r"^\d{1,2}:\d{2}$"), other=""), errors="coerce")
    return base + parsed.fillna(pd.Timedelta(0))


def resolve_path(cfg: dict[str, Any], key: str) -> Path:
    raw = cfg.get("paths", {}).get(key)
    if not raw:
        raise KeyError(f"Missing paths.{key} in PCORnet configuration")
    return Path(raw).expanduser()

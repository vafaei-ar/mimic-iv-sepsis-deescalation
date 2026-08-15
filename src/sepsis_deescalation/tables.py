from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd

from .stats import smd_binary, smd_continuous, weighted_mean


def table1(
    df: pd.DataFrame,
    continuous_vars: Iterable[str],
    binary_vars: Iterable[str],
    treatment_col: str = "A",
    weight_col: str | None = None,
) -> pd.DataFrame:
    rows: list[dict] = []
    t = df[treatment_col] == 1
    c = df[treatment_col] == 0
    wt = df.loc[t, weight_col] if weight_col else None
    wc = df.loc[c, weight_col] if weight_col else None

    for var in continuous_vars:
        if var not in df:
            continue
        x = pd.to_numeric(df[var], errors="coerce")
        fill = x.median() if x.notna().any() else 0.0
        x = x.fillna(fill)
        mt = float(x[t].mean()) if weight_col is None else weighted_mean(x[t], wt)
        mc = float(x[c].mean()) if weight_col is None else weighted_mean(x[c], wc)
        smd = smd_continuous(x[t], x[c], wt, wc)
        rows.append({"variable": var, "type": "continuous", "treated": mt, "continued": mc, "SMD": smd, "abs_SMD": abs(smd)})

    for var in binary_vars:
        if var not in df:
            continue
        x = pd.to_numeric(df[var], errors="coerce").fillna(0)
        pt = float(x[t].mean()) if weight_col is None else weighted_mean(x[t], wt)
        pc = float(x[c].mean()) if weight_col is None else weighted_mean(x[c], wc)
        smd = smd_binary(x[t], x[c], wt, wc)
        rows.append({"variable": var, "type": "binary", "treated": pt, "continued": pc, "SMD": smd, "abs_SMD": abs(smd)})
    return pd.DataFrame(rows)


def cohort_flow_row(step: str, n: int) -> dict:
    return {"step": step, "n": int(n)}

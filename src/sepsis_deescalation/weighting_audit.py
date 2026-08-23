from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .stats import balance_table, effective_sample_size, risks, weighted_mean


def add_overlap_weights(
    df: pd.DataFrame,
    ps_col: str = "ps_den",
    treatment_col: str = "A",
    out_col: str = "OW_A",
) -> pd.DataFrame:
    """Add overlap weights using an already fitted propensity score.

    These weights target the overlap population, not the ATE targeted by the
    primary stabilized IPTW analysis. They are therefore a sensitivity analysis.
    """
    d = df.copy()
    ps = pd.to_numeric(d[ps_col], errors="coerce").clip(0.001, 0.999)
    a = pd.to_numeric(d[treatment_col], errors="coerce")
    d[out_col] = np.where(a == 1, 1.0 - ps, ps)
    return d


def weight_tail_diagnostics(
    df: pd.DataFrame,
    weight_col: str = "SW_A",
    ps_col: str = "ps_den",
    treatment_col: str = "A",
) -> pd.DataFrame:
    rows: list[dict] = []
    for a, label in [(1, "deescalated_stopped"), (0, "continued_broad")]:
        sub = df.loc[df[treatment_col] == a].copy()
        w = pd.to_numeric(sub[weight_col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        ps = pd.to_numeric(sub[ps_col], errors="coerce")
        ordered = w.sort_values(ascending=False)
        total = float(ordered.sum())
        row = {
            "group": label,
            "n": len(sub),
            "effective_sample_size": effective_sample_size(w),
            "max_weight": float(w.max()) if len(w) else np.nan,
            "p99_weight": float(w.quantile(0.99)) if len(w) else np.nan,
            "ps_min": float(ps.min()) if ps.notna().any() else np.nan,
            "ps_p01": float(ps.quantile(0.01)) if ps.notna().any() else np.nan,
            "ps_p05": float(ps.quantile(0.05)) if ps.notna().any() else np.nan,
            "n_ps_lt_0_01": int((ps < 0.01).sum()),
            "n_ps_lt_0_025": int((ps < 0.025).sum()),
            "n_ps_lt_0_05": int((ps < 0.05).sum()),
        }
        for k in [1, 5, 10, 20, 50, 100]:
            row[f"top_{k}_weight_share"] = float(ordered.head(k).sum() / total) if total > 0 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def exact_duplicate_covariates(df: pd.DataFrame, variables: Iterable[str]) -> pd.DataFrame:
    """Find exact duplicate numeric covariates after the same simple missing handling used for auditing."""
    vars_present = [v for v in variables if v in df.columns]
    prepared: dict[str, pd.Series] = {}
    for var in vars_present:
        x = pd.to_numeric(df[var], errors="coerce")
        fill = 0.0 if set(x.dropna().unique()).issubset({0, 1, 0.0, 1.0}) else float(x.mean())
        prepared[var] = x.fillna(fill)
    rows = []
    for i, left in enumerate(vars_present):
        for right in vars_present[i + 1 :]:
            if prepared[left].equals(prepared[right]):
                rows.append({"variable_1": left, "variable_2": right, "relationship": "exact_duplicate"})
    return pd.DataFrame(rows, columns=["variable_1", "variable_2", "relationship"])


def design_matrix_summary(df: pd.DataFrame, ps_variables: Iterable[str]) -> pd.DataFrame:
    vars_present = [v for v in ps_variables if v in df.columns]
    if not vars_present:
        return pd.DataFrame([{"n_rows": len(df), "n_ps_columns": 0, "rank_with_intercept": 0, "rank_deficiency": 0}])
    x = df[vars_present].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float)
    x = np.column_stack([np.ones(len(x)), x])
    rank = int(np.linalg.matrix_rank(x))
    ncols = int(x.shape[1])
    return pd.DataFrame(
        [{
            "n_rows": int(x.shape[0]),
            "n_ps_columns": len(vars_present),
            "n_columns_with_intercept": ncols,
            "rank_with_intercept": rank,
            "rank_deficiency": ncols - rank,
        }]
    )


def _weighted_summary(
    df: pd.DataFrame,
    variables: Iterable[str],
    weight_col: str,
    label: str,
    estimand: str,
) -> dict:
    rt, rc, rd, rr = risks(df, "death_by_horizon", weight_col)
    bal = balance_table(df, variables, weight_col=weight_col)
    t = df["A"] == 1
    c = df["A"] == 0
    return {
        "analysis": label,
        "estimand": estimand,
        "risk_deescalated_stopped": rt,
        "risk_continued": rc,
        "risk_difference": rd,
        "risk_ratio": rr,
        "max_post_smd": float(bal["after"].max()) if len(bal) else np.nan,
        "worst_balanced_variable": str(bal.sort_values("after", ascending=False).iloc[0]["variable"]) if len(bal) else "",
        "ess_deescalated_stopped": effective_sample_size(df.loc[t, weight_col]),
        "ess_continued": effective_sample_size(df.loc[c, weight_col]),
        "max_weight": float(pd.to_numeric(df[weight_col], errors="coerce").max()),
    }


def run_weighting_audit(
    cohort_w: pd.DataFrame,
    candidate_vars: Iterable[str],
    out_dir: str | Path | None = None,
) -> dict[str, pd.DataFrame]:
    """Audit positivity and weight sensitivity without changing the primary estimand.

    The audit compares the original ATE IPTW result with deterministic weight
    truncation and an overlap-weight sensitivity. It does not select a preferred
    result based on effect size.
    """
    d = cohort_w.copy()
    rows = [_weighted_summary(d, candidate_vars, "SW_A", "Primary stabilized IPTW", "ATE")]

    for low, high in [(1.0, 99.0), (2.5, 97.5)]:
        lo = float(d["SW_A"].quantile(low / 100.0))
        hi = float(d["SW_A"].quantile(high / 100.0))
        col = f"SW_A_trunc_{low:g}_{high:g}"
        tmp = d.copy()
        tmp[col] = tmp["SW_A"].clip(lo, hi)
        rows.append(_weighted_summary(tmp, candidate_vars, col, f"IPTW truncated {low:g}/{high:g}", "ATE, truncated weights"))

    overlap = add_overlap_weights(d)
    rows.append(_weighted_summary(overlap, candidate_vars, "OW_A", "Overlap weighting", "overlap population"))

    summary = pd.DataFrame(rows)
    balance = balance_table(d, candidate_vars).sort_values("after", ascending=False).reset_index(drop=True)
    tails = weight_tail_diagnostics(d)
    duplicates = exact_duplicate_covariates(d, candidate_vars)

    if out_dir is not None:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out / "weighting_audit_summary.csv", index=False)
        balance.head(25).to_csv(out / "top_imbalanced_covariates.csv", index=False)
        tails.to_csv(out / "weight_tail_diagnostics.csv", index=False)
        duplicates.to_csv(out / "exact_duplicate_covariates.csv", index=False)

    return {
        "summary": summary,
        "top_imbalanced": balance.head(25),
        "tails": tails,
        "duplicates": duplicates,
    }

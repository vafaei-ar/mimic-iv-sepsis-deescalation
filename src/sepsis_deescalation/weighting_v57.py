from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import statsmodels.api as sm

from .stats import (
    balance_table,
    effective_sample_size,
    fit_stabilized_iptw,
    risks,
)


def deduplicate_covariates(df: pd.DataFrame, variables: Iterable[str]) -> tuple[list[str], pd.DataFrame]:
    """Remove exact duplicate analysis covariates while preserving the first occurrence."""
    kept: list[str] = []
    prepared: dict[str, pd.Series] = {}
    removed: list[dict] = []
    for var in variables:
        if var not in df.columns:
            continue
        x = pd.to_numeric(df[var], errors="coerce")
        if x.isna().all() or x.nunique(dropna=True) < 2:
            continue
        vals = set(x.dropna().unique().tolist())
        fill = 0.0 if vals.issubset({0, 1, 0.0, 1.0}) else float(x.mean())
        x2 = x.fillna(fill)
        duplicate_of = next((old for old in kept if prepared[old].equals(x2)), None)
        if duplicate_of is None:
            kept.append(var)
            prepared[var] = x2
        else:
            removed.append({"removed_variable": var, "retained_variable": duplicate_of, "reason": "exact_duplicate"})
    return kept, pd.DataFrame(removed, columns=["removed_variable", "retained_variable", "reason"])


def add_overlap_weights(df: pd.DataFrame, ps_col: str = "ps_den", treatment_col: str = "A") -> pd.DataFrame:
    d = df.copy()
    ps = pd.to_numeric(d[ps_col], errors="coerce").clip(0.001, 0.999)
    a = pd.to_numeric(d[treatment_col], errors="coerce")
    d["OW_A"] = np.where(a == 1, 1.0 - ps, ps)
    return d


def _summary(df: pd.DataFrame, vars_: Sequence[str], weight_col: str, label: str, estimand: str) -> dict:
    rt, rc, rd, rr = risks(df, "death_by_horizon", weight_col)
    bal = balance_table(df, vars_, weight_col=weight_col)
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


def fit_quadratic_ps(
    df: pd.DataFrame,
    base_vars: Sequence[str],
    nonlinear_vars: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """PS sensitivity with prespecified quadratic terms for selected continuous covariates."""
    base_w, fit, diag = fit_stabilized_iptw(df, base_vars)
    d = base_w.copy()
    exog = pd.DataFrame(np.asarray(fit.model.exog, dtype=float), index=d.index)
    names = list(getattr(fit.model, "exog_names", [f"x{i}" for i in range(exog.shape[1])]))
    exog.columns = names
    added = []
    for var in nonlinear_vars:
        z = f"{var}_z"
        if z in d.columns:
            col = f"{z}_sq"
            exog[col] = np.square(pd.to_numeric(d[z], errors="coerce").fillna(0.0))
            added.append({"variable": var, "term": col, "type": "quadratic"})
    model = sm.GLM(pd.to_numeric(d["A"], errors="coerce"), exog, family=sm.families.Binomial())
    try:
        qfit = model.fit(maxiter=200, disp=0)
        if not bool(getattr(qfit, "converged", True)):
            qfit = model.fit_regularized(alpha=0.001, L1_wt=0.0, maxiter=200)
    except Exception:
        qfit = model.fit_regularized(alpha=0.001, L1_wt=0.0, maxiter=200)
    params = np.asarray(qfit.params, dtype=float)
    pden = 1 / (1 + np.exp(-np.clip(exog.to_numpy(dtype=float) @ params, -35, 35)))
    pden = np.clip(pden, 0.001, 0.999)
    pnum = float(d["A"].mean())
    d["ps_den_quadratic"] = pden
    d["SW_A_quadratic"] = np.where(d["A"] == 1, pnum / pden, (1 - pnum) / (1 - pden))
    return d, pd.DataFrame(added)


def _bootstrap_strategy(
    df: pd.DataFrame,
    variables: Sequence[str],
    strategy: str,
    reps: int,
    seed: int,
    nonlinear_vars: Sequence[str] | None = None,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    n = len(df)
    for rep in range(reps):
        s = df.iloc[rng.integers(0, n, size=n)].copy()
        if s["A"].nunique() != 2:
            continue
        try:
            w, _, _ = fit_stabilized_iptw(s, variables)
            if strategy == "overlap":
                w = add_overlap_weights(w)
                weight_col = "OW_A"
            elif strategy.startswith("truncate_"):
                _, low, high = strategy.split("_")
                low_f, high_f = float(low), float(high)
                lo = w["SW_A"].quantile(low_f / 100.0)
                hi = w["SW_A"].quantile(high_f / 100.0)
                weight_col = "SW_A_boot_trunc"
                w[weight_col] = w["SW_A"].clip(lo, hi)
            elif strategy == "quadratic":
                w, _ = fit_quadratic_ps(s, variables, nonlinear_vars or [])
                weight_col = "SW_A_quadratic"
            else:
                raise ValueError(strategy)
            rt, rc, rd, rr = risks(w, "death_by_horizon", weight_col)
            rows.append({"rep": rep, "risk_treated": rt, "risk_control": rc, "risk_difference": rd, "risk_ratio": rr})
        except Exception:
            continue
    return pd.DataFrame(rows)


def _ci(boot: pd.DataFrame, label: str, requested: int) -> pd.DataFrame:
    rows = []
    for est in ["risk_difference", "risk_ratio"]:
        vals = pd.to_numeric(boot.get(est), errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        rows.append({"analysis": label, "estimand": est, "lower_95": vals.quantile(0.025), "upper_95": vals.quantile(0.975), "n_success": len(vals), "n_requested": requested})
    return pd.DataFrame(rows)


def run_v57_weighting_sensitivities(
    cohort: pd.DataFrame,
    candidate_vars: Sequence[str],
    out_dir: str | Path,
    reps: int = 1000,
    seed: int = 20260426,
    nonlinear_vars: Sequence[str] = ("temperature_48_72h", "broad_abx_hours_pre72", "systemic_abx_hours_pre72"),
) -> dict[str, pd.DataFrame]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    dedup_vars, dedup = deduplicate_covariates(cohort, candidate_vars)
    dedup.to_csv(out / "ps_exact_duplicates_removed.csv", index=False)

    base, _, _ = fit_stabilized_iptw(cohort, dedup_vars)
    rows = [_summary(base, dedup_vars, "SW_A", "Deduplicated primary PS", "ATE")]

    overlap = add_overlap_weights(base)
    rows.append(_summary(overlap, dedup_vars, "OW_A", "Overlap weighting", "overlap population"))

    for low, high in [(1.0, 99.0), (2.5, 97.5)]:
        tmp = base.copy()
        col = f"SW_A_trunc_{low:g}_{high:g}"
        tmp[col] = tmp["SW_A"].clip(tmp["SW_A"].quantile(low / 100.0), tmp["SW_A"].quantile(high / 100.0))
        rows.append(_summary(tmp, dedup_vars, col, f"IPTW truncated {low:g}/{high:g}", "ATE, truncated weights"))

    q, qterms = fit_quadratic_ps(cohort, dedup_vars, nonlinear_vars)
    rows.append(_summary(q, dedup_vars, "SW_A_quadratic", "Quadratic PS sensitivity", "ATE"))
    qterms.to_csv(out / "quadratic_ps_terms.csv", index=False)

    summary = pd.DataFrame(rows)
    summary.to_csv(out / "v57_weighting_point_estimates.csv", index=False)

    specs = [
        ("overlap", "Overlap weighting", seed + 1),
        ("truncate_1_99", "IPTW truncated 1/99", seed + 2),
        ("truncate_2.5_97.5", "IPTW truncated 2.5/97.5", seed + 3),
        ("quadratic", "Quadratic PS sensitivity", seed + 4),
    ]
    ci_parts = []
    for strategy, label, sseed in specs:
        boot = _bootstrap_strategy(cohort, dedup_vars, strategy, reps, sseed, nonlinear_vars)
        boot.to_csv(out / f"bootstrap_{strategy.replace('.', '_')}.csv", index=False)
        ci_parts.append(_ci(boot, label, reps))
    ci = pd.concat(ci_parts, ignore_index=True)
    ci.to_csv(out / "v57_weighting_bootstrap_ci.csv", index=False)
    return {"summary": summary, "ci": ci, "duplicates": dedup, "quadratic_terms": qterms}

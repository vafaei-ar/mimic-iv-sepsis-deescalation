from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


def weighted_mean(x, w) -> float:
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    mask = np.isfinite(x) & np.isfinite(w)
    if not mask.any() or w[mask].sum() <= 0:
        return np.nan
    return float(np.sum(x[mask] * w[mask]) / np.sum(w[mask]))


def weighted_var(x, w) -> float:
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    mask = np.isfinite(x) & np.isfinite(w)
    if not mask.any() or w[mask].sum() <= 0:
        return np.nan
    m = weighted_mean(x[mask], w[mask])
    return float(np.sum(w[mask] * (x[mask] - m) ** 2) / np.sum(w[mask]))


def smd_continuous(x_t, x_c, w_t=None, w_c=None) -> float:
    x_t = np.asarray(x_t, dtype=float)
    x_c = np.asarray(x_c, dtype=float)
    mt = np.nanmean(x_t) if w_t is None else weighted_mean(x_t, w_t)
    mc = np.nanmean(x_c) if w_c is None else weighted_mean(x_c, w_c)
    vt = np.nanvar(x_t) if w_t is None else weighted_var(x_t, w_t)
    vc = np.nanvar(x_c) if w_c is None else weighted_var(x_c, w_c)
    sd = np.sqrt((vt + vc) / 2.0)
    return 0.0 if not np.isfinite(sd) or sd == 0 else float((mt - mc) / sd)


def smd_binary(x_t, x_c, w_t=None, w_c=None) -> float:
    pt = np.nanmean(x_t) if w_t is None else weighted_mean(x_t, w_t)
    pc = np.nanmean(x_c) if w_c is None else weighted_mean(x_c, w_c)
    den = np.sqrt((pt * (1 - pt) + pc * (1 - pc)) / 2.0)
    return 0.0 if not np.isfinite(den) or den == 0 else float((pt - pc) / den)


def _is_binary_01(x: pd.Series) -> bool:
    values = set(pd.Series(x).dropna().unique().tolist())
    return values.issubset({0, 1, 0.0, 1.0})


def prepare_ps_covariates(df: pd.DataFrame, candidate_vars: Iterable[str]) -> tuple[pd.DataFrame, list[str], pd.DataFrame]:
    d = df.copy()
    used: list[str] = []
    rows: list[dict] = []
    for var in candidate_vars:
        if var not in d.columns:
            continue
        x = pd.to_numeric(d[var], errors="coerce")
        if x.notna().sum() == 0:
            continue
        missing = x.isna().astype(int)
        if missing.any() and not var.endswith("_missing"):
            miss_name = f"{var}_ps_missing"
            d[miss_name] = missing
            if d[miss_name].nunique() > 1:
                used.append(miss_name)
        if _is_binary_01(x.dropna()) and x.dropna().nunique() <= 2:
            fill = float(x.median()) if x.notna().any() else 0.0
            d[var] = x.fillna(fill).astype(float)
            if d[var].nunique() > 1:
                used.append(var)
            rows.append({"variable": var, "mean": np.nan, "sd": np.nan, "type": "binary", "fill": fill})
            continue
        fill = float(x.median()) if x.notna().any() else 0.0
        x = x.fillna(fill)
        mean = float(x.mean())
        sd = float(x.std(ddof=0))
        if not np.isfinite(sd) or sd <= 0:
            continue
        ps_name = f"{var}_ps_z"
        d[ps_name] = (x - mean) / sd
        used.append(ps_name)
        rows.append({"variable": var, "ps_variable": ps_name, "mean": mean, "sd": sd, "type": "continuous", "fill": fill})
    return d, used, pd.DataFrame(rows)


def fit_binomial_formula(formula: str, data: pd.DataFrame):
    try:
        return smf.glm(formula=formula, data=data, family=sm.families.Binomial()).fit()
    except Exception as exc:
        warnings.warn(f"Standard GLM failed ({exc}); using regularized binomial fit.")
        model = smf.glm(formula=formula, data=data, family=sm.families.Binomial())
        return model.fit_regularized(alpha=1e-6, L1_wt=0.0, maxiter=1000)


def fit_stabilized_iptw(
    cohort: pd.DataFrame,
    candidate_vars: Iterable[str],
    treatment_col: str = "A",
) -> tuple[pd.DataFrame, object, dict]:
    d, used, std_table = prepare_ps_covariates(cohort, candidate_vars)
    if d[treatment_col].nunique() != 2:
        raise ValueError("Treatment must have two levels.")
    formula = treatment_col + " ~ " + (" + ".join(used) if used else "1")
    fit = fit_binomial_formula(formula, d)
    ps = np.asarray(fit.predict(d), dtype=float)
    ps = np.clip(ps, 0.01, 0.99)
    p_treated = float(d[treatment_col].mean())
    d["ps_den"] = ps
    d["SW_A"] = np.where(d[treatment_col] == 1, p_treated / ps, (1 - p_treated) / (1 - ps))
    diagnostics = {
        "formula": formula,
        "used_vars": used,
        "standardization_table": std_table,
        "treated_probability": p_treated,
    }
    return d, fit, diagnostics


def effective_sample_size(weights) -> float:
    w = np.asarray(weights, dtype=float)
    w = w[np.isfinite(w)]
    if len(w) == 0 or np.sum(w * w) == 0:
        return np.nan
    return float((w.sum() ** 2) / np.sum(w * w))


def weight_summary(df: pd.DataFrame, weight_col: str = "SW_A", treatment_col: str = "A") -> pd.DataFrame:
    rows = []
    for label, sub in [("overall", df), ("deescalated_stopped", df.loc[df[treatment_col] == 1]), ("continued", df.loc[df[treatment_col] == 0])]:
        w = pd.to_numeric(sub[weight_col], errors="coerce").dropna()
        rows.append(
            {
                "group": label,
                "n": len(sub),
                "ess": effective_sample_size(w),
                "mean": w.mean(),
                "median": w.median(),
                "min": w.min(),
                "p01": w.quantile(0.01),
                "p05": w.quantile(0.05),
                "p95": w.quantile(0.95),
                "p99": w.quantile(0.99),
                "max": w.max(),
            }
        )
    return pd.DataFrame(rows)


def risks(df: pd.DataFrame, outcome: str, weight_col: str | None = None, treatment_col: str = "A") -> tuple[float, float, float, float]:
    t = df.loc[df[treatment_col] == 1]
    c = df.loc[df[treatment_col] == 0]
    rt = float(t[outcome].mean()) if weight_col is None else weighted_mean(t[outcome], t[weight_col])
    rc = float(c[outcome].mean()) if weight_col is None else weighted_mean(c[outcome], c[weight_col])
    rd = rt - rc
    rr = rt / rc if rc > 0 else np.nan
    return rt, rc, rd, rr


def mean_difference(df: pd.DataFrame, outcome: str, weight_col: str | None = None, treatment_col: str = "A") -> tuple[float, float, float]:
    t = df.loc[df[treatment_col] == 1]
    c = df.loc[df[treatment_col] == 0]
    mt = float(t[outcome].mean()) if weight_col is None else weighted_mean(t[outcome], t[weight_col])
    mc = float(c[outcome].mean()) if weight_col is None else weighted_mean(c[outcome], c[weight_col])
    return mt, mc, mt - mc


def balance_table(
    df: pd.DataFrame,
    variables: Iterable[str],
    treatment_col: str = "A",
    weight_col: str = "SW_A",
) -> pd.DataFrame:
    rows = []
    for var in variables:
        if var not in df:
            continue
        x = pd.to_numeric(df[var], errors="coerce")
        if x.notna().sum() == 0:
            continue
        fill = float(x.median()) if x.notna().any() else 0.0
        x = x.fillna(fill)
        tmask = df[treatment_col] == 1
        cmask = ~tmask
        binary = _is_binary_01(x)
        fn = smd_binary if binary else smd_continuous
        before = fn(x[tmask], x[cmask])
        after = fn(x[tmask], x[cmask], df.loc[tmask, weight_col], df.loc[cmask, weight_col])
        rows.append({"variable": var, "before": abs(before), "after": abs(after), "signed_before": before, "signed_after": after})
    return pd.DataFrame(rows).sort_values("before", ascending=False).reset_index(drop=True)


def bootstrap_iptw_ci(
    df: pd.DataFrame,
    candidate_vars: Iterable[str],
    outcome: str,
    kind: str,
    reps: int,
    seed: int,
    treatment_col: str = "A",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    records: list[dict] = []
    n = len(df)
    for b in range(reps):
        sample = df.iloc[rng.integers(0, n, size=n)].copy()
        if sample[treatment_col].nunique() < 2:
            continue
        try:
            w, _, _ = fit_stabilized_iptw(sample, candidate_vars, treatment_col=treatment_col)
            if kind == "risk":
                rt, rc, rd, rr = risks(w, outcome, "SW_A", treatment_col)
                records.append({"rep": b, "risk_treated": rt, "risk_control": rc, "risk_difference": rd, "risk_ratio": rr})
            elif kind == "mean":
                mt, mc, md = mean_difference(w, outcome, "SW_A", treatment_col)
                records.append({"rep": b, "mean_treated": mt, "mean_control": mc, "mean_difference": md})
            else:
                raise ValueError(f"Unknown bootstrap kind: {kind}")
        except Exception:
            continue
    boot = pd.DataFrame(records)
    if boot.empty:
        return pd.DataFrame(), boot
    estimands = ["risk_difference", "risk_ratio"] if kind == "risk" else ["mean_difference"]
    ci_rows = []
    for est in estimands:
        vals = pd.to_numeric(boot[est], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        ci_rows.append(
            {
                "estimand": est,
                "lower_95": vals.quantile(0.025),
                "upper_95": vals.quantile(0.975),
                "n_success": len(vals),
                "n_requested": reps,
            }
        )
    return pd.DataFrame(ci_rows), boot


def truncate_weights(df: pd.DataFrame, low: float, high: float, weight_col: str = "SW_A") -> pd.DataFrame:
    d = df.copy()
    lo = d[weight_col].quantile(low / 100.0)
    hi = d[weight_col].quantile(high / 100.0)
    new_col = f"{weight_col}_trunc_{low:g}_{high:g}"
    d[new_col] = d[weight_col].clip(lo, hi)
    return d

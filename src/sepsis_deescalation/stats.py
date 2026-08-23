from __future__ import annotations

import warnings
from typing import Iterable

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.special import expit


def weighted_mean(x, w) -> float:
    x = np.asarray(x, dtype=float); w = np.asarray(w, dtype=float)
    mask = np.isfinite(x) & np.isfinite(w)
    if not mask.any() or w[mask].sum() <= 0: return np.nan
    return float(np.sum(x[mask] * w[mask]) / np.sum(w[mask]))


def weighted_var(x, w) -> float:
    x = np.asarray(x, dtype=float); w = np.asarray(w, dtype=float)
    mask = np.isfinite(x) & np.isfinite(w)
    if not mask.any() or w[mask].sum() <= 0: return np.nan
    m = weighted_mean(x[mask], w[mask]); return float(np.sum(w[mask] * (x[mask] - m) ** 2) / np.sum(w[mask]))


def smd_continuous(x_t, x_c, w_t=None, w_c=None) -> float:
    x_t = np.asarray(x_t, dtype=float); x_c = np.asarray(x_c, dtype=float)
    mt = np.nanmean(x_t) if w_t is None else weighted_mean(x_t, w_t); mc = np.nanmean(x_c) if w_c is None else weighted_mean(x_c, w_c)
    vt = np.nanvar(x_t) if w_t is None else weighted_var(x_t, w_t); vc = np.nanvar(x_c) if w_c is None else weighted_var(x_c, w_c)
    sd = np.sqrt((vt + vc) / 2.0); return 0.0 if not np.isfinite(sd) or sd == 0 else float((mt - mc) / sd)


def smd_binary(x_t, x_c, w_t=None, w_c=None) -> float:
    pt = np.nanmean(x_t) if w_t is None else weighted_mean(x_t, w_t); pc = np.nanmean(x_c) if w_c is None else weighted_mean(x_c, w_c)
    den = np.sqrt((pt * (1 - pt) + pc * (1 - pc)) / 2.0); return 0.0 if not np.isfinite(den) or den == 0 else float((pt - pc) / den)


def _is_binary_01(x: pd.Series) -> bool:
    vals = pd.Series(x).dropna().unique()
    return len(vals) > 0 and set(vals).issubset({0, 1, 0.0, 1.0})


def prepare_ps_covariates(df: pd.DataFrame, candidate_vars: Iterable[str]) -> tuple[pd.DataFrame, list[str], pd.DataFrame]:
    """Exact scripted equivalent of the v5.5 stable PS preparation."""
    d = df.copy(); ps_vars = []; rows = []
    for var in candidate_vars:
        if var not in d.columns: continue
        x = pd.to_numeric(d[var], errors="coerce")
        if x.isna().all() or x.nunique(dropna=True) < 2: continue
        if _is_binary_01(x):
            d[var] = x.fillna(0).astype(float); ps_vars.append(var)
            rows.append({"original_variable": var, "ps_variable": var, "type": "binary", "mean": np.nan, "sd": np.nan})
        else:
            mu = float(x.mean()); sd = float(x.std(ddof=0))
            if not np.isfinite(sd) or sd == 0: continue
            zname = var + "_z"; d[zname] = ((x.fillna(mu) - mu) / sd).clip(lower=-8, upper=8); ps_vars.append(zname)
            rows.append({"original_variable": var, "ps_variable": zname, "type": "continuous_standardized", "mean": mu, "sd": sd})
    return d, ps_vars, pd.DataFrame(rows)


def _build_formula(df: pd.DataFrame, treatment_col: str, variables: Iterable[str]) -> tuple[str, list[str]]:
    keep = [v for v in variables if v in df.columns and not df[v].isna().all() and df[v].nunique(dropna=True) >= 2]
    return treatment_col + " ~ " + (" + ".join(keep) if keep else "1"), keep


def fit_binomial_formula(formula: str, data: pd.DataFrame):
    model = smf.glm(formula=formula, data=data, family=sm.families.Binomial())
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="overflow encountered in exp", category=RuntimeWarning)
        try:
            fit = model.fit(maxiter=200, disp=0)
            if not bool(getattr(fit, "converged", True)):
                fit = model.fit_regularized(alpha=0.001, L1_wt=0.0, maxiter=200); return fit, "regularized_glm"
            return fit, "glm"
        except Exception as exc:
            warnings.warn(f"Standard GLM failed ({exc}); using regularized binomial fit.")
            return model.fit_regularized(alpha=0.001, L1_wt=0.0, maxiter=200), "regularized_glm"


def _safe_predict_training_data(fit) -> np.ndarray:
    exog = np.asarray(fit.model.exog, dtype=float); params = np.asarray(fit.params, dtype=float)
    linpred = np.clip(exog @ params, -35, 35); return expit(linpred)


def fit_stabilized_iptw(cohort: pd.DataFrame, candidate_vars: Iterable[str], treatment_col: str = "A") -> tuple[pd.DataFrame, object, dict]:
    d, ps_vars, std_table = prepare_ps_covariates(cohort, candidate_vars)
    if d[treatment_col].nunique() != 2: raise ValueError("Treatment must have two levels.")
    formula, used = _build_formula(d, treatment_col, ps_vars)
    num_fit, num_method = fit_binomial_formula(f"{treatment_col} ~ 1", d); den_fit, den_method = fit_binomial_formula(formula, d)
    pnum = np.clip(_safe_predict_training_data(num_fit), 0.001, 0.999); pden = np.clip(_safe_predict_training_data(den_fit), 0.001, 0.999)
    d["ps_den"] = pden; d["ps_num"] = pnum; d["SW_A"] = np.where(d[treatment_col] == 1, pnum / pden, (1 - pnum) / (1 - pden))
    return d, den_fit, {"ps_formula": formula, "formula": formula, "used_vars": used, "num_method": num_method, "den_method": den_method, "standardization_table": std_table}


def effective_sample_size(weights) -> float:
    w = np.asarray(weights, dtype=float); w = w[np.isfinite(w)]
    return np.nan if len(w) == 0 or np.sum(w * w) == 0 else float((w.sum() ** 2) / np.sum(w * w))


def weight_summary(df: pd.DataFrame, weight_col: str = "SW_A", treatment_col: str = "A") -> pd.DataFrame:
    rows = []
    for label, sub in [("overall", df), ("deescalated_stopped", df.loc[df[treatment_col] == 1]), ("continued_broad", df.loc[df[treatment_col] == 0])]:
        w = pd.to_numeric(sub[weight_col], errors="coerce").dropna()
        rows.append({"group": label, "n": len(sub), "effective_sample_size": effective_sample_size(w), "mean_weight": w.mean(), "median_weight": w.median(), "min_weight": w.min(), "p01_weight": w.quantile(.01), "p05_weight": w.quantile(.05), "p95_weight": w.quantile(.95), "p99_weight": w.quantile(.99), "max_weight": w.max()})
    return pd.DataFrame(rows)


def risks(df: pd.DataFrame, outcome: str, weight_col: str | None = None, treatment_col: str = "A") -> tuple[float, float, float, float]:
    t = df[treatment_col] == 1; c = df[treatment_col] == 0
    rt = float(df.loc[t, outcome].mean()) if weight_col is None else weighted_mean(df.loc[t, outcome], df.loc[t, weight_col]); rc = float(df.loc[c, outcome].mean()) if weight_col is None else weighted_mean(df.loc[c, outcome], df.loc[c, weight_col])
    return rt, rc, rt - rc, rt / rc if rc > 0 else np.nan


def mean_difference(df: pd.DataFrame, outcome: str, weight_col: str | None = None, treatment_col: str = "A") -> tuple[float, float, float]:
    t = df[treatment_col] == 1; c = df[treatment_col] == 0
    mt = float(df.loc[t, outcome].mean()) if weight_col is None else weighted_mean(df.loc[t, outcome], df.loc[t, weight_col]); mc = float(df.loc[c, outcome].mean()) if weight_col is None else weighted_mean(df.loc[c, outcome], df.loc[c, weight_col])
    return mt, mc, mt - mc


def balance_table(df: pd.DataFrame, variables: Iterable[str], treatment_col: str = "A", weight_col: str = "SW_A") -> pd.DataFrame:
    rows = []; t = df[treatment_col] == 1; c = df[treatment_col] == 0
    for var in variables:
        if var not in df: continue
        x = pd.to_numeric(df[var], errors="coerce")
        if x.isna().all() or x.nunique(dropna=True) < 2: continue
        if _is_binary_01(x):
            xx = x.fillna(0); fn = smd_binary
        else:
            xx = x.fillna(x.mean()); fn = smd_continuous
        before = fn(xx[t], xx[c]); after = fn(xx[t], xx[c], df.loc[t, weight_col], df.loc[c, weight_col])
        rows.append({"variable": var, "before": abs(before), "after": abs(after), "signed_before": before, "signed_after": after})
    return pd.DataFrame(rows).sort_values("before", ascending=False).reset_index(drop=True) if rows else pd.DataFrame(columns=["variable", "before", "after", "signed_before", "signed_after"])


def bootstrap_iptw_ci(df: pd.DataFrame, candidate_vars: Iterable[str], outcome: str, kind: str, reps: int, seed: int, treatment_col: str = "A") -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed); records = []; failures = 0; n = len(df)
    for b in range(reps):
        sample = df.iloc[rng.integers(0, n, size=n)].copy()
        if sample[treatment_col].nunique() < 2: failures += 1; continue
        try:
            w, _, _ = fit_stabilized_iptw(sample, candidate_vars, treatment_col=treatment_col)
            if kind == "risk":
                rt, rc, rd, rr = risks(w, outcome, "SW_A", treatment_col); records.append({"rep": b, "risk_treated": rt, "risk_control": rc, "risk_difference": rd, "risk_ratio": rr})
            elif kind == "mean":
                mt, mc, md = mean_difference(w, outcome, "SW_A", treatment_col); records.append({"rep": b, "mean_treated": mt, "mean_control": mc, "mean_difference": md})
            else: raise ValueError(f"Unknown bootstrap kind: {kind}")
        except Exception: failures += 1
    boot = pd.DataFrame(records)
    if boot.empty: return pd.DataFrame(), boot
    ci_rows = []
    for est in [c for c in boot.columns if c != "rep"]:
        vals = pd.to_numeric(boot[est], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        ci_rows.append({"estimand": est, "bootstrap_successes": len(vals), "bootstrap_failures": failures, "n_success": len(vals), "n_requested": reps, "lower_95": vals.quantile(.025), "upper_95": vals.quantile(.975)})
    return pd.DataFrame(ci_rows), boot


def truncate_weights(df: pd.DataFrame, low: float, high: float, weight_col: str = "SW_A") -> pd.DataFrame:
    d = df.copy(); lo = d[weight_col].quantile(low / 100.0); hi = d[weight_col].quantile(high / 100.0); d[f"{weight_col}_trunc_{low:g}_{high:g}"] = d[weight_col].clip(lo, hi); return d

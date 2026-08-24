from __future__ import annotations

import time
import warnings
from typing import Sequence

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit

from .stats import prepare_ps_covariates, risks


def _fit_strategy(
    df: pd.DataFrame,
    variables: Sequence[str],
    glm_maxiter: int,
) -> tuple[np.ndarray, np.ndarray, str, float, float, float]:
    """Fit GLM up to glm_maxiter then use the historical ridge fallback."""
    d, ps_vars, _ = prepare_ps_covariates(df, variables)
    y = pd.to_numeric(d["A"], errors="coerce").to_numpy(dtype=float)
    x = d[ps_vars].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    x = np.column_stack([np.ones(len(d), dtype=float), x])
    model = sm.GLM(y, x, family=sm.families.Binomial())

    started = time.perf_counter()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="overflow encountered in exp", category=RuntimeWarning)
        try:
            fit = model.fit(maxiter=int(glm_maxiter), disp=0)
            if bool(getattr(fit, "converged", True)):
                method = "glm"
            else:
                fit = model.fit_regularized(alpha=0.001, L1_wt=0.0, maxiter=200, opt_method="bfgs")
                method = "ridge_fallback"
        except Exception:
            fit = model.fit_regularized(alpha=0.001, L1_wt=0.0, maxiter=200, opt_method="bfgs")
            method = "ridge_fallback"
    elapsed = time.perf_counter() - started

    params = np.asarray(fit.params, dtype=float)
    pden = np.clip(expit(np.clip(x @ params, -35.0, 35.0)), 0.001, 0.999)
    pnum = float(np.clip(np.mean(y), 0.001, 0.999))
    weights = np.where(y == 1, pnum / pden, (1.0 - pnum) / (1.0 - pden))
    tmp = pd.DataFrame({"A": y, "death_by_horizon": pd.to_numeric(d["death_by_horizon"], errors="coerce").to_numpy(), "SW": weights})
    _, _, rd, rr = risks(tmp, "death_by_horizon", "SW")
    return pden, weights, method, elapsed, rd, rr


def benchmark_capped_glm(
    cohort: pd.DataFrame,
    variables: Sequence[str],
    reps: int = 20,
    seed: int = 20260426,
    caps: Sequence[int] = (10, 20, 30, 50),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare capped GLM attempts with the historical 200-iteration strategy.

    The goal is to find a lower iteration cap that reproduces the historical
    estimator exactly: samples that converge early remain ordinary GLM; samples
    that do not converge by the cap go directly to the same ridge fallback.
    """
    rng = np.random.default_rng(seed)
    detail: list[dict] = []
    n = len(cohort)
    for rep in range(reps):
        idx = rng.integers(0, n, size=n, dtype=np.int32)
        sample = cohort.iloc[idx].copy()
        if sample["A"].nunique() != 2:
            continue

        ref_ps, ref_w, ref_method, ref_sec, ref_rd, ref_rr = _fit_strategy(sample, variables, 200)
        for cap in caps:
            ps, w, method, sec, rd, rr = _fit_strategy(sample, variables, int(cap))
            detail.append(
                {
                    "rep": rep,
                    "cap": int(cap),
                    "reference_method": ref_method,
                    "capped_method": method,
                    "reference_seconds": ref_sec,
                    "capped_seconds": sec,
                    "speedup": ref_sec / sec if sec > 0 else np.nan,
                    "reference_rd": ref_rd,
                    "capped_rd": rd,
                    "rd_difference": rd - ref_rd,
                    "reference_rr": ref_rr,
                    "capped_rr": rr,
                    "rr_difference": rr - ref_rr,
                    "max_abs_ps_difference": float(np.max(np.abs(ref_ps - ps))),
                    "max_abs_weight_difference": float(np.max(np.abs(ref_w - w))),
                    "method_match": bool(ref_method == method),
                }
            )

    df = pd.DataFrame(detail)
    summaries: list[dict] = []
    for cap, sub in df.groupby("cap", sort=True):
        summaries.append(
            {
                "cap": int(cap),
                "reps": int(len(sub)),
                "median_reference_seconds": float(sub["reference_seconds"].median()),
                "median_capped_seconds": float(sub["capped_seconds"].median()),
                "median_speedup": float(sub["speedup"].median()),
                "method_matches": int(sub["method_match"].sum()),
                "max_abs_rd_difference": float(sub["rd_difference"].abs().max()),
                "median_abs_rd_difference": float(sub["rd_difference"].abs().median()),
                "max_abs_rr_difference": float(sub["rr_difference"].abs().max()),
                "max_abs_ps_difference": float(sub["max_abs_ps_difference"].max()),
                "max_abs_weight_difference": float(sub["max_abs_weight_difference"].max()),
            }
        )
    return df, pd.DataFrame(summaries)

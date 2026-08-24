from __future__ import annotations

import time
import warnings
from typing import Sequence

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit

from .stats import prepare_ps_covariates, risks


def _fit_weights(df: pd.DataFrame, variables: Sequence[str], force_ridge: bool) -> tuple[pd.DataFrame, str, float]:
    d, ps_vars, _ = prepare_ps_covariates(df, variables)
    # Consolidate the prepared frame before adding benchmark-only columns. The
    # PS preparation intentionally adds many standardized columns one by one and
    # can leave pandas with a highly fragmented internal block layout.
    d = d.copy()
    y = pd.to_numeric(d["A"], errors="coerce").to_numpy(dtype=float)
    x = d[ps_vars].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    x = np.column_stack([np.ones(len(d), dtype=float), x])
    model = sm.GLM(y, x, family=sm.families.Binomial())
    started = time.perf_counter()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="overflow encountered in exp", category=RuntimeWarning)
        if force_ridge:
            fit = model.fit_regularized(alpha=0.001, L1_wt=0.0, maxiter=200, opt_method="bfgs")
            method = "ridge"
        else:
            try:
                fit = model.fit(maxiter=200, disp=0)
                if not bool(getattr(fit, "converged", True)):
                    fit = model.fit_regularized(alpha=0.001, L1_wt=0.0, maxiter=200, opt_method="bfgs")
                    method = "ridge_fallback"
                else:
                    method = "glm"
            except Exception:
                fit = model.fit_regularized(alpha=0.001, L1_wt=0.0, maxiter=200, opt_method="bfgs")
                method = "ridge_fallback"
    elapsed = time.perf_counter() - started
    params = np.asarray(fit.params, dtype=float)
    pden = np.clip(expit(np.clip(x @ params, -35.0, 35.0)), 0.001, 0.999)
    pnum = float(np.clip(np.mean(y), 0.001, 0.999))
    weights = np.where(y == 1, pnum / pden, (1.0 - pnum) / (1.0 - pden))
    d = pd.concat(
        [
            d,
            pd.DataFrame(
                {"ps_den_bench": pden, "SW_bench": weights},
                index=d.index,
            ),
        ],
        axis=1,
    )
    return d, method, elapsed


def benchmark_ridge_vs_auto(
    cohort: pd.DataFrame,
    variables: Sequence[str],
    reps: int = 20,
    seed: int = 20260426,
) -> pd.DataFrame:
    """Compare historical auto GLM/fallback with direct ridge on identical resamples.

    This is a diagnostic only. It does not change the production estimator.
    """
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    n = len(cohort)
    for rep in range(reps):
        idx = rng.integers(0, n, size=n, dtype=np.int32)
        s = cohort.iloc[idx].copy()
        if s["A"].nunique() != 2:
            continue
        auto, auto_method, auto_seconds = _fit_weights(s, variables, force_ridge=False)
        ridge, _, ridge_seconds = _fit_weights(s, variables, force_ridge=True)
        _, _, auto_rd, auto_rr = risks(auto, "death_by_horizon", "SW_bench")
        _, _, ridge_rd, ridge_rr = risks(ridge, "death_by_horizon", "SW_bench")
        rows.append({
            "rep": rep,
            "auto_method": auto_method,
            "auto_seconds": auto_seconds,
            "ridge_seconds": ridge_seconds,
            "speedup": auto_seconds / ridge_seconds if ridge_seconds > 0 else np.nan,
            "auto_rd": auto_rd,
            "ridge_rd": ridge_rd,
            "rd_difference": ridge_rd - auto_rd,
            "auto_rr": auto_rr,
            "ridge_rr": ridge_rr,
            "rr_difference": ridge_rr - auto_rr,
            "max_abs_ps_difference": float(np.max(np.abs(auto["ps_den_bench"].to_numpy() - ridge["ps_den_bench"].to_numpy()))),
            "max_abs_weight_difference": float(np.max(np.abs(auto["SW_bench"].to_numpy() - ridge["SW_bench"].to_numpy()))),
        })
    return pd.DataFrame(rows)

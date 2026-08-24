from __future__ import annotations

import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from multiprocessing import get_context
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit
from threadpoolctl import threadpool_limits

from .stats import mean_difference, prepare_ps_covariates, risks


# Benchmarked against the historical 200-iteration rule on matched MIMIC
# bootstrap samples. A cap of 150 reproduced the historical fit method,
# propensity scores, weights, RD, and RR exactly in the validation benchmark,
# while lower caps (<=125) did not. Ridge fallback remains unchanged.
BOOTSTRAP_GLM_MAXITER = 150


@dataclass(frozen=True)
class OutcomeSpec:
    label: str
    column: str
    kind: str  # "risk" or "mean"


_BOOT_DF: pd.DataFrame | None = None
_BOOT_VARS: list[str] = []
_BOOT_OUTCOMES: list[OutcomeSpec] = []
_BOOT_TRUNCATIONS: list[tuple[float, float]] = []


def _set_single_thread_env() -> None:
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"


def _init_worker(
    df: pd.DataFrame,
    variables: Sequence[str],
    outcomes: Sequence[OutcomeSpec] | None = None,
    truncations: Sequence[Sequence[float]] | None = None,
) -> None:
    global _BOOT_DF, _BOOT_VARS, _BOOT_OUTCOMES, _BOOT_TRUNCATIONS
    _set_single_thread_env()
    _BOOT_DF = df
    _BOOT_VARS = list(variables)
    _BOOT_OUTCOMES = list(outcomes or [])
    _BOOT_TRUNCATIONS = [(float(x[0]), float(x[1])) for x in (truncations or [])]


def resolve_jobs(jobs: int | str | None) -> int:
    if jobs is None or str(jobs).lower() == "auto":
        cpu = os.cpu_count() or 2
        return max(1, min(cpu - 1 if cpu > 1 else 1, 8))
    return max(1, int(jobs))


def fit_stabilized_iptw_fast(
    cohort: pd.DataFrame,
    candidate_vars: Iterable[str],
    treatment_col: str = "A",
) -> tuple[pd.DataFrame, dict]:
    """Numerically stable matrix-form equivalent of the scripted PS fit.

    The same covariate preparation, standardization, clipping and ridge fallback
    used by the primary implementation are preserved. Formula/Patsy parsing and
    the intercept-only numerator GLM are avoided for speed. The stabilized
    numerator probability is exactly the observed treatment prevalence.

    For bootstrap execution only, the ordinary GLM attempt is capped at
    BOOTSTRAP_GLM_MAXITER. The current value (150) was selected only after a
    matched-sample parity benchmark showed exact agreement with the historical
    200-iteration rule; lower caps failed that parity check.
    """
    d, ps_vars, std_table = prepare_ps_covariates(cohort, candidate_vars)
    y = pd.to_numeric(d[treatment_col], errors="coerce").to_numpy(dtype=float)
    if len(np.unique(y[np.isfinite(y)])) != 2:
        raise ValueError("Treatment must have two levels.")

    if ps_vars:
        x = d[ps_vars].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        x = np.column_stack([np.ones(len(d), dtype=float), x])
    else:
        x = np.ones((len(d), 1), dtype=float)

    model = sm.GLM(y, x, family=sm.families.Binomial())
    method = "glm_matrix"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="overflow encountered in exp", category=RuntimeWarning)
        try:
            fit = model.fit(maxiter=BOOTSTRAP_GLM_MAXITER, disp=0)
            if not bool(getattr(fit, "converged", True)):
                fit = model.fit_regularized(alpha=0.001, L1_wt=0.0, maxiter=200)
                method = "regularized_glm_matrix"
        except Exception:
            fit = model.fit_regularized(alpha=0.001, L1_wt=0.0, maxiter=200)
            method = "regularized_glm_matrix"

    params = np.asarray(fit.params, dtype=float)
    pden = expit(np.clip(x @ params, -35.0, 35.0))
    pden = np.clip(pden, 0.001, 0.999)
    pnum_scalar = float(np.clip(np.nanmean(y), 0.001, 0.999))
    pnum = np.full(len(d), pnum_scalar, dtype=float)
    d["ps_den"] = pden
    d["ps_num"] = pnum
    d["SW_A"] = np.where(y == 1, pnum / pden, (1.0 - pnum) / (1.0 - pden))
    return d, {
        "used_vars": list(ps_vars),
        "den_method": method,
        "standardization_table": std_table,
        "glm_maxiter": BOOTSTRAP_GLM_MAXITER,
    }


def _bootstrap_index_batches(n: int, reps: int, seed: int) -> list[tuple[int, np.ndarray]]:
    rng = np.random.default_rng(seed)
    return [(rep, rng.integers(0, n, size=n, dtype=np.int32)) for rep in range(reps)]


def _run_tasks(worker, tasks, jobs: int, initargs: tuple) -> list:
    if jobs <= 1:
        _init_worker(*initargs)
        return [worker(task) for task in tasks]
    try:
        ctx = get_context("fork")
    except ValueError:
        ctx = get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=jobs,
        mp_context=ctx,
        initializer=_init_worker,
        initargs=initargs,
    ) as ex:
        return list(ex.map(worker, tasks, chunksize=1))


def _multi_outcome_worker(task: tuple[int, np.ndarray]) -> tuple[list[dict], str | None]:
    rep, idx = task
    if _BOOT_DF is None:
        raise RuntimeError("Bootstrap worker was not initialized")
    s = _BOOT_DF.iloc[idx].copy()
    if s["A"].nunique() != 2:
        return [], "single_treatment_level"
    try:
        with threadpool_limits(limits=1):
            w, diag = fit_stabilized_iptw_fast(s, _BOOT_VARS)
        rows: list[dict] = []
        for spec in _BOOT_OUTCOMES:
            if spec.kind == "risk":
                rt, rc, rd, rr = risks(w, spec.column, "SW_A")
                rows.append({
                    "rep": rep,
                    "analysis": spec.label,
                    "outcome": spec.column,
                    "kind": "risk",
                    "risk_treated": rt,
                    "risk_control": rc,
                    "risk_difference": rd,
                    "risk_ratio": rr,
                    "ps_method": diag["den_method"],
                })
            elif spec.kind == "mean":
                mt, mc, md = mean_difference(w, spec.column, "SW_A")
                rows.append({
                    "rep": rep,
                    "analysis": spec.label,
                    "outcome": spec.column,
                    "kind": "mean",
                    "mean_treated": mt,
                    "mean_control": mc,
                    "mean_difference": md,
                    "ps_method": diag["den_method"],
                })
            else:
                raise ValueError(f"Unknown outcome kind: {spec.kind}")
        return rows, None
    except Exception as exc:
        return [], type(exc).__name__


def bootstrap_multi_outcome_iptw(
    df: pd.DataFrame,
    candidate_vars: Sequence[str],
    outcomes: Sequence[OutcomeSpec],
    reps: int,
    seed: int,
    jobs: int | str | None = "auto",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Bootstrap several outcomes with one PS fit per replicate."""
    n_jobs = resolve_jobs(jobs)
    tasks = _bootstrap_index_batches(len(df), reps, seed)
    results = _run_tasks(
        _multi_outcome_worker,
        tasks,
        n_jobs,
        (df, list(candidate_vars), list(outcomes), []),
    )
    rows = [row for rep_rows, _ in results for row in rep_rows]
    boot = pd.DataFrame(rows)
    failures = [err for _, err in results if err is not None]
    methods = boot["ps_method"].value_counts().to_dict() if len(boot) and "ps_method" in boot else {}
    diag = pd.DataFrame([
        {
            "n_requested": reps,
            "n_successful_replicates": reps - len(failures),
            "n_failed_replicates": len(failures),
            "jobs": n_jobs,
            "glm_matrix_replicates": int(methods.get("glm_matrix", 0) / max(1, len(outcomes))),
            "regularized_glm_matrix_replicates": int(methods.get("regularized_glm_matrix", 0) / max(1, len(outcomes))),
            "failure_types": ";".join(sorted(set(failures))),
            "glm_maxiter": BOOTSTRAP_GLM_MAXITER,
        }
    ])
    return boot, diag


def bootstrap_ci_long(boot: pd.DataFrame, requested: int) -> pd.DataFrame:
    rows: list[dict] = []
    if boot.empty:
        return pd.DataFrame(columns=["analysis", "estimand", "lower_95", "upper_95", "n_success", "n_requested"])
    for label, sub in boot.groupby("analysis", sort=False):
        kind = str(sub["kind"].iloc[0]) if "kind" in sub else "risk"
        estimands = ["risk_treated", "risk_control", "risk_difference", "risk_ratio"] if kind == "risk" else ["mean_treated", "mean_control", "mean_difference"]
        for est in estimands:
            vals = pd.to_numeric(sub.get(est), errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
            rows.append({
                "analysis": label,
                "estimand": est,
                "lower_95": vals.quantile(0.025),
                "upper_95": vals.quantile(0.975),
                "n_success": len(vals),
                "n_requested": requested,
            })
    return pd.DataFrame(rows)


def _weighting_worker(task: tuple[int, np.ndarray]) -> tuple[list[dict], str | None]:
    rep, idx = task
    if _BOOT_DF is None:
        raise RuntimeError("Bootstrap worker was not initialized")
    s = _BOOT_DF.iloc[idx].copy()
    if s["A"].nunique() != 2:
        return [], "single_treatment_level"
    try:
        with threadpool_limits(limits=1):
            w, diag = fit_stabilized_iptw_fast(s, _BOOT_VARS)
        ps = pd.to_numeric(w["ps_den"], errors="coerce").clip(0.001, 0.999)
        a = pd.to_numeric(w["A"], errors="coerce")
        w["OW_A"] = np.where(a == 1, 1.0 - ps, ps)
        rows: list[dict] = []
        rt, rc, rd, rr = risks(w, "death_by_horizon", "OW_A")
        rows.append({
            "rep": rep,
            "analysis": "Overlap weighting",
            "risk_treated": rt,
            "risk_control": rc,
            "risk_difference": rd,
            "risk_ratio": rr,
            "ps_method": diag["den_method"],
        })
        for low, high in _BOOT_TRUNCATIONS:
            col = "SW_A_truncated"
            tmp = w.copy()
            lo = float(tmp["SW_A"].quantile(low / 100.0))
            hi = float(tmp["SW_A"].quantile(high / 100.0))
            tmp[col] = tmp["SW_A"].clip(lo, hi)
            rt, rc, rd, rr = risks(tmp, "death_by_horizon", col)
            rows.append({
                "rep": rep,
                "analysis": f"IPTW truncated {low:g}/{high:g}",
                "risk_treated": rt,
                "risk_control": rc,
                "risk_difference": rd,
                "risk_ratio": rr,
                "ps_method": diag["den_method"],
            })
        return rows, None
    except Exception as exc:
        return [], type(exc).__name__


def bootstrap_weighting_strategies(
    df: pd.DataFrame,
    candidate_vars: Sequence[str],
    truncation_percentiles: Sequence[Sequence[float]],
    reps: int,
    seed: int,
    jobs: int | str | None = "auto",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit the PS once per bootstrap replicate and evaluate all weight strategies."""
    n_jobs = resolve_jobs(jobs)
    tasks = _bootstrap_index_batches(len(df), reps, seed)
    results = _run_tasks(
        _weighting_worker,
        tasks,
        n_jobs,
        (df, list(candidate_vars), [], list(truncation_percentiles)),
    )
    rows = [row for rep_rows, _ in results for row in rep_rows]
    boot = pd.DataFrame(rows)
    failures = [err for _, err in results if err is not None]
    methods = boot["ps_method"].value_counts().to_dict() if len(boot) and "ps_method" in boot else {}
    n_strategies = 1 + len(truncation_percentiles)
    diag = pd.DataFrame([
        {
            "n_requested": reps,
            "n_successful_replicates": reps - len(failures),
            "n_failed_replicates": len(failures),
            "jobs": n_jobs,
            "glm_matrix_replicates": int(methods.get("glm_matrix", 0) / max(1, n_strategies)),
            "regularized_glm_matrix_replicates": int(methods.get("regularized_glm_matrix", 0) / max(1, n_strategies)),
            "failure_types": ";".join(sorted(set(failures))),
            "glm_maxiter": BOOTSTRAP_GLM_MAXITER,
        }
    ])
    return boot, diag

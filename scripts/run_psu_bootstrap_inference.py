#!/usr/bin/env python3
"""Bootstrap inference for the frozen PSU modified external replication.

Why this exists
---------------
The PSU analysis was deliberately frozen in stages: cohort/exposure/covariates first,
then outcomes, then treatment-effect estimation. This script performs inference only;
it must not change any of those scientific definitions.

Implementation note
-------------------
The primary point-estimate script reuses the already-audited propensity-score source
verbatim and appends the frozen outcome/effect block after the PS is fit. We keep the
same parity-preserving mechanism here so the bootstrap refits the *same* PS
specification inside every resample instead of maintaining a second hand-copied model.
This is intentionally conservative for reproducibility, even though it is less elegant
than a full software refactor. Any future refactor must first demonstrate numerical
parity with the frozen publication results.

The bootstrap is encounter-level, nonparametric, uses a fixed seed, refits the frozen
PS in every replicate, and recomputes the primary stabilized ATE IPTW plus the
prespecified overlap and weight-truncation sensitivities. Only aggregate outputs are
written; patient identifiers and row-level bootstrap samples never leave the process.
"""
from __future__ import annotations

import json
import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit
from threadpoolctl import threadpool_limits

from run_psu_point_estimates import INJECT as POINT_INJECT

# Publication inference settings. These are fixed rather than CLI-tunable so a reviewer
# cannot accidentally produce a different inferential result with a different seed or
# bootstrap count while believing it is the frozen analysis.
REPS = 1000
SEED = 20260826
GLM_MAXITER = 150
_BOOT_DF = None
_BOOT_BINARY = None
_BOOT_CONTINUOUS = None
_BOOT_OUTCOMES = None


def _init_worker(df, binary, continuous, outcomes):
    global _BOOT_DF, _BOOT_BINARY, _BOOT_CONTINUOUS, _BOOT_OUTCOMES
    _BOOT_DF = df
    _BOOT_BINARY = binary
    _BOOT_CONTINUOUS = continuous
    _BOOT_OUTCOMES = outcomes
    # One BLAS/OpenMP thread per process avoids oversubscription when several bootstrap
    # workers run in parallel. This changes runtime only, not the estimand or result.
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"


def _wmean(z, w):
    m = np.isfinite(z) & np.isfinite(w)
    return float(np.sum(z[m] * w[m]) / np.sum(w[m]))


def _worker(task):
    rep, idx = task
    s = _BOOT_DF.iloc[idx].copy()
    y = s["A"].astype(int).to_numpy()
    if np.unique(y).size != 2:
        return [], "single_treatment_level"

    # Recreate the frozen PS preprocessing *inside each bootstrap sample*.
    # Binary covariates: missing -> 0, then clip to [0,1].
    # Continuous covariates: sample median imputation, standardize within replicate,
    # and clip standardized values to +/-8 for numerical stability.
    X = s[_BOOT_BINARY + _BOOT_CONTINUOUS].copy()
    for col in _BOOT_BINARY:
        X[col] = pd.to_numeric(X[col], errors="coerce").fillna(0).clip(0, 1)
    for col in _BOOT_CONTINUOUS:
        v = pd.to_numeric(X[col], errors="coerce")
        med = float(v.median()) if v.notna().any() else 0.0
        v = v.fillna(med)
        mu = float(v.mean())
        sd = float(v.std(ddof=0))
        X[col] = ((v - mu) / (sd if sd > 0 else 1.0)).clip(-8, 8)

    Xc = sm.add_constant(X, has_constant="add").to_numpy(float)
    method = "glm_matrix"
    try:
        with threadpool_limits(limits=1), warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            model = sm.GLM(y, Xc, family=sm.families.Binomial())
            fit = model.fit(maxiter=GLM_MAXITER, disp=0)
            # This fallback is the same prespecified numerical safeguard used by the
            # frozen analysis. It is used only when the ordinary GLM does not converge.
            if not bool(getattr(fit, "converged", True)):
                fit = model.fit_regularized(alpha=0.001, L1_wt=0.0, maxiter=200)
                method = "regularized_glm_matrix"

        eta = Xc @ np.asarray(fit.params, float)
        ps = np.clip(expit(np.clip(eta, -35, 35)), 0.001, 0.999)
        pa = float(np.mean(y))
        sw = np.where(y == 1, pa / ps, (1 - pa) / (1 - ps))
        ow = np.where(y == 1, 1 - ps, ps)
        q01, q99 = np.quantile(sw, [0.01, 0.99])
        w199 = np.clip(sw, q01, q99)
        q025, q975 = np.quantile(sw, [0.025, 0.975])
        w25975 = np.clip(sw, q025, q975)
        methods = {
            "stabilized_ate_iptw": sw,
            "overlap_weighting": ow,
            "trunc_1_99": w199,
            "trunc_2p5_97p5": w25975,
        }

        rows = []
        for mname, w in methods.items():
            for outcome, kind in _BOOT_OUTCOMES:
                z = pd.to_numeric(s[outcome], errors="coerce").to_numpy(float)
                mu1 = _wmean(z[y == 1], w[y == 1])
                mu0 = _wmean(z[y == 0], w[y == 0])
                rows.append(
                    {
                        "rep": rep,
                        "method": mname,
                        "outcome": outcome,
                        "type": kind,
                        "deescalated_mean_or_risk": mu1,
                        "continued_mean_or_risk": mu0,
                        "difference_A1_minus_A0": mu1 - mu0,
                        "risk_ratio_A1_over_A0": (mu1 / mu0 if kind == "binary" and mu0 > 0 else np.nan),
                        "ps_method": method,
                    }
                )
        return rows, None
    except Exception as exc:
        return [], type(exc).__name__


def run_bootstrap(df, binary, continuous, effect_outcomes, outdir):
    rng = np.random.default_rng(SEED)
    n = len(df)
    tasks = [(r, rng.integers(0, n, size=n, dtype=np.int32)) for r in range(REPS)]
    jobs = max(1, min((os.cpu_count() or 2) - 1, 8))
    try:
        ctx = get_context("fork")
    except ValueError:
        ctx = get_context("spawn")

    with ProcessPoolExecutor(
        max_workers=jobs,
        mp_context=ctx,
        initializer=_init_worker,
        initargs=(df, binary, continuous, effect_outcomes),
    ) as ex:
        results = list(ex.map(_worker, tasks, chunksize=1))

    rows = [row for rr, _ in results for row in rr]
    failures = [e for _, e in results if e is not None]
    boot = pd.DataFrame(rows)
    ci = []
    for (method, outcome), sub in boot.groupby(["method", "outcome"], sort=False):
        kind = str(sub["type"].iloc[0])
        metrics = ["deescalated_mean_or_risk", "continued_mean_or_risk", "difference_A1_minus_A0"]
        if kind == "binary":
            metrics.append("risk_ratio_A1_over_A0")
        for metric in metrics:
            vals = pd.to_numeric(sub[metric], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
            ci.append(
                {
                    "method": method,
                    "outcome": outcome,
                    "estimand": metric,
                    "lower_95": float(vals.quantile(0.025)),
                    "upper_95": float(vals.quantile(0.975)),
                    "n_success": int(len(vals)),
                    "n_requested": REPS,
                }
            )
    pd.DataFrame(ci).to_csv(outdir / "bootstrap_ci.csv", index=False)

    method_counts = boot.drop_duplicates(["rep"])["ps_method"].value_counts().to_dict() if len(boot) else {}
    diag = {
        "n_requested": REPS,
        "n_successful_replicates": REPS - len(failures),
        "n_failed_replicates": len(failures),
        "jobs": jobs,
        "seed": SEED,
        "glm_maxiter": GLM_MAXITER,
        "glm_matrix_replicates": int(method_counts.get("glm_matrix", 0)),
        "regularized_glm_matrix_replicates": int(method_counts.get("regularized_glm_matrix", 0)),
        "failure_types": sorted(set(failures)),
    }
    (outdir / "bootstrap_diagnostics.json").write_text(json.dumps(diag, indent=2))


BOOT_INJECT = r'''
    run_bootstrap(df, binary, continuous, effect_outcomes, args.output_dir)
'''


def main():
    here = Path(__file__).resolve().parent
    source_path = here / "audit_psu_ps_balance.py"
    source = source_path.read_text(encoding="utf-8")
    marker = "    max_pre=float(bdf.abs_pre_smd.max()); max_post=float(bdf.abs_post_smd.max()); worst=str(bdf.iloc[0].variable)\n"
    if source.count(marker) != 1:
        raise RuntimeError("Frozen PS source marker missing or ambiguous")

    # The marker is a fail-closed guardrail: if the frozen PS source is edited in a way
    # that makes the insertion point ambiguous, this script stops rather than silently
    # running a different model.
    source = source.replace(marker, POINT_INJECT + "\n" + BOOT_INJECT + "\n" + marker, 1)
    ns = {"__name__": "__main__", "__file__": str(source_path), "run_bootstrap": run_bootstrap}
    exec(compile(source, str(source_path), "exec"), ns, ns)


if __name__ == "__main__":
    main()

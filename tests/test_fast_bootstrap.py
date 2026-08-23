import numpy as np
import pandas as pd

from sepsis_deescalation.fast_bootstrap import (
    OutcomeSpec,
    bootstrap_multi_outcome_iptw,
    fit_stabilized_iptw_fast,
)
from sepsis_deescalation.stats import fit_stabilized_iptw


def _toy_data(n=300, seed=7):
    rng = np.random.default_rng(seed)
    x1 = rng.normal(size=n)
    x2 = rng.binomial(1, 0.35, size=n)
    p = 1 / (1 + np.exp(-(-0.4 + 0.7 * x1 - 0.3 * x2)))
    a = rng.binomial(1, p)
    y = rng.binomial(1, 1 / (1 + np.exp(-(-1.5 + 0.2 * a + 0.5 * x1))))
    z = 5 + 0.5 * a + x1 + rng.normal(size=n)
    return pd.DataFrame({"A": a, "x1": x1, "x2": x2, "y": y, "z": z})


def test_fast_ps_matches_formula_ps_on_well_behaved_data():
    d = _toy_data()
    slow, _, _ = fit_stabilized_iptw(d, ["x1", "x2"])
    fast, diag = fit_stabilized_iptw_fast(d, ["x1", "x2"])
    assert diag["den_method"] in {"glm_matrix", "regularized_glm_matrix"}
    assert np.allclose(slow["ps_den"], fast["ps_den"], atol=1e-8, rtol=1e-7)
    assert np.allclose(slow["SW_A"], fast["SW_A"], atol=1e-8, rtol=1e-7)


def test_shared_bootstrap_returns_each_outcome_per_successful_rep():
    d = _toy_data()
    outcomes = [OutcomeSpec("binary", "y", "risk"), OutcomeSpec("continuous", "z", "mean")]
    boot, diag = bootstrap_multi_outcome_iptw(d, ["x1", "x2"], outcomes, reps=8, seed=11, jobs=1)
    assert set(boot["analysis"]) == {"binary", "continuous"}
    counts = boot.groupby("analysis")["rep"].nunique().to_dict()
    assert counts["binary"] == counts["continuous"]
    assert int(diag.iloc[0]["n_requested"]) == 8
    assert int(diag.iloc[0]["n_failed_replicates"]) == 0

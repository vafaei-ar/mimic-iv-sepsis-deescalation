import numpy as np
import pandas as pd

from sepsis_deescalation.stats import balance_table, fit_stabilized_iptw, risks, weight_summary


def test_iptw_smoke():
    rng = np.random.default_rng(42)
    n = 400
    x = rng.normal(size=n)
    p = 1 / (1 + np.exp(-0.7 * x))
    a = rng.binomial(1, p)
    y = rng.binomial(1, 1 / (1 + np.exp(-(-1.5 + 0.3 * a + 0.4 * x))))
    d = pd.DataFrame({"A": a, "x": x, "y": y})
    w, _, diag = fit_stabilized_iptw(d, ["x"])
    assert np.isfinite(w["SW_A"]).all()
    assert len(diag["used_vars"]) >= 1
    rt, rc, rd, rr = risks(w, "y", "SW_A")
    assert all(np.isfinite(v) for v in [rt, rc, rd, rr])
    summary = weight_summary(w)
    assert set(summary["group"]) == {"overall", "deescalated_stopped", "continued_broad"}
    bal = balance_table(w, ["x"])
    assert "after" in bal

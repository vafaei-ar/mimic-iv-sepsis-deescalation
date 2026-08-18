import numpy as np
import pandas as pd

from sepsis_deescalation.weighting_audit import add_overlap_weights, run_weighting_audit, weight_tail_diagnostics


def test_overlap_weights_and_audit_smoke():
    rng = np.random.default_rng(7)
    n = 300
    x = rng.normal(size=n)
    ps = 1 / (1 + np.exp(-0.8 * x))
    a = rng.binomial(1, ps)
    y = rng.binomial(1, 0.15 + 0.05 * a)
    p_treated = a.mean()
    sw = np.where(a == 1, p_treated / np.clip(ps, 0.001, 0.999), (1 - p_treated) / (1 - np.clip(ps, 0.001, 0.999)))
    d = pd.DataFrame({"A": a, "x": x, "death_by_horizon": y, "ps_den": ps, "SW_A": sw})

    ow = add_overlap_weights(d)
    assert np.isfinite(ow["OW_A"]).all()
    assert (ow["OW_A"] > 0).all()

    tails = weight_tail_diagnostics(d)
    assert set(tails["group"]) == {"deescalated_stopped", "continued_broad"}

    result = run_weighting_audit(d, ["x"])
    assert set(result["summary"]["analysis"]) == {
        "Primary stabilized IPTW",
        "IPTW truncated 1/99",
        "IPTW truncated 2.5/97.5",
        "Overlap weighting",
    }

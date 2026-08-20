import numpy as np
import pandas as pd

from sepsis_deescalation.weighting_final import add_overlap_weights, run_final_weighting_sensitivities


def test_overlap_weights_are_bounded():
    d = pd.DataFrame({"A": [1, 1, 0, 0], "ps_den": [0.2, 0.8, 0.2, 0.8]})
    out = add_overlap_weights(d)
    assert np.allclose(out["OW_A"], [0.8, 0.2, 0.2, 0.8])
    assert out["OW_A"].between(0, 1).all()


def test_final_weighting_sensitivities_smoke(tmp_path):
    rng = np.random.default_rng(7)
    n = 250
    x = rng.normal(size=n)
    p = 1 / (1 + np.exp(-0.6 * x))
    a = rng.binomial(1, p)
    y = rng.binomial(1, 1 / (1 + np.exp(-(-1.3 + 0.2 * a + 0.3 * x))))
    d = pd.DataFrame({"A": a, "x": x, "death_by_horizon": y})
    result = run_final_weighting_sensitivities(
        d,
        ["x"],
        out_dir=tmp_path,
        truncation_percentiles=((1.0, 99.0),),
        reps=20,
        seed=123,
    )
    assert set(result["summary"]["analysis"]) == {
        "Primary stabilized IPTW",
        "Overlap weighting",
        "IPTW truncated 1/99",
    }
    assert not result["ci"].empty
    assert (tmp_path / "final_weighting_point_estimates.csv").exists()

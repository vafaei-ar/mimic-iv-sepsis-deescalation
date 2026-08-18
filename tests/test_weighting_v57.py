import numpy as np
import pandas as pd

from sepsis_deescalation.weighting_v57 import add_overlap_weights, deduplicate_covariates


def test_deduplicate_covariates_keeps_first_exact_duplicate():
    d = pd.DataFrame({"x": [0.0, 1.0, 2.0], "x_copy": [0.0, 1.0, 2.0], "z": [1.0, 0.0, 1.0]})
    kept, removed = deduplicate_covariates(d, ["x", "x_copy", "z"])
    assert kept == ["x", "z"]
    assert removed.iloc[0]["removed_variable"] == "x_copy"
    assert removed.iloc[0]["retained_variable"] == "x"


def test_overlap_weights_are_bounded_and_treatment_specific():
    d = pd.DataFrame({"A": [1, 1, 0, 0], "ps_den": [0.2, 0.8, 0.2, 0.8]})
    out = add_overlap_weights(d)
    assert np.allclose(out["OW_A"].to_numpy(), [0.8, 0.2, 0.2, 0.8])
    assert ((out["OW_A"] >= 0) & (out["OW_A"] <= 1)).all()

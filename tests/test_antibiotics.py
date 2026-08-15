import pandas as pd

from sepsis_deescalation.antibiotics import prepare_coverage, systemic_antibiotic_mask
from sepsis_deescalation.pcornet_antibiotics import classify_broad_drug, is_observed_systemic_antibiotic


def test_mimic_broad_requires_iv():
    d = pd.DataFrame({"drug_lower": ["vancomycin", "vancomycin", "cefepime"], "route_lower": ["iv", "po", "intravenous"]})
    mask = systemic_antibiotic_mask(d, broad=True)
    assert mask.tolist() == [True, False, True]


def test_mimic_missing_stop_fill():
    d = pd.DataFrame({"starttime": pd.to_datetime(["2020-01-01 00:00"]), "stoptime": [pd.NaT]})
    out = prepare_coverage(d, 24)
    assert (out.loc[0, "coverage_stop"] - out.loc[0, "coverage_start"]).total_seconds() == 24 * 3600


def test_pcornet_broad_mapping_and_non_iv_exclusion():
    assert classify_broad_drug("20481", "cefepime", "UN") == "cefepime"
    assert classify_broad_drug("313570", "vancomycin capsule", "UN") is None
    assert is_observed_systemic_antibiotic("", "ceftriaxone injection", "UN")
    assert not is_observed_systemic_antibiotic("", "vancomycin oral capsule", "ORAL")

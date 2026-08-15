import pandas as pd

from sepsis_deescalation.microbiology import eligible_microbiology


def _cohort():
    return pd.DataFrame({
        "hadm_id": [1, 2, 3],
        "first_broad_time": pd.to_datetime(["2020-01-01"] * 3),
        "decision_time": pd.to_datetime(["2020-01-04"] * 3),
    })


def test_result_available_by_decision_is_primary_rule():
    micro = pd.DataFrame({
        "hadm_id": [1, 2, 3],
        "subject_id": [11, 12, 13],
        "specimen_time": pd.to_datetime(["2020-01-02"] * 3),
        "result_available_time": pd.to_datetime(["2020-01-03", "2020-01-05", "2020-01-03"]),
        "org_name": ["E coli", "E coli", ""],
        "spec_type_desc": ["BLOOD CULTURE"] * 3,
        "test_name": ["Blood Culture"] * 3,
        "spec_type_desc_lower": ["blood culture"] * 3,
        "test_name_lower": ["blood culture"] * 3,
    })
    eligible, sets, _ = eligible_microbiology(_cohort(), micro)
    assert set(eligible["hadm_id"]) == {2, 3}
    assert sets["positive_available_by_72"] == {1}
    assert sets["positive_eventual"] == {1, 2}
    assert sets["eventual_culture_negative"] == {3}

from __future__ import annotations

import pandas as pd

SCREEN_PATTERN = (
    r"screen|surveillance|mrsa|vre|rectal|nares|nasal|stool ova|parasite|"
    r"viral culture|covid|sars|influenza|respiratory viral|rsv|staph aureus swab|"
    r"legionella urinary antigen|strep pneumo antigen|clostridioides|c difficile|c\. difficile|toxin|antigen"
)
NON_CULTURE_TEST_PATTERN = (
    r"pcr|polymerase|antigen|antibody|serolog|immunolog|igm|igg|ifa|eia|elisa|"
    r"gram stain|smear|acid fast smear|afb smear|koh|potassium hydroxide|"
    r"rapid|screen|panel|probe|assay|toxin|immunofluorescent|immunofluorescence|dfa|"
    r"pneumocystis|direct fluorescent|molecular"
)
CULTURE_TEST_NAME_PATTERN = (
    r"culture|cultures|blood/fungal culture|blood/afb culture|"
    r"fluid culture|urine culture|blood culture|respiratory culture|wound culture|"
    r"fungal culture|anaerobic culture|acid fast culture|legionella culture|"
    r"campylobacter culture|fecal culture|catheter tip culture|tissue culture|"
    r"aerobic culture|bronchial culture|sputum culture|csf culture"
)
CLINICAL_SPECIMEN_PATTERN = (
    r"blood|urine|sputum|tracheal|bronch|lavage|respiratory|wound|abscess|fluid|"
    r"pleural|peritoneal|csf|cerebrospinal|tissue|catheter|bile|synovial|sterile"
)
FLORA_PATTERN = r"normal flora|mixed bacterial flora|skin flora|mixed flora"


def classify_microbiology(micro: pd.DataFrame) -> pd.DataFrame:
    """Replicate the v5.5 primary and strict microbiology classifier."""
    d = micro.copy()
    spec = d["spec_type_desc_lower"].fillna("")
    test = d["test_name_lower"].fillna("")
    org = d["org_name"].fillna("").str.strip()
    text = spec + " " + test

    d["screen_or_surveillance"] = text.str.contains(SCREEN_PATTERN, na=False, regex=True).astype(int)
    d["non_culture_test"] = test.str.contains(NON_CULTURE_TEST_PATTERN, na=False, regex=True).astype(int)
    d["culture_test_name"] = test.str.contains(CULTURE_TEST_NAME_PATTERN, na=False, regex=True).astype(int)
    d["clinical_specimen_type"] = spec.str.contains(CLINICAL_SPECIMEN_PATTERN, na=False, regex=True).astype(int)
    d["clinical_micro"] = (d["screen_or_surveillance"] == 0).astype(int)
    d["true_culture_micro"] = (
        (d["clinical_micro"] == 1) & (d["culture_test_name"] == 1) & (d["non_culture_test"] == 0)
    ).astype(int)

    d["positive_clinical_culture"] = ((d["clinical_micro"] == 1) & org.ne("")).astype(int)
    d["positive_true_culture"] = ((d["true_culture_micro"] == 1) & org.ne("")).astype(int)
    flora = d["org_name"].fillna("").str.lower().str.contains(FLORA_PATTERN, na=False, regex=True)
    d["positive_clinical_culture_no_flora"] = ((d["positive_clinical_culture"] == 1) & ~flora).astype(int)
    d["positive_true_culture_no_flora"] = ((d["positive_true_culture"] == 1) & ~flora).astype(int)
    d["positive_organism_row"] = d["positive_clinical_culture"]
    d["result_time_missing_for_positive_row"] = (
        (d["positive_organism_row"] == 1) & d["result_available_time"].isna()
    ).astype(int)

    d["strict_culture_exclusion_reason"] = "not_excluded_from_strict_culture"
    d.loc[d["screen_or_surveillance"] == 1, "strict_culture_exclusion_reason"] = "screen_or_surveillance"
    d.loc[(d["screen_or_surveillance"] == 0) & (d["non_culture_test"] == 1), "strict_culture_exclusion_reason"] = "non_culture_diagnostic_test"
    d.loc[(d["screen_or_surveillance"] == 0) & (d["non_culture_test"] == 0) & (d["culture_test_name"] == 0), "strict_culture_exclusion_reason"] = "test_name_not_culture"
    d["micro_category"] = "clinical_non_culture_or_uncertain"
    d.loc[d["screen_or_surveillance"] == 1, "micro_category"] = "screen_or_surveillance_excluded"
    d.loc[(d["clinical_micro"] == 1) & (d["non_culture_test"] == 1), "micro_category"] = "non_culture_test_excluded_from_strict_sensitivity"
    d.loc[(d["clinical_micro"] == 1) & (d["non_culture_test"] == 0) & (d["culture_test_name"] == 0), "micro_category"] = "clinical_record_test_name_not_culture"
    d.loc[d["true_culture_micro"] == 1, "micro_category"] = "strict_test_name_culture_included"
    d["strict_culture_record"] = d["true_culture_micro"]
    d["strict_positive_culture"] = d["positive_true_culture"]
    return d


def eligible_microbiology(cohort: pd.DataFrame, micro: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, set[int]], pd.DataFrame]:
    """Primary eligibility uses specimen time for sampling and result-availability time for positivity."""
    d = classify_microbiology(micro)
    d = d.merge(cohort[["hadm_id", "first_broad_time", "decision_time"]], on="hadm_id", how="inner")
    d["micro_time"] = d["specimen_time"]
    win = d.loc[
        d["specimen_time"].notna()
        & (d["specimen_time"] >= d["first_broad_time"] - pd.Timedelta(hours=24))
        & (d["specimen_time"] <= d["decision_time"])
        & (d["clinical_micro"] == 1)
    ].copy()

    sampled = set(win["hadm_id"].dropna().astype(int))
    available_positive = set(win.loc[
        (win["positive_clinical_culture"] == 1)
        & win["result_available_time"].notna()
        & (win["result_available_time"] <= win["decision_time"]), "hadm_id"
    ].dropna().astype(int))
    eventual_positive = set(win.loc[win["positive_clinical_culture"] == 1, "hadm_id"].dropna().astype(int))
    missing_result_time_positive = set(win.loc[win["result_time_missing_for_positive_row"] == 1, "hadm_id"].dropna().astype(int))

    strict_win = win.loc[win["true_culture_micro"] == 1].copy()
    strict_sampled = set(strict_win["hadm_id"].dropna().astype(int))
    strict_available_positive = set(strict_win.loc[
        (strict_win["positive_true_culture"] == 1)
        & strict_win["result_available_time"].notna()
        & (strict_win["result_available_time"] <= strict_win["decision_time"]), "hadm_id"
    ].dropna().astype(int))

    out = cohort.loc[cohort["hadm_id"].astype(int).isin(sampled - available_positive)].copy()
    sets = {
        "sampled": sampled,
        "positive_available_by_72": available_positive,
        "positive_eventual": eventual_positive,
        "positive_missing_result_time": missing_result_time_positive,
        "strict_sampled": strict_sampled,
        "strict_positive_available_by_72": strict_available_positive,
        "strict_eligible": strict_sampled - strict_available_positive,
        "eventual_culture_negative": sampled - eventual_positive,
    }
    return out, sets, win


def audit_counts(sets: dict[str, set[int]]) -> pd.DataFrame:
    return pd.DataFrame([{"metric": key, "n_admissions": len(value)} for key, value in sets.items()])

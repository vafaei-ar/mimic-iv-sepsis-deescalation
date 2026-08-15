from __future__ import annotations

import pandas as pd

SCREEN_PATTERN = r"screen|surveillance|mrsa|vre|rectal|nares|nasal|stool ova|parasite|viral culture"
STRICT_CULTURE_TEST_PATTERN = (
    r"culture|blood culture|urine culture|respiratory culture|sputum culture|wound culture|"
    r"body fluid culture|tissue culture|csf culture|fungal culture|afb culture|mycobacterial culture"
)
STRICT_EXCLUSION_PATTERN = (
    r"pcr|polymerase|antigen|antibody|serolog|immuno|gram stain|smear|koh|dfa|"
    r"pneumocystis|toxin|assay|panel|probe|screen|surveillance"
)
FLORA_PATTERN = r"normal flora|mixed bacterial flora|skin flora|mixed flora"


def classify_microbiology(micro: pd.DataFrame) -> pd.DataFrame:
    d = micro.copy()
    spec = d["spec_type_desc_lower"].fillna("")
    test = d["test_name_lower"].fillna("")
    combined = spec + " " + test
    screen = combined.str.contains(SCREEN_PATTERN, na=False, regex=True)
    d["clinical_micro"] = (~screen).astype(int)

    org = d["org_name"].fillna("").str.strip()
    positive = (d["clinical_micro"] == 1) & org.ne("")
    d["positive_clinical_culture"] = positive.astype(int)
    d["positive_organism_row"] = positive.astype(int)
    flora = d["org_name"].fillna("").str.lower().str.contains(FLORA_PATTERN, na=False, regex=True)
    d["positive_clinical_culture_no_flora"] = (positive & ~flora).astype(int)
    d["result_time_missing_for_positive_row"] = (positive & d["result_available_time"].isna()).astype(int)

    strict_test = test.str.contains(STRICT_CULTURE_TEST_PATTERN, na=False, regex=True)
    strict_excluded = test.str.contains(STRICT_EXCLUSION_PATTERN, na=False, regex=True)
    d["strict_culture_record"] = ((d["clinical_micro"] == 1) & strict_test & ~strict_excluded).astype(int)
    d["strict_positive_culture"] = ((d["strict_culture_record"] == 1) & org.ne("")).astype(int)
    return d


def eligible_microbiology(
    cohort: pd.DataFrame,
    micro: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, set[int]], pd.DataFrame]:
    """Return primary eligible cohort and microbiology ID sets.

    Sampling is based on specimen time from 24 h before first broad-spectrum exposure through
    the 72-h decision. Primary positivity uses result availability time <= decision. Eventual
    positivity from qualifying specimens is retained for sensitivity analysis.
    """
    d = classify_microbiology(micro)
    d = d.merge(
        cohort[["hadm_id", "first_broad_time", "decision_time"]],
        on="hadm_id",
        how="inner",
    )
    d["micro_time"] = d["specimen_time"]
    win = d.loc[
        d["specimen_time"].notna()
        & (d["specimen_time"] >= d["first_broad_time"] - pd.Timedelta(hours=24))
        & (d["specimen_time"] <= d["decision_time"])
        & (d["clinical_micro"] == 1)
    ].copy()

    sampled = set(win["hadm_id"].dropna().astype(int))
    available_positive = set(
        win.loc[
            (win["positive_organism_row"] == 1)
            & win["result_available_time"].notna()
            & (win["result_available_time"] <= win["decision_time"]),
            "hadm_id",
        ].dropna().astype(int)
    )
    eventual_positive = set(
        win.loc[win["positive_organism_row"] == 1, "hadm_id"].dropna().astype(int)
    )
    missing_result_time_positive = set(
        win.loc[win["result_time_missing_for_positive_row"] == 1, "hadm_id"].dropna().astype(int)
    )

    strict_sampled = set(
        win.loc[win["strict_culture_record"] == 1, "hadm_id"].dropna().astype(int)
    )
    strict_available_positive = set(
        win.loc[
            (win["strict_positive_culture"] == 1)
            & win["result_available_time"].notna()
            & (win["result_available_time"] <= win["decision_time"]),
            "hadm_id",
        ].dropna().astype(int)
    )

    eligible_ids = sampled - available_positive
    out = cohort.loc[cohort["hadm_id"].astype(int).isin(eligible_ids)].copy()
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


def diagnostic_intensity(cohort: pd.DataFrame, micro_win: pd.DataFrame) -> pd.DataFrame:
    d = cohort.copy()
    counts = micro_win.groupby("hadm_id").size()
    clinical_counts = micro_win.loc[micro_win["clinical_micro"] == 1].groupby("hadm_id").size()
    specimens = micro_win.groupby("hadm_id")["spec_type_desc_lower"].nunique()
    tests = micro_win.groupby("hadm_id")["test_name_lower"].nunique()
    d["micro_records_pre72"] = d["hadm_id"].map(counts).fillna(0).astype(float)
    d["clinical_micro_records_pre72"] = d["hadm_id"].map(clinical_counts).fillna(0).astype(float)
    d["micro_specimen_types_pre72"] = d["hadm_id"].map(specimens).fillna(0).astype(float)
    d["micro_test_names_pre72"] = d["hadm_id"].map(tests).fillna(0).astype(float)
    return d


def audit_counts(sets: dict[str, set[int]]) -> pd.DataFrame:
    return pd.DataFrame(
        [{"metric": key, "n_admissions": len(value)} for key, value in sets.items()]
    )

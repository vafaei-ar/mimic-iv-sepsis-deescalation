from __future__ import annotations

CONTINUOUS_VARS = [
    "age", "hours_admit_to_icu", "bmi_pre72",
    "vasopressor_hours_pre72",
    "lactate_last_pre72", "lactate_change_pre72",
    "creatinine_last_pre72", "creatinine_change_pre72",
    "white_blood_cells_last_pre72", "white_blood_cells_change_pre72",
    "hr_max_pre72", "rr_max_pre72", "spo2_min_pre72", "map_min_pre72", "temp_max_pre72",
    "severity_pre72",
    "lactate_early_last_0_24h", "lactate_late_last_48_72h", "lactate_change_early_to_late",
    "wbc_early_last_0_24h", "wbc_late_last_48_72h", "wbc_change_early_to_late",
    "creatinine_early_last_0_24h", "creatinine_late_last_48_72h", "creatinine_change_early_to_late",
    "platelet_late_worst_48_72h", "platelet_change_early_to_late",
    "bilirubin_late_worst_48_72h", "bilirubin_change_early_to_late",
    "heart_rate_48_72h", "resp_rate_48_72h", "spo2_48_72h", "map_48_72h", "temperature_48_72h",
    "gcs_total_48_72h", "fio2_48_72h",
    "urine_output_ml_48_72h", "urine_output_change_pre72",
    "sofa_like_0_24h", "sofa_like_48_72h", "sofa_like_change_pre72",
    "systemic_abx_hours_pre72", "systemic_abx_agents_pre72",
    "broad_abx_hours_pre72", "broad_abx_agents_pre72",
    "micro_records_pre72", "clinical_micro_records_pre72", "strict_culture_records_pre72", "distinct_specimen_types_pre72",
]

BINARY_VARS = [
    "sex_male", "race_white", "comorb", "heart_failure", "chronic_kidney", "vent_proc",
    "micu", "sicu", "cardiac_icu", "neuro_icu",
    "vasopressor_any_pre72",
    "lactate_missing_pre72", "creatinine_missing_pre72", "white_blood_cells_missing_pre72",
    "hr_max_pre72_missing", "rr_max_pre72_missing", "spo2_min_pre72_missing", "map_min_pre72_missing", "temp_max_pre72_missing",
    "lactate_rising_pre72", "wbc_rising_pre72", "creatinine_rising_pre72", "platelet_falling_pre72", "bilirubin_rising_pre72",
    "fever_last12h_pre72", "map_improved_pre72", "spo2_improved_pre72", "temp_improved_pre72", "rr_improved_pre72",
    "vasopressor_any_0_24h", "vasopressor_any_48_72h", "vasopressor_stopped_before_72h",
    "low_urine_output_48_72h", "sofa_like_improved_pre72",
    "anti_mrsa_pre72", "antipseudomonal_pre72", "carbapenem_pre72", "anaerobic_coverage_pre72",
    "repeat_micro_48_72h", "blood_culture_pre72", "respiratory_culture_pre72", "urine_culture_pre72", "sterile_fluid_culture_pre72",
    "steroid_any_pre72", "hydrocortisone_any_pre72",
]

CANDIDATE_PS_VARS = CONTINUOUS_VARS + BINARY_VARS

_BASE = ["age", "sex_male", "race_white", "comorb", "heart_failure", "chronic_kidney", "hours_admit_to_icu", "micu", "sicu", "cardiac_icu", "neuro_icu"]

PROGRESSIVE_MODELS = [
    {"model": "M1 demographics/comorbidity only", "vars": _BASE},
    {
        "model": "M2 + baseline severity/labs",
        "vars": _BASE + ["vent_proc", "vasopressor_any_pre72", "vasopressor_hours_pre72", "lactate_last_pre72", "creatinine_last_pre72", "white_blood_cells_last_pre72", "severity_pre72"],
    },
    {
        "model": "M3 + near-decision clinical status",
        "vars": _BASE + [
            "vent_proc", "vasopressor_any_pre72", "vasopressor_hours_pre72", "lactate_last_pre72", "creatinine_last_pre72", "white_blood_cells_last_pre72", "severity_pre72",
            "lactate_late_last_48_72h", "wbc_late_last_48_72h", "creatinine_late_last_48_72h", "platelet_late_worst_48_72h", "bilirubin_late_worst_48_72h",
            "heart_rate_48_72h", "resp_rate_48_72h", "spo2_48_72h", "map_48_72h", "temperature_48_72h", "urine_output_ml_48_72h", "sofa_like_48_72h",
        ],
    },
    {"model": "M4 + improvement trajectories and intensity", "vars": CANDIDATE_PS_VARS},
]

CLINICAL_IMPROVEMENT_VARS = [
    "lactate_change_early_to_late", "wbc_change_early_to_late", "creatinine_change_early_to_late",
    "platelet_change_early_to_late", "bilirubin_change_early_to_late", "map_improved_pre72", "spo2_improved_pre72",
    "temp_improved_pre72", "rr_improved_pre72", "fever_last12h_pre72", "vasopressor_any_0_24h", "vasopressor_any_48_72h",
    "vasopressor_stopped_before_72h", "sofa_like_0_24h", "sofa_like_48_72h", "sofa_like_change_pre72",
    "sofa_like_improved_pre72", "systemic_abx_agents_pre72", "broad_abx_agents_pre72", "micro_records_pre72", "strict_culture_records_pre72",
]

# PSU / PCORnet harmonization crosswalk

This file must be updated after inspecting the actual Penn State datamart. Do not infer that a standard PCORnet field has the needed clinical semantics simply because it exists.

| Construct | MIMIC-IV | PSU/PCORnet planned source | Harmonization status | Required audit before freeze |
|---|---|---|---|---|
| Patient ID | `subject_id` | `PATID` | Direct | Uniqueness/linkage |
| Hospital encounter | `hadm_id` | `ENCOUNTERID` | Direct/approximate | Inpatient/EI encounter mapping, duplicate encounters |
| ICU episode | `icustays.intime/outtime` | Local ADT/ICU table | Local extension required | Exact ICU entry/exit timestamps and multiple ICU episodes |
| First broad-spectrum exposure | `PRESCRIPTIONS.starttime` | `PRESCRIBING.RX_START_DATE` plus time field/local order timestamp | Approximate | Start-time completeness, order semantics, route completeness |
| Actual antibiotic administration | Not primary MIMIC exposure | `MED_ADMIN` | PSU enrichment | RxNorm/NDC mapping, route, single-dose versus infusion stop time |
| Specimen time | `microbiologyevents.charttime/chartdate` | `LAB_RESULT_CM.SPECIMEN_DATE/TIME` or local micro | Approximate | Collection-time semantics |
| Result availability | `microbiologyevents.storetime/storedate` | `LAB_RESULT_CM.RESULT_DATE/TIME` or local micro | Approximate | **Critical:** confirm ETL does not substitute specimen date for unavailable result date |
| Organism positivity | `org_name` | Prefer local microbiology organism field | Local extension likely required | Positive/no-growth/flora/contaminant encoding |
| Vasopressors | `inputevents` | `MED_ADMIN` | Approximate | Name/code list, interval construction, route |
| HR/BP/weight/height | `chartevents` | `VITAL` and/or `OBS_CLIN` | Approximate | Time precision, source, repeated readings |
| Temperature/RR/SpO2/GCS/FiO2 | `chartevents` | `OBS_CLIN`/local EHR extracts | Approximate/local | Code/value mappings and completeness |
| Labs | `labevents` | `LAB_RESULT_CM` | Approximate | LOINC/name mapping, specimen time, units |
| Diagnoses | `diagnoses_icd` | `DIAGNOSIS` | Approximate | Lookback timing and diagnosis source |
| Procedures | `procedures_icd` | `PROCEDURES` | Approximate | Timing and mechanical ventilation/source-control mapping |
| Death | `deathtime`, `dod` | `DEATH`, encounter death disposition | Approximate | Out-of-system death completeness |
| Hospital discharge | `admissions.dischtime` | `ENCOUNTER.DISCHARGE_DATE/TIME` | Approximate | Time completeness |

## Medication measurement plan

The primary external replication should use the PSU source that best matches the MIMIC prescription/order construct. `PRESCRIBING` is the default. `MED_ADMIN` should then repeat the exposure classification using delivered medications. Report a reclassification matrix and agreement proportion between order-defined and administration-defined strategies.

## Microbiology warning

PCORnet CDM 6.1 permits `RESULT_DATE` to be populated with `SPECIMEN_DATE` when the true result date is unavailable. Therefore, the presence of `RESULT_DATE/TIME` is not sufficient evidence that it represents clinician-visible result availability. Penn State ETL documentation or a local microbiology source must resolve this before the primary external validation is run.

## Variables that may strengthen PSU sensitivity analyses

If reliably available, add source of infection, source control, immunosuppression, infectious-disease consultation, code status/goals of care, discharge antibiotics, and actual administration data as PSU-only enrichment analyses. Do not insert these silently into the primary harmonized model.

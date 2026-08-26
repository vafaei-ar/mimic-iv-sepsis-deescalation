# PSU / PCORnet harmonization crosswalk

This file records the current frozen harmonization decisions for the Penn State external replication. The PSU analysis is a modified external replication, not an exact reproduction of the MIMIC-IV target trial, because exact ICU timing and faithful culture-positivity linkage are not available in the current PCORnet extract.

| Construct | MIMIC-IV | PSU/PCORnet source | Harmonization status | Frozen decision / limitation |
|---|---|---|---|---|
| Patient ID | `subject_id` | `PATID` | Direct | Use for internal linkage only; never export identifiers. |
| Hospital encounter | `hadm_id` | `ENCOUNTERID` | Direct/approximate | Use encounter linkage; exact ICU episode reconstruction is unavailable from the current extract. |
| ICU episode | `icustays.intime/outtime` | No validated exact local ICU/ADT source identified | Not harmonizable exactly | PSU cannot be described as an exact ICU-clock replication. Any PSU analysis must use a clearly labeled modified hospital-encounter clock unless a new validated ICU source becomes available. |
| First broad-spectrum exposure | `PRESCRIPTIONS.starttime` | `PRESCRIBING` | Approximate, primary PSU source | Use order-based `PRESCRIBING` as primary because it best matches the MIMIC prescription/order construct. Use corrected legacy RxNorm mapping plus broad-spectrum medication-name fallback, with explicit exclusion of recognized oral and inhaled formulations. Because most PSU route values are `UN`, this is a systemic broad-spectrum proxy rather than a fully verified IV-only phenotype. |
| Actual antibiotic administration | Not primary MIMIC exposure | `MED_ADMIN` | PSU sensitivity source | Repeat exposure classification as a sensitivity analysis. Encounter concordance with corrected legacy mapping was very high: all 112,638 MED_ADMIN mapped encounters were contained within 113,657 PRESCRIBING mapped encounters, with 1,019 PRESCRIBING-only encounters. |
| Ambiguous antibiotic formulations | Route/formulation available in source medication records | `PRESCRIBING` and `MED_ADMIN` names/codes | Partially harmonizable | Explicit oral/inhaled signals are excluded. Vancomycin, linezolid, and aztreonam still have many route-unspecified rows, so PSU should not claim verified IV exposure for those records. |
| Specimen time | `microbiologyevents.charttime/chartdate` | `LAB_RESULT_CM.SPECIMEN_DATE/TIME` | Strong timing availability | Culture collection timing is available and can be used descriptively. |
| Result availability | `microbiologyevents.storetime/storedate` | `LAB_RESULT_CM.RESULT_DATE/TIME` | Timing fields present, semantic linkage unresolved | Result timestamps are populated, but the current extract does not preserve a defensible structural link from parent culture rows to organism/positivity results. |
| Organism positivity | `org_name` | Current PSU `LAB_RESULT_CM` extract | Not faithfully recoverable | Multiple linkage strategies failed or produced clear false-positive contamination. Therefore the MIMIC eligibility criterion "no positive clinical culture result available by day 3" cannot be reproduced faithfully in PSU from the current extract. Do not create a pseudo-positive phenotype by further keyword tuning. |
| Vasopressors | `inputevents` | `MED_ADMIN` | Approximate | Requires medication mapping and interval construction if included in PSU covariates. |
| HR/BP/weight/height | `chartevents` | `VITAL` and/or `OBS_CLIN` | Approximate | Audit time precision and source before model construction. |
| Temperature/RR/SpO2/GCS/FiO2 | `chartevents` | `OBS_CLIN`/local EHR extracts | Approximate/local | Code/value mappings and completeness must be validated before use. |
| Labs | `labevents` | `LAB_RESULT_CM` | Approximate | Use validated LOINC/name mapping, timing, and units. |
| Diagnoses | `diagnoses_icd` | `DIAGNOSIS` | Approximate | Harmonize lookback timing and diagnosis source. |
| Procedures | `procedures_icd` | `PROCEDURES` | Approximate | Harmonize timing and mechanical ventilation/source-control mapping. |
| Death | `deathtime`, `dod` | `DEATH`, encounter death disposition | Approximate | PSU death timing is less precise; exact death-through-96h exclusion may not be reproducible. |
| Hospital discharge | `admissions.dischtime` | `ENCOUNTER` | Approximate | Current sepsis encounter extract has date-level discharge information; exact 96h landmark exclusion may require approximation or a more precise source. |

## Frozen medication measurement plan

Primary PSU medication source: `PRESCRIBING`.

Primary PSU broad-spectrum phenotype: corrected legacy RxNorm mapping plus the frozen broad-spectrum medication-name fallback, while excluding recognized oral and inhaled formulations. Do not require the raw PCORnet route field because it is overwhelmingly coded `UN` and would discard most clinically plausible systemic therapy. This phenotype must be described as a **systemic broad-spectrum antibiotic proxy** rather than a verified IV-only phenotype.

Sensitivity analysis: repeat the same phenotype using `MED_ADMIN`. Report encounter-level reclassification/agreement between order-defined and administration-defined strategies.

For the three route-ambiguous agents specifically:

- Vancomycin: retain reviewed include codes `11124`, `202368`, `239209`; explicitly exclude oral product codes `313570`, `313571`, `2000134`; exclude medication-name rows with explicit oral wording.
- Linezolid: include reviewed codes `190376`, `261710`; exclude rows with explicit oral wording. Route-unspecified rows remain a measurement limitation.
- Aztreonam: include reviewed codes `1272`, `202561`; exclude explicit inhaled/Cayston rows. Route-unspecified rows remain a measurement limitation.

## Frozen microbiology decision

The current PSU PCORnet extract reliably identifies culture tests and their collection/result timestamps, but repeated audits did not identify a defensible structural representation of culture positivity linked to those parent culture tests. Exact specimen timestamp, panel-prefix, same-encounter temporal sibling, and refined component-linkage strategies either failed completely or produced obvious unrelated laboratory contamination.

Therefore, do **not** claim that PSU reproduces the MIMIC criterion "no positive clinical culture result available by day 3." Unless a validated local microbiology source becomes available, PSU should either:

1. be presented as a modified external replication that omits this eligibility restriction, with this difference stated prominently; or
2. be used only for secondary/generalizability analyses that do not require exact culture-negative eligibility.

The first option is preferred if the remaining clock and landmark definitions can be made sufficiently transparent.

## Variables that may strengthen PSU sensitivity analyses

If reliably available, add source of infection, source control, immunosuppression, infectious-disease consultation, code status/goals of care, discharge antibiotics, and actual administration data as PSU-only enrichment analyses. Do not insert these silently into the primary harmonized model.

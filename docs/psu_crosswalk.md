# PSU / PCORnet harmonization crosswalk

This file records the frozen harmonization decisions for the Penn State external replication. The PSU analysis is a **modified external replication**, not an exact reproduction of the MIMIC-IV target trial, because exact ICU timing and faithful culture-positivity/result-availability linkage are not available in the current PCORnet extract.

| Construct | MIMIC-IV | PSU/PCORnet source | Harmonization status | Frozen decision / limitation |
|---|---|---|---|---|
| Patient ID | `subject_id` | `PATID` | Direct | Use for internal linkage only; never export identifiers. |
| Hospital encounter | `hadm_id` | `ENCOUNTERID` | Direct/approximate | Use encounter linkage; exact ICU episode reconstruction is unavailable from the current extract. |
| ICU episode | `icustays.intime/outtime` | No validated exact local ICU/ADT source identified | Not harmonizable exactly | PSU uses the hospital encounter clock and cannot be described as an exact ICU-clock replication. |
| First broad-spectrum exposure | `PRESCRIPTIONS.starttime` | `PRESCRIBING` | Approximate, primary PSU source | Use order-based `PRESCRIBING` as primary. Corrected legacy RxNorm mapping plus broad-spectrum medication-name fallback is used, with explicit recognized oral/inhaled/non-systemic exclusions. Because route is usually unspecified, this is a systemic broad-spectrum proxy rather than a verified IV-only phenotype. |
| Actual antibiotic administration | Not primary MIMIC exposure | `MED_ADMIN` | Prespecified PSU sensitivity | Reclassify the 72-96 h exposure window as a measurement sensitivity. Agreement with PRESCRIBING was extremely high in the frozen strict cohort, with only 2 encounters changing treatment classification. |
| Ambiguous antibiotic formulations | Route/formulation available in source medication records | `PRESCRIBING` and `MED_ADMIN` names/codes | Partially harmonizable | Explicit oral/inhaled signals are excluded. Vancomycin, linezolid, and aztreonam have many route-unspecified rows, so PSU must not claim verified IV exposure for those records. |
| Specimen time | `microbiologyevents.charttime/chartdate` | `LAB_RESULT_CM.SPECIMEN_DATE/TIME` | Usable for physiologic/culture collection timing | The audited specimen clock is used where collection timing is needed. |
| Result availability | `microbiologyevents.storetime/storedate` | `LAB_RESULT_CM.RESULT_DATE/TIME` | Not independently usable | In the audited PSU extract, specimen and result timestamps were identical row-by-row. Therefore RESULT_DATE/TIME cannot be treated as an independent result-availability clock. |
| Organism positivity | `org_name` | Current PSU `LAB_RESULT_CM` extract | Not faithfully recoverable | Multiple linkage strategies failed or produced false-positive contamination. The MIMIC criterion `no positive clinical culture result available by day 3` is therefore omitted from the PSU primary cohort rather than approximated with an unvalidated pseudo-phenotype. |
| Vasopressors | `inputevents` | `MED_ADMIN` | Approximate | Raw MED_ADMIN TIME fields were not usable. The frozen date-span rule treats start date as 00:00:00 and stop date as 23:59:59 for interval-overlap construction. |
| HR/RR/SpO2/temperature | `chartevents` | `OBS_CLIN`/validated local extracts | Approximate | Use audited codes, time fields, units, and physiologic-range filters. |
| Mean arterial pressure | `chartevents` | PSU vital/observation components | Derived/approximate | Use the validated derived MAP implementation retained in the frozen PSU covariate set. |
| GCS/FiO2/ventilation/urine output | `chartevents`/`outputevents`/procedure sources | Current PSU extract | Unavailable/unreliable for frozen replication | Omitted from the frozen PSU PS rather than replaced with unvalidated proxies. |
| Labs | `labevents` | `LAB_RESULT_CM` | Approximate | Use validated LOINC/name mapping and specimen date/time. Result timestamps are not used as independent availability times. |
| Diagnoses/comorbidity | `diagnoses_icd` | `DIAGNOSIS` | Approximate | Use the frozen diagnosis composite plus heart-failure and chronic-kidney indicators. |
| Death | `deathtime`, `dod` | `DEATH` | Date-level approximation | Primary mortality is death from the 96-h landmark date through +30 days inclusive. The strict cohort excludes deaths through the landmark using the date-level rule. |
| Hospital discharge | `admissions.dischtime` | sepsis encounter discharge date | Date-level approximation | Strict primary requires discharge date > landmark date; lenient sensitivity allows discharge on the landmark date with >=. |
| Antibiotic burden after landmark | hour-level MIMIC medication coverage | `PRESCRIBING` date intervals | Modified | Measures ordered calendar-day treatment, not verified administration. Missing prescription end date is treated as same-day. |
| Late recurrent/persistent course | in-hospital antibiotic trajectory | `PRESCRIBING` | Exploratory modified outcome | >=3 consecutive systemic-proxy calendar days beginning day 7+ after landmark. Strongly vulnerable to differential discharge/observation. |

## Frozen medication measurement plan

Primary PSU medication source: `PRESCRIBING`.

Primary PSU broad-spectrum phenotype: corrected legacy RxNorm mapping plus the frozen broad-spectrum medication-name fallback, while excluding recognized oral, inhaled, topical/non-systemic formulations. Do not require the raw PCORnet route field because it is overwhelmingly unspecified and would discard most clinically plausible systemic therapy. The phenotype must be described as a **systemic broad-spectrum antibiotic proxy** or **ordered systemic broad-spectrum antibiotic proxy**, not a verified IV-only phenotype.

Prespecified sensitivity: repeat day-3 exposure classification using `MED_ADMIN` and the validated date-span interval rule. This sensitivity changes only exposure measurement; the frozen covariate, PS, and outcome definitions remain unchanged.

For the three route-ambiguous agents specifically:

- Vancomycin: retain reviewed include codes `11124`, `202368`, `239209`; explicitly exclude oral product codes `313570`, `313571`, `2000134`; exclude medication-name rows with explicit oral wording.
- Linezolid: include reviewed codes `190376`, `261710`; exclude rows with explicit oral wording. Route-unspecified rows remain a measurement limitation.
- Aztreonam: include reviewed codes `1272`, `202561`; exclude explicit inhaled/Cayston rows. Route-unspecified rows remain a measurement limitation.

## Frozen microbiology decision

The current PSU extract can identify culture tests and specimen timing, but repeated audits did not identify a defensible structural representation of culture positivity linked to those parent culture tests. In addition, the nominal result timestamp duplicated the specimen timestamp row-by-row and therefore cannot stand in for clinical result availability.

Therefore, do **not** claim that PSU reproduces the MIMIC criterion `no positive clinical culture result available by day 3`, and do not create a pseudo-positive phenotype by further keyword tuning. The frozen PSU analysis omits this eligibility restriction and is explicitly labeled a modified external replication.

## Frozen primary covariate set

Retained constructs are age, sex, race, diagnosis/comorbidity indicators, heart failure, chronic kidney disease, lactate, creatinine, white blood cells, platelets, total bilirubin, heart rate, respiratory rate, SpO2, temperature, derived MAP, early vasopressor exposure, late vasopressor exposure, and vasopressor stopped before 72 h, including the frozen early/late/delta forms where specified in code.

Direct GCS, FiO2, ventilation procedure, urine output, exact admission-to-ICU hours, ICU type, microbiology intensity, broad-antibiotic hours before 72 h, and BMI are unavailable or insufficiently reliable in the current extract and are omitted from the frozen primary PSU PS.

The primary PSU analysis uses the MIMIC-harmonized missing-data convention: binary fill with 0 and continuous median imputation/standardization, without adding missingness indicators. Late laboratory missingness differs between exposure groups and remains a limitation; it was not used to tune the primary model after treatment effects were observed.

## Publication freeze rule

No new PSU-only enrichment variables, subgroup analyses, alternative propensity models, new trimming thresholds, or revised outcome definitions should be added to the publication analysis after seeing the treatment effects unless explicitly labeled as a new post hoc analysis. The accepted publication sensitivity set is limited to the prespecified weighting checks, MED_ADMIN exposure reclassification, and lenient 96-h landmark analysis documented in the reviewer walkthrough.

# Frozen target-trial / landmark specification

This document is the common scientific contract for MIMIC-IV and the PSU external replication. Site-specific code may differ, but it must not silently change these constructs.

## Population and time axis

- Adult patients with a qualifying ICU stay.
- Trial clock `t0`: first qualifying systemic IV broad-spectrum antibiotic exposure.
- Qualifying first exposure: from 6 h before ICU entry through 24 h after ICU entry.
- Clinical microbiology specimen collected from 24 h before `t0` through the 72-h decision time.
- Primary culture-negative rule: no qualifying positive organism result **available by the 72-h decision time**.
- Alive and hospitalized through the 96-h landmark.
- No active vasopressor overlap during the 6 h before the 72-h decision.
- Broad-spectrum therapy still present during 48-72 h.

## Treatment strategies

Treatment is observed during 72-96 h after `t0`.

- `A=1`: no systemic IV broad-spectrum coverage during the 72-96-h classification window.
  - `narrowed_or_non_broad_only`: systemic antibiotic therapy remains but no broad-spectrum therapy remains.
  - `stopped_all_observed_systemic_antibiotics`: no observed systemic antibiotic coverage.
- `A=0`: any systemic IV broad-spectrum coverage overlaps the 72-96-h window.

This is an operational treatment phenotype. It should not be described as a clinician's intention when only order/prescription data are available.

## Follow-up

- Follow-up begins at 96 h after first broad-spectrum exposure.
- Primary safety outcome: 30-day all-cause mortality from the landmark.
- Patient-centered secondary outcome: hospital-free days through 30 days, death assigned 0.
- Stewardship/treatment-separation outcomes: antibiotic-free days, normalized systemic antibiotic exposure, normalized broad-spectrum exposure.
- Late recurrent/persistent antibiotic-course use is exploratory because inpatient observation time differs with discharge.

## Confounding adjustment

Only information available before the 72-h decision may enter the primary propensity model. The final model includes demographics, comorbidity, baseline and near-decision physiology, organ-dysfunction trajectories, vasopressor trajectories, diagnostic intensity, antibiotic intensity, steroids, and BMI proxy when available.

The progressive M1-M4 sequence is retained because attenuation after clinical-improvement adjustment is a central scientific result.

## External replication rules

1. PSU uses the same conceptual `t0`, 72-h decision, 72-96-h classification, and 96-h landmark.
2. PSU does not use hospital admission as `t0`.
3. A coded/pre-filtered sepsis cohort is not required for the primary replication unless the MIMIC population is changed to require it.
4. Differences in data representation must be documented in `docs/psu_crosswalk.md`.
5. MIMIC prescription/order exposure and PSU PRESCRIBING should form the closest primary measurement match when feasible.
6. PSU MED_ADMIN is analyzed separately as an administration-based measurement-validation sensitivity.
7. PSU microbiology result timestamps must be validated against local ETL semantics before treating RESULT_DATE/TIME as result availability.

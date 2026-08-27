# Frozen target-trial / landmark specification

This document defines the scientific contract for the MIMIC-IV primary analysis and identifies which parts can and cannot be reproduced exactly in the Penn State (PSU) modified external replication.

## MIMIC-IV population and time axis

- Adult patients with a qualifying ICU stay.
- Trial clock `t0`: first qualifying systemic IV broad-spectrum antibiotic exposure.
- Qualifying first exposure: from 6 h before ICU entry through 24 h after ICU entry.
- Clinical microbiology specimen collected from 24 h before `t0` through the 72-h decision time.
- Primary culture-negative rule: no qualifying positive clinical culture result **available by the 72-h decision time**.
- Alive and hospitalized through the 96-h landmark.
- No active vasopressor overlap during the 6 h before the 72-h decision.
- Broad-spectrum therapy still present during 48-72 h.

## MIMIC-IV treatment strategies

Treatment is observed during 72-96 h after `t0`.

- `A=1`: no systemic IV broad-spectrum coverage during the 72-96-h classification window.
  - `narrowed_or_non_broad_only`: systemic antibiotic therapy remains but no broad-spectrum therapy remains.
  - `stopped_all_observed_systemic_antibiotics`: no observed systemic antibiotic coverage.
- `A=0`: any systemic IV broad-spectrum coverage overlaps the 72-96-h window.

This is an operational treatment phenotype based on structured medication records. It should not be described as a clinician's intention, and prescription/order data should not be described as confirmed bedside administration.

## Follow-up and outcomes

- Follow-up begins at 96 h after first broad-spectrum exposure.
- Primary safety outcome: 30-day all-cause mortality from the landmark.
- Hospital-free days through 30 days: death assigned 0.
- Antibiotic-free days through 30 days: death assigned 0.
- Normalized systemic antibiotic exposure and normalized broad-spectrum exposure are treatment-burden outcomes.
- Late recurrent/persistent antibiotic-course use is exploratory because inpatient observation time differs with discharge.

## Confounding adjustment

Only information available before the 72-h decision may enter the primary propensity model. The final corrected MIMIC model includes demographics, comorbidity, baseline and near-decision physiology, organ-dysfunction/clinical trajectories, vasopressor trajectories, diagnostic intensity, antibiotic intensity, steroids, and BMI proxy where validly measured.

Direct GCS and FiO2 covariates were removed after a targeted measurement audit showed that their saved implementations were non-informative/unreliable. This was an audit-driven correction, not a prespecified modeling change, and is documented in `docs/mimic_v57_freeze_review.md`.

The progressive M1-M4 sequence is retained because attenuation after clinical-improvement adjustment is a central scientific result.

## PSU modified external replication

PSU preserves the **conceptual** treatment clock and estimand structure but cannot reproduce all MIMIC eligibility/measurement details. These differences are intentional and must remain explicit.

1. PSU anchors on the first qualifying broad-spectrum `PRESCRIBING` order in the first 24 h of the hospital encounter because a validated exact ICU `intime` clock is unavailable in the current extract.
2. PSU preserves decision time = anchor +72 h, exposure classification over 72-96 h, and follow-up beginning at the 96-h landmark.
3. PSU does **not** impose the MIMIC criterion `no positive clinical culture result available by 72 h`. Culture specimen timing can be identified, but the audited result timestamp is not an independent result-availability clock and faithful positivity linkage is unavailable.
4. PSU primary medication exposure is an **ordered systemic broad-spectrum antibiotic proxy** from `PRESCRIBING`, not verified IV administration. Explicit recognized oral/inhaled formulations are excluded.
5. PSU `MED_ADMIN` is a prespecified exposure-measurement sensitivity. Its raw time fields were not usable, so validated date-span intervals are used for this sensitivity.
6. PSU discharge, death, and medication intervals include date-level semantics. The strict date-level 96-h landmark rule is primary; the lenient calendar-date rule is prespecified sensitivity.
7. The frozen PSU covariate model uses only constructs that could be mapped defensibly from the available extract. Missing/unreliable MIMIC constructs are omitted rather than replaced with unvalidated proxies.
8. All site-specific differences and reasons are documented in `docs/psu_crosswalk.md` and `docs/psu_analysis_walkthrough.md`.

Therefore PSU should be described as a **modified external replication**, not as an exact validation of the MIMIC culture-negative ICU target trial.

# PSU modified external replication: reviewer walkthrough

This document is the code-review map for the Penn State (PSU) analysis. It is written for a collaborator who was not involved in the original implementation and needs to understand both **what the code does** and **why each non-obvious decision was made**.

## 1. What this analysis is, and is not

The PSU analysis is a **modified external replication** of the MIMIC-IV target-trial emulation. It is not an exact replication.

The conceptual treatment clock is preserved:

- anchor: first qualifying broad-spectrum antibiotic order in the first 24 hours of the hospital encounter;
- decision time: anchor + 72 hours;
- treatment-classification window: 72-96 hours;
- follow-up start: the 96-hour landmark;
- primary outcome: 30-day all-cause mortality from the landmark.

However, the available PSU extract differs from MIMIC in four important ways:

1. **No validated day-3 culture-result-availability phenotype.** The PSU lab extract identifies culture tests and specimen timing, but the audited result timestamp is not an independent result-availability clock. The PSU primary cohort therefore does not impose MIMIC's `no positive clinical culture result available by 72 h` restriction.
2. **No validated exact ICU clock.** PSU uses the hospital encounter clock. Therefore the paper must not imply that the PSU cohort is an exact ICU-timed replication.
3. **PRESCRIBING is ordered therapy, not verified administration.** The primary exposure is a systemic broad-spectrum antibiotic proxy. It should not be described as verified IV administration.
4. **Several event times are date-level.** Discharge, death, prescription intervals, and some administration intervals require calendar-date approximations. The strict landmark rule is primary; a lenient landmark rule is prespecified sensitivity analysis.

These differences are scientific limitations, not software defects. They are intentionally visible in the code and manuscript rather than hidden behind forced harmonization.

## 2. Data safety

The PSU source data remain on the approved Penn State machine. Scripts may use patient/encounter identifiers internally for joins, but only aggregate outputs are declared as RunRelay artifacts or copied outside the local analysis environment.

Do not add patient-level CSV/Parquet outputs to this repository or to RunRelay artifact lists.

## 3. Analysis stages and entry points

The PSU work was frozen in stages to prevent outcome-driven tuning.

### Stage A. Source and semantic audits

These scripts established what the PSU fields actually mean before treatment effects were estimated:

- `audit_psu_target_trial_crosswalk.py`
- `audit_psu_microbiology_semantics.py` and related culture-linkage audits
- `audit_psu_antibiotic_source_mapping.py`
- `audit_psu_antibiotic_route_rxnorm.py`
- `audit_psu_antibiotic_rxnorm_specificity.py`
- `audit_psu_lab_clock_semantics.py`
- `audit_psu_medadmin_clock_semantics.py`

The purpose of these scripts is semantic validation, not effect estimation.

### Stage B. Frozen cohort/covariates and propensity score

Primary code:

- `audit_psu_final_covariate_freeze.py`
- `audit_psu_ps_balance.py`

`audit_psu_final_covariate_freeze.py` re-exports the frozen base definitions and applies the validated MED_ADMIN date-span timing rule for vasopressor covariates. The raw MED_ADMIN TIME fields were audited and found unusable; start dates are therefore interpreted at 00:00:00 and stop dates at 23:59:59.

`audit_psu_ps_balance.py` constructs the strict cohort, primary PRESCRIBING exposure, frozen covariates, propensity model, overlap/weight diagnostics, missingness diagnostics, and SMDs. It intentionally does **not** estimate treatment effects.

Frozen strict primary cohort parity target:

- cohort: 19,841 encounters;
- de-escalated: 5,346;
- continued broad-spectrum: 14,495.

Frozen primary PS diagnostics:

- maximum post-weighting absolute SMD: approximately 0.02262;
- treated ESS: approximately 4,513;
- continued ESS: approximately 14,222;
- maximum stabilized weight: approximately 7.20.

If a software cleanup changes these values materially, treat that as a parity failure until explained.

### Stage C. Frozen outcomes

Primary code:

- `audit_psu_final_outcome_freeze.py`

This stage was deliberately run **without exposure groups or treatment-effect estimation**. It froze outcome definitions before the treatment-effect run.

Important outcome decisions:

- `death_30d`: death from the 96-hour landmark date through +30 days, inclusive;
- hospital-free days: death-to-zero, otherwise 30 minus landmark-to-index-discharge calendar days;
- antibiotic-free days: death-to-zero, otherwise 30 minus calendar days with the frozen systemic-antibiotic proxy;
- normalized antibiotic exposure: antibiotic calendar days divided by days alive through 30 days;
- late recurrent/persistent antibiotic course: at least 3 consecutive systemic-antibiotic calendar days beginning day 7 or later after the landmark.

The late-course outcome is exploratory because discharge shortens observation and can mechanically reduce the opportunity to observe later treatment.

Frozen overall outcome parity target includes 2,381 deaths within the post-landmark 30-day horizon.

### Stage D. Point estimates

Primary code:

- `run_psu_point_estimates.py`

Primary estimand:

- stabilized inverse-probability-treatment weighted ATE;
- effect direction is always **de-escalation minus continuation** for differences and **de-escalation / continuation** for risk ratios.

Prespecified weighting sensitivities:

- overlap weighting, which targets the overlap population (ATO) rather than the ATE;
- stabilized IPTW truncated at 1st/99th percentiles;
- stabilized IPTW truncated at 2.5th/97.5th percentiles.

Do not substitute the overlap estimate for the primary ATE after seeing results simply because its balance is better.

### Stage E. Bootstrap inference

Primary code:

- `run_psu_bootstrap_inference.py`

The bootstrap uses 1,000 encounter-level nonparametric resamples with a fixed seed. The propensity score is **refit within every bootstrap sample**. This is intentional: a fixed-weight bootstrap would condition on the originally estimated PS and understate one component of uncertainty.

The same binary/continuous preprocessing and the same GLM fallback are applied inside each replicate.

### Stage F. Prespecified robustness analyses

Primary code:

- `run_psu_prespecified_robustness.py`
- `run_psu_prespecified_robustness_bootstrap.py`

Only two PSU robustness variants were prespecified before effect inspection:

1. **MED_ADMIN exposure sensitivity.** The strict primary cohort is retained, but day-3 broad-spectrum exposure is reclassified using MED_ADMIN with the validated date-span timing rule.
2. **Lenient 96-hour landmark sensitivity.** Because PSU discharge/death are date-level, patients with discharge/death on the landmark calendar date are allowed by `>=` instead of excluded by the strict `>` rule.

These variants change only the intended exposure/landmark component. Covariates, PS specification, outcomes, and weighting remain frozen.

## 4. Why PRESCRIBING is the primary PSU exposure

MIMIC's primary exposure is based on structured prescribing/order information and inferred coverage intervals. PSU PRESCRIBING is therefore the closest available source representation.

MED_ADMIN is not automatically 'better' for harmonization because its timing semantics differ and its raw TIME fields were not usable. It is used as a measurement sensitivity rather than silently replacing the primary exposure.

Because route data are incomplete, the PSU manuscript language must be **systemic broad-spectrum antibiotic proxy** or **ordered broad-spectrum therapy**, not verified IV treatment.

## 5. Antibiotic mapping decisions

The PSU mapping combines:

- the audited legacy RxNorm broad-spectrum code list;
- broad-spectrum medication-name matching;
- explicit exclusions from the legacy mapping;
- removal of clearly non-systemic names/routes (for example oral vancomycin and inhaled aztreonam where detectable).

This mapping was frozen only after route/RxNorm contamination audits. Do not broaden the pattern based on whether an added drug makes the treatment effect larger or smaller.

## 6. Covariate decisions

The final PSU PS uses the subset of MIMIC concepts that could be defensibly mapped in the local data. Retained families include:

- demographics;
- diagnosis/comorbidity proxies;
- early and near-decision laboratory trajectories;
- early and near-decision vital-sign trajectories;
- vasopressor windows and stopping before the decision.

Unavailable or rejected concepts include GCS, FiO2, ventilation procedure, urine output, exact ICU timing/type, microbiology intensity, broad-antibiotic hours before 72 h, and BMI.

These were excluded because the required source semantics were unavailable or not validated, not because of their relationship to the outcome.

Primary missing-data handling intentionally mirrors the frozen harmonized approach:

- binary covariates: missing -> 0;
- continuous covariates: median imputation;
- no missingness indicators in the primary model;
- continuous variables standardized and clipped to +/-8 before PS fitting.

The missingness-by-exposure report is retained because some late laboratory measurements differ in availability between treatment groups. That is a limitation to report, not a reason to tune the primary model after seeing the effect.

## 7. Strict versus lenient landmark

The strict cohort requires discharge and death dates to be strictly later than the calendar date corresponding to anchor + 96 h. This avoids knowingly including people whose event may have occurred before the exact 96-hour landmark, at the cost of excluding some events that may actually have occurred later on the same calendar day.

The lenient sensitivity uses `>=` and therefore brackets this date-resolution ambiguity from the other side.

## 8. Why some scripts use source injection

`run_psu_point_estimates.py` and the bootstrap/robustness wrappers insert small, named blocks into the already-audited `audit_psu_ps_balance.py` source at literal markers.

This is unusual. It was chosen as a **parity-preserving publication-freeze strategy**: the frozen PS construction is executed verbatim instead of being copied into several independent scripts. Each marker must appear exactly once; otherwise execution fails.

For long-term maintenance, extracting the frozen pipeline into reusable functions would be cleaner. Do that only on a dedicated refactor branch and prove numerical parity against the frozen cohort, PS balance, ESS, weights, point estimates, and bootstrap intervals before merging.

## 9. Reproducing the frozen PSU analyses locally

Create the environment described in the root README, then define the local data root without committing it:

```bash
export PSU_DATA_ROOT=/path/to/approved/Hwang_Bonavia/data/root
```

Run the frozen stages in order:

```bash
python scripts/audit_psu_final_covariate_freeze.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_final_covariate_freeze/latest

python scripts/audit_psu_ps_balance.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_ps_balance/latest

python scripts/audit_psu_final_outcome_freeze.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_final_outcome_freeze/latest

python scripts/run_psu_point_estimates.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_point_estimates/latest

python scripts/run_psu_bootstrap_inference.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_bootstrap_inference/latest

python scripts/run_psu_prespecified_robustness.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_prespecified_robustness/latest

python scripts/run_psu_prespecified_robustness_bootstrap.py "$PSU_DATA_ROOT" \
  --output-dir outputs/psu_prespecified_robustness_bootstrap/latest
```

The publication RunRelay jobs used the same entry points on the approved machine, with exact git commits and manual approval.

## 10. Reviewer checklist

Before accepting a code change, verify:

1. no patient-level outputs were added;
2. strict cohort remains 19,841 unless the change intentionally targets cohort construction;
3. primary treatment counts remain 5,346 vs 14,495;
4. primary post-weighting max absolute SMD remains approximately 0.02262;
5. treated/continued ESS and maximum weight reproduce within numerical tolerance;
6. frozen overall death count remains 2,381;
7. point-estimate direction is A=1 de-escalation minus A=0 continuation;
8. bootstrap refits the PS within replicate;
9. MED_ADMIN and lenient-landmark code change only their prespecified component;
10. manuscript terminology still describes PSU as a modified external replication and the exposure as an ordered/systemic proxy.

## 11. What not to change without a new scientific amendment

Do not add post hoc subgroup analyses, new PS covariates, new trimming thresholds, alternative outcome definitions, machine-learning PS models, or new exposure windows merely because the current mortality estimates are surprising or differ between MIMIC and PSU.

Software refactoring is welcome. Scientific-definition changes require explicit documentation as a new analysis rather than being hidden inside cleanup work.

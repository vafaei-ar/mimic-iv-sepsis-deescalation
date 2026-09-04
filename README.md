# MIMIC-IV Sepsis Day-3 De-escalation

Reproducible analysis code for the MIMIC-IV day-3 broad-spectrum antibiotic de-escalation target-trial emulation and the Penn State (PSU) **modified external replication**.

## Start here

For a scientific/code review, read these files in order:

1. `docs/target_trial_spec.md` — frozen target-trial contract.
2. `docs/mimic_analysis_walkthrough.md` — reviewer-oriented MIMIC workflow, rationale, corrected-vital publication path, and parity targets.
3. `docs/psu_crosswalk.md` — MIMIC-to-PSU data-source crosswalk.
4. `docs/psu_analysis_walkthrough.md` — reviewer-oriented PSU pipeline, rationale, parity targets, and exact reproduction commands.
5. `docs/mimic_v57_freeze_review.md` — historical MIMIC publication-freeze audit and implementation notes.

The comments in the final PSU entry-point scripts intentionally explain non-obvious scientific decisions next to the code that implements them. Historical audit scripts are retained because they document how source semantics were validated before the final definitions were frozen.

## Scientific design

### MIMIC-IV primary analysis

The primary MIMIC analysis is anchored at the first qualifying systemic IV broad-spectrum antibiotic exposure. The treatment decision time is 72 h later, treatment is classified over 72-96 h, and follow-up starts at the 96-h landmark. The primary culture-negative definition requires qualifying microbiology sampling and no positive clinical culture result available by the 72-h decision time.

The primary estimand is the ATE estimated using stabilized inverse probability treatment weighting. Covariates are measured before the 72-h treatment decision. The final model includes clinical status and treatment/intensity trajectories because confounding by clinical improvement is central to the scientific question.

### PSU modified external replication

PSU preserves the conceptual first-broad-spectrum anchor, 72-h decision, 72-96-h classification window, and 96-h landmark, but it is **not an exact MIMIC replication**. The available PSU data use different source semantics and granularity for several constructs, including hospital/ICU representation, medication exposure, microbiology, and date-level event timing.

The primary PSU exposure therefore uses the closest defensible structured source: `PRESCRIBING`, interpreted as an **ordered systemic broad-spectrum antibiotic proxy**, not verified IV administration. `MED_ADMIN` is retained as a prespecified measurement sensitivity. These differences are documented explicitly in `docs/psu_analysis_walkthrough.md` and `docs/psu_crosswalk.md`.

The frozen PSU analytic cohort is drawn from the upstream `sepsis_encounter` source. The manuscript describes that source cohort as meeting an adapted sepsis definition supplied by the Penn State team: suspected or confirmed infection plus an absolute modified SOFA score >=2, with the neurologic/GCS component omitted because GCS could not be reliably mapped in PCORI. This is not identical to the full Sepsis-3 >=2-point change definition.

## Repository layout

```text
config/                 analysis configuration and site mapping templates
src/sepsis_deescalation reusable MIMIC analysis package
scripts/                audit and final command-line entry points
tests/                  unit/smoke and publication-contract tests
docs/                   scientific contracts, freeze reviews, and reviewer guides
outputs/                 generated local results/caches; ignored by git
.runrelay/               approved-machine execution manifest
```

## Data policy

**Never commit patient-level MIMIC-IV, PSU, PCORnet, or derived analytic data.** Patient-level CSV/Parquet files and `outputs/` are ignored by git. The repository should contain code, configuration templates, tests, documentation, and non-sensitive aggregate outputs only.

The public repository does not contain credentialed MIMIC-IV data or restricted PSU data. Reviewers who do not have data access can still inspect all phenotype, statistical, and documentation logic.

## Installation

Use a repository-local virtual environment:

```bash
/usr/bin/python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev]'
pytest -q
ruff check .
```

Python dependencies are constrained by major version in `pyproject.toml`. Publication runs also record the exact git commit and execution metadata through RunRelay.

## MIMIC-IV reproduction

Set the local MIMIC-IV 3.1 root without committing it:

```bash
export MIMIC_SOURCE=/path/to/mimiciv/3.1
```

Validate first:

```bash
python scripts/validate_mimic.py --config config/mimic.yaml
```

### Fast development mode

```bash
python scripts/run_mimic.py \
  --config config/mimic.yaml \
  --mode fast \
  --jobs auto
```

### Final publication mode

The final manuscript source is the corrected-vital inference rerun, not the pre-repair point estimates from the source-dependent base run.

```bash
python scripts/run_mimic.py \
  --config config/mimic.yaml \
  --mode final \
  --jobs auto

export RUN_DIR=outputs/mimic/mimic_iv_v5_7_final_YYYYMMDDTHHMMSSZ

python scripts/repair_v57_vital_covariates.py "$RUN_DIR" \
  --config config/mimic.yaml

python scripts/rerun_inference.py "$RUN_DIR" \
  --config config/mimic.yaml \
  --mode final \
  --jobs auto \
  --cohort-path "$RUN_DIR/audits/vital_repair/analysis_cohort_vital_corrected.csv" \
  --label vital_corrected_final
```

The manuscript primary/secondary outcomes, progressive-adjustment sequence, and final weighting diagnostics come from the corrected `final_vital_corrected_final_*` inference rerun. Source-dependent microbiology and missing-stop-time sensitivities still require the complete base run.

`--jobs auto` uses up to eight worker processes. BLAS/OpenMP libraries are restricted to one thread inside each bootstrap worker to avoid CPU oversubscription. This affects runtime, not the estimand.

The optimized bootstrap engine uses a numeric design matrix rather than rebuilding formulas for every replicate. The primary/secondary outcomes share one propensity-score fit per bootstrap replicate. Point-estimate definitions and target estimands are unchanged.

### Inference-only resume mode

A complete run can write a local patient-level checkpoint under `outputs/cache/mimic/<analysis_version>/`. This checkpoint never belongs in git and remains subject to MIMIC data-use restrictions.

Inference-only mode does not replace a complete source-dependent final run for microbiology and other source-level sensitivities.

## PSU reproduction

The final PSU publication analysis is not produced by the generic `run_pcornet.py` template. It is produced by the frozen PSU audit/inference sequence below. This distinction matters because the PSU extract required several source-semantic audits before a defensible modified replication could be frozen.

Set the approved local data root:

```bash
export PSU_DATA_ROOT=/path/to/approved/psu/data/root
```

Run the stages in order:

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

Expected frozen parity targets and the rationale for each stage are listed in `docs/psu_analysis_walkthrough.md`.

## Publication figures

The preferred manuscript-facing figure set is generated by:

```bash
python scripts/build_nature_figures.py
```

This writes the combined Figure 1, progressive-adjustment Figure 2, cross-dataset Figure 3, and ESM diagnostic figures under `outputs/publication_integration/nature_figures/`. The figure code reads frozen manuscript-facing aggregate values where available; raw/row-level data remain local and are never committed.

## Reproducibility rules

1. Work from a named git commit and a clean repository state.
2. Use a repository-local environment and do not rely on ad hoc global packages.
3. Do not edit generated effect tables by hand.
4. Keep bootstrap seeds and publication replicate counts fixed in code/configuration.
5. Freeze scientific definitions before looking at treatment-effect changes.
6. Site-specific differences must be documented rather than silently forced into MIMIC semantics.
7. Patient-level checkpoints and source data stay local and out of git/artifact transport.
8. A software refactor of frozen PSU code must demonstrate parity on cohort counts, PS balance, ESS/weights, outcomes, point estimates, and bootstrap intervals before replacing the publication implementation.
9. Overlap weighting is a different estimand (ATO) and must not silently replace the primary ATE.
10. PSU should be described as a **modified external replication**, and its primary medication exposure as an ordered/systemic proxy rather than verified IV administration.

## Why historical audit scripts remain in the repository

The many `audit_psu_*` scripts are intentional provenance. They record how ICU timing, microbiology, antibiotic mapping, route coding, laboratory clocks, MED_ADMIN timing, covariate availability, missingness, and outcome observability were evaluated before the final definitions were selected.

They are not all required for every reproduction run, but removing them would erase the evidence for several non-obvious data decisions. The reviewer walkthrough identifies which scripts are historical audits and which are final analysis entry points.

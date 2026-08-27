# MIMIC-IV Sepsis Day-3 De-escalation

Reproducible analysis code for the MIMIC-IV day-3 broad-spectrum antibiotic de-escalation target-trial emulation and the Penn State (PSU) **modified external replication**.

## Start here

For a scientific/code review, read these files in order:

1. `docs/target_trial_spec.md` — frozen target-trial contract and explicit PSU deviations.
2. `docs/mimic_analysis_walkthrough.md` — reviewer-oriented MIMIC pipeline, final vital-sign correction provenance, and exact publication reproduction sequence.
3. `docs/mimic_v57_freeze_review.md` — final corrected MIMIC results, diagnostics, and reporting rules.
4. `docs/psu_crosswalk.md` — MIMIC-to-PSU data-source crosswalk and frozen limitations.
5. `docs/psu_analysis_walkthrough.md` — reviewer-oriented PSU pipeline, rationale, parity targets, and exact reproduction commands.

Comments in the maintained publication entry points intentionally explain non-obvious scientific decisions next to the code that implements them. Historical audit scripts are retained because they document how source semantics were validated before final definitions were frozen.

## Scientific design

### MIMIC-IV primary analysis

The primary MIMIC analysis is anchored at the first qualifying systemic IV broad-spectrum antibiotic exposure. The treatment decision time is 72 h later, treatment is classified over 72-96 h, and follow-up starts at the 96-h landmark. The primary culture-negative definition requires qualifying microbiology sampling and no positive clinical culture result available by the 72-h decision time.

The primary estimand is the ATE estimated using stabilized inverse probability treatment weighting. Covariates are measured before the 72-h treatment decision. The final model includes clinical status and treatment/intensity trajectories because confounding by clinical improvement is central to the scientific question.

The final manuscript numbers come from the **audited vital-corrected inference rerun**, not directly from the pre-repair point estimates produced by the initial source-dependent v5.7 run. The exact three-stage sequence is documented below and in `docs/mimic_analysis_walkthrough.md`.

### PSU modified external replication

PSU preserves the conceptual 72-h decision, 72-96-h classification window, and 96-h landmark, but it is **not an exact MIMIC replication**. The available PSU extract does not provide a validated day-3 culture-result-availability phenotype or an exact ICU clock, medication route is incomplete, and several timestamps are date-level.

The primary PSU exposure therefore uses the closest defensible structured source: `PRESCRIBING`, interpreted as an **ordered systemic broad-spectrum antibiotic proxy**, not verified IV administration. `MED_ADMIN` is retained as a prespecified measurement sensitivity. These differences are documented explicitly in `docs/psu_analysis_walkthrough.md` and `docs/psu_crosswalk.md`.

## Repository layout

```text
config/                 analysis configuration and site mapping templates
src/sepsis_deescalation reusable MIMIC analysis package
scripts/                audit, correction, and final command-line entry points
tests/                  unit/smoke/publication-contract tests
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
ruff check src scripts tests
```

Python dependencies are constrained in `pyproject.toml`. Publication execution is tied to a named git commit, and approved-machine runs record their execution provenance through RunRelay.

## MIMIC-IV publication reproduction

Set the local MIMIC-IV 3.1 root without committing it:

```bash
export MIMIC_SOURCE=/path/to/mimiciv/3.1
```

Validate first:

```bash
python scripts/validate_mimic.py --config config/mimic.yaml
```

### 1. Source-dependent v5.7 run

Run the complete source-dependent analysis. This establishes cohort membership, exposure classification, microbiology eligibility, outcomes, and source-dependent sensitivities.

```bash
python scripts/run_mimic.py \
  --config config/mimic.yaml \
  --mode final \
  --jobs auto
```

Save the timestamped run directory printed by the command:

```bash
export RUN_DIR=outputs/mimic/mimic_iv_v5_7_final_YYYYMMDDTHHMMSSZ
```

### 2. Apply the audited vital-sign correction

The final v5.7 audit found that baseline temperature required reading-level Fahrenheit/Celsius normalization and that routine GCS/FiO2 measurement needed targeted reconstruction/audit. The accepted correction is kept as an explicit provenance stage rather than being hidden inside the historical base run.

```bash
python scripts/repair_v57_vital_covariates.py "$RUN_DIR" \
  --config config/mimic.yaml
```

This writes the corrected patient-level analytic cohort locally under:

```text
$RUN_DIR/audits/vital_repair/analysis_cohort_vital_corrected.csv
```

That file is restricted analysis data and must not be committed or transported as a public artifact.

### 3. Rerun final inference from the corrected cohort

```bash
python scripts/rerun_inference.py "$RUN_DIR" \
  --config config/mimic.yaml \
  --mode final \
  --jobs auto \
  --cohort-path "$RUN_DIR/audits/vital_repair/analysis_cohort_vital_corrected.csv" \
  --label vital_corrected_final
```

**The manuscript primary/secondary, progressive-adjustment, and final weighting results come from this corrected `final_vital_corrected_final_*` inference rerun.** Do not substitute the pre-repair effect estimates in the base run directory.

The base source-dependent run remains necessary for microbiology-membership and missing-stop-time sensitivities that cannot be reconstructed from an inference-only checkpoint.

### Fast development mode

For software development, the initial source-dependent run may be exercised with reduced bootstrap counts:

```bash
python scripts/run_mimic.py \
  --config config/mimic.yaml \
  --mode fast \
  --jobs auto
```

`--jobs auto` uses up to eight worker processes. BLAS/OpenMP libraries are restricted to one thread inside each bootstrap worker to avoid CPU oversubscription. This changes runtime only, not the estimand. Fast mode is not publication-quality inference.

The optimized bootstrap engine uses a numeric design matrix rather than rebuilding formulas for every replicate. The primary/secondary outcomes share one propensity-score fit per bootstrap replicate, and the PS is refit within each bootstrap sample. See `src/sepsis_deescalation/fast_bootstrap.py` for the documented runtime-parity guardrails.

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

## Reproducibility rules

1. Work from a named git commit and a clean repository state.
2. Use a repository-local environment and do not rely on ad hoc global packages.
3. Do not edit generated effect tables by hand.
4. Keep bootstrap seeds and publication replicate counts fixed in code/configuration.
5. Freeze scientific definitions before looking at treatment-effect changes.
6. Site-specific differences must be documented rather than silently forced into MIMIC semantics.
7. Patient-level checkpoints and source data stay local and out of git/artifact transport.
8. A software refactor of frozen publication code must demonstrate parity on the relevant cohort, feature, PS, outcome, point-estimate, and bootstrap outputs before replacing the accepted implementation.
9. Overlap weighting is a different estimand (ATO) and must not silently replace the primary ATE.
10. PSU should be described as a **modified external replication**, and its primary medication exposure as an ordered/systemic proxy rather than verified IV administration.

## Why historical audit scripts remain in the repository

The `audit_mimic_*` and `audit_psu_*` scripts are intentional provenance. They record how vital-unit semantics, ICU timing, microbiology, antibiotic mapping, route coding, laboratory clocks, MED_ADMIN timing, covariate availability, missingness, and outcome observability were evaluated before final definitions were selected.

They are not all required for every reproduction run, but removing them would erase the evidence for several non-obvious data decisions. The reviewer walkthroughs identify which scripts are historical audits and which are maintained publication entry points.

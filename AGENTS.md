# mimic-iv-sepsis-deescalation agent instructions

## Project scope

Treat `vafaei-ar/mimic-iv-sepsis-deescalation` as the only active scientific project when working from this repository unless the user explicitly expands scope.

Do not inspect, infer from, or mix in other user repositories to supply project-specific scientific assumptions or files. The private RunRelay control repository `vafaei-ar/RunRelay` may be accessed only for RunRelay job orchestration.

## RunRelay execution

This repository uses RunRelay for workstation execution.

- Project id: `mimic-iv-sepsis-deescalation`
- RunRelay control repository: `vafaei-ar/RunRelay`
- Execution manifest: `.runrelay/project.yaml`
- Bound RunRelay machine: `pshjl4vf24`
- Local checkout: `/home/asadr/works/repos/mimic-iv-sepsis-deescalation`
- Treat `.runrelay/project.yaml` as the authoritative list of named tasks and execution policy.
- For every RunRelay job for this repository, set `requested_machine_id` to `pshjl4vf24`. Never leave it null or infer a machine from another project.
- This is a public repository, so use RunRelay safe mode, exact commits, named tasks only, and Telegram manual approval through the private RunRelay control repository.
- Prefer RunRelay over asking the user to execute shell commands manually when an equivalent named task exists.

When a new execution operation is needed, add a narrowly scoped named task to `.runrelay/project.yaml`; do not introduce arbitrary shell execution.

## Standard execution flow

1. Modify code/configuration in this repository only.
2. Commit the intended state and obtain the exact full commit SHA.
3. Re-read committed `AGENTS.md` and `.runrelay/project.yaml` before creating a RunRelay job.
4. Create a unique immutable job JSON under `jobs/` in `vafaei-ar/RunRelay` with project id `mimic-iv-sepsis-deescalation`, the exact project commit, an allowed named task, and `requested_machine_id: "pshjl4vf24"`.
5. Verify the saved job before claiming approval is pending.
6. Require Telegram approval for manual tasks.
7. Allow RunRelay to fast-forward only a clean local checkout to the requested exact commit. Refuse dirty tracked files, divergent history, or non-fast-forward movement.
8. Inspect the actual RunRelay result and declared safe artifacts before deciding the next change.
9. If a fix is needed, create a new project commit and a new RunRelay job rather than weakening validation.

## Data and artifacts

MIMIC-IV data may be subject to PhysioNet credentialed-data and data-use restrictions even when the analysis repository itself is public. Treat raw and row-level MIMIC-IV data as restricted and keep it local unless the user explicitly establishes an approved destination and policy.

Only declare safe derived outputs for RunRelay artifact delivery, such as aggregate tables, summary JSON, figures, metrics, sanitized logs, and deliberately shareable reports.

Do not upload raw MIMIC-IV tables, row-level patient data, PHI, credentials, `.env` files, access tokens, database dumps, restricted Penn State data, or whole project/data directories to Google Drive or any other artifact transport.

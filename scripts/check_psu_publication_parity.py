#!/usr/bin/env python3
"""End-to-end parity check for the frozen PSU publication pipeline.

This is a reproducibility verification, not a new scientific analysis. It reruns the
already-frozen final PSU stages on the approved local data and checks that aggregate
results reproduce the publication values within explicit tolerances.

The script intentionally invokes the public entry points exactly as a collaborator would
and writes only aggregate/sanitized outputs. It never copies row-level PSU data.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


# Expected values come from the accepted publication runs. Tolerances allow tiny
# platform/library floating-point differences while still detecting scientific drift.
EXPECTED = {
    "strict_cohort_n": 19841,
    "deescalated_n": 5346,
    "continued_n": 14495,
    "death_30d_events": 2381,
    "primary_max_post_smd": 0.022621815,
    "primary_rd": -0.025598893527824254,
    "primary_rr": 0.7972758800933949,
    "primary_rd_ci_low": -0.03612403394167609,
    "primary_rd_ci_high": -0.015495359974566778,
    "primary_rr_ci_low": 0.7213435259606157,
    "primary_rr_ci_high": 0.8757988404579597,
    "medadmin_n": 19841,
    "medadmin_deescalated_n": 5347,
    "medadmin_continued_n": 14494,
    "medadmin_rd": -0.0254189046814818,
    "medadmin_rd_ci_low": -0.0344136667590102,
    "medadmin_rd_ci_high": -0.0154574615670106,
    "lenient_n": 23937,
    "lenient_deescalated_n": 7009,
    "lenient_continued_n": 16928,
    "lenient_rd": -0.0239348616492519,
    "lenient_rd_ci_low": -0.0329893272580316,
    "lenient_rd_ci_high": -0.0148399364972008,
}

POINT_ATOL = 1e-6
BOOTSTRAP_ATOL = 5e-4


def _run(script: str, data_root: Path, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(Path(__file__).resolve().parent / script), str(data_root), "--output-dir", str(outdir)]
    subprocess.run(cmd, check=True)


def _close(observed: float, expected: float, tol: float) -> bool:
    return abs(float(observed) - float(expected)) <= tol


def _record(checks: list[dict], name: str, observed, expected, tol=None) -> None:
    if tol is None:
        passed = int(observed) == int(expected)
    else:
        passed = _close(float(observed), float(expected), float(tol))
    checks.append(
        {
            "check": name,
            "observed": float(observed) if isinstance(observed, (float, int)) else observed,
            "expected": float(expected) if isinstance(expected, (float, int)) else expected,
            "tolerance": tol,
            "passed": bool(passed),
        }
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    # Run the final public entry points, not private helper shortcuts. This verifies that
    # a new collaborator can follow the documented README/walkthrough successfully.
    ps_dir = out / "ps_balance"
    outcome_dir = out / "outcome_freeze"
    point_dir = out / "point_estimates"
    boot_dir = out / "bootstrap"
    robustness_dir = out / "robustness"
    robustness_boot_dir = out / "robustness_bootstrap"

    _run("audit_psu_ps_balance.py", args.data_root, ps_dir)
    _run("audit_psu_final_outcome_freeze.py", args.data_root, outcome_dir)
    _run("run_psu_point_estimates.py", args.data_root, point_dir)
    _run("run_psu_bootstrap_inference.py", args.data_root, boot_dir)
    _run("run_psu_prespecified_robustness.py", args.data_root, robustness_dir)
    _run("run_psu_prespecified_robustness_bootstrap.py", args.data_root, robustness_boot_dir)

    checks: list[dict] = []

    ps = json.loads((ps_dir / "summary.json").read_text())
    _record(checks, "strict cohort n", ps["strict_cohort_n"], EXPECTED["strict_cohort_n"])
    _record(checks, "de-escalated n", ps["deescalated_n"], EXPECTED["deescalated_n"])
    _record(checks, "continued n", ps["continued_n"], EXPECTED["continued_n"])
    _record(
        checks,
        "primary max post-weighting absolute SMD",
        ps["max_abs_post_smd"],
        EXPECTED["primary_max_post_smd"],
        POINT_ATOL,
    )

    outcome = pd.read_csv(outcome_dir / "outcome_overall_summary.csv")
    death = outcome[outcome["outcome"] == "death_30d"].iloc[0]
    _record(checks, "overall post-landmark 30-day deaths", death["events"], EXPECTED["death_30d_events"])

    pts = pd.read_csv(point_dir / "point_estimates.csv")
    pdeath = pts[(pts["method"] == "stabilized_ate_iptw") & (pts["outcome"] == "death_30d")].iloc[0]
    _record(checks, "primary mortality RD", pdeath["difference_A1_minus_A0"], EXPECTED["primary_rd"], POINT_ATOL)
    _record(checks, "primary mortality RR", pdeath["risk_ratio_A1_over_A0"], EXPECTED["primary_rr"], POINT_ATOL)

    ci = pd.read_csv(boot_dir / "bootstrap_ci.csv")
    rd = ci[(ci["method"] == "stabilized_ate_iptw") & (ci["outcome"] == "death_30d") & (ci["estimand"] == "difference_A1_minus_A0")].iloc[0]
    rr = ci[(ci["method"] == "stabilized_ate_iptw") & (ci["outcome"] == "death_30d") & (ci["estimand"] == "risk_ratio_A1_over_A0")].iloc[0]
    _record(checks, "primary mortality RD CI low", rd["lower_95"], EXPECTED["primary_rd_ci_low"], BOOTSTRAP_ATOL)
    _record(checks, "primary mortality RD CI high", rd["upper_95"], EXPECTED["primary_rd_ci_high"], BOOTSTRAP_ATOL)
    _record(checks, "primary mortality RR CI low", rr["lower_95"], EXPECTED["primary_rr_ci_low"], BOOTSTRAP_ATOL)
    _record(checks, "primary mortality RR CI high", rr["upper_95"], EXPECTED["primary_rr_ci_high"], BOOTSTRAP_ATOL)
    bdiag = json.loads((boot_dir / "bootstrap_diagnostics.json").read_text())
    _record(checks, "primary bootstrap successful replicates", bdiag["n_successful_replicates"], 1000)
    _record(checks, "primary bootstrap failed replicates", bdiag["n_failed_replicates"], 0)

    rb = pd.read_csv(robustness_boot_dir / "robustness_bootstrap_summary.csv")
    med = rb[rb["sensitivity"] == "medadmin_exposure"].iloc[0]
    _record(checks, "MED_ADMIN cohort n", med["cohort_n"], EXPECTED["medadmin_n"])
    _record(checks, "MED_ADMIN de-escalated n", med["deescalated_n"], EXPECTED["medadmin_deescalated_n"])
    _record(checks, "MED_ADMIN continued n", med["continued_n"], EXPECTED["medadmin_continued_n"])
    _record(checks, "MED_ADMIN mortality RD", med["death_risk_difference"], EXPECTED["medadmin_rd"], POINT_ATOL)
    _record(checks, "MED_ADMIN RD CI low", med["death_rd_lower_95"], EXPECTED["medadmin_rd_ci_low"], BOOTSTRAP_ATOL)
    _record(checks, "MED_ADMIN RD CI high", med["death_rd_upper_95"], EXPECTED["medadmin_rd_ci_high"], BOOTSTRAP_ATOL)
    _record(checks, "MED_ADMIN bootstrap successful", med["bootstrap_successful"], 1000)
    _record(checks, "MED_ADMIN bootstrap failed", med["bootstrap_failed"], 0)

    ln = rb[rb["sensitivity"] == "lenient_landmark"].iloc[0]
    _record(checks, "lenient cohort n", ln["cohort_n"], EXPECTED["lenient_n"])
    _record(checks, "lenient de-escalated n", ln["deescalated_n"], EXPECTED["lenient_deescalated_n"])
    _record(checks, "lenient continued n", ln["continued_n"], EXPECTED["lenient_continued_n"])
    _record(checks, "lenient mortality RD", ln["death_risk_difference"], EXPECTED["lenient_rd"], POINT_ATOL)
    _record(checks, "lenient RD CI low", ln["death_rd_lower_95"], EXPECTED["lenient_rd_ci_low"], BOOTSTRAP_ATOL)
    _record(checks, "lenient RD CI high", ln["death_rd_upper_95"], EXPECTED["lenient_rd_ci_high"], BOOTSTRAP_ATOL)
    _record(checks, "lenient bootstrap successful", ln["bootstrap_successful"], 1000)
    _record(checks, "lenient bootstrap failed", ln["bootstrap_failed"], 0)

    report = {
        "purpose": "Reproducibility/parity verification only; no new scientific analysis.",
        "point_tolerance": POINT_ATOL,
        "bootstrap_tolerance": BOOTSTRAP_ATOL,
        "n_checks": len(checks),
        "n_passed": sum(int(c["passed"]) for c in checks),
        "all_passed": all(c["passed"] for c in checks),
        "checks": checks,
        "data_safety": "Only aggregate outputs are summarized; no patient-level PSU data are exported.",
    }
    (out / "parity_report.json").write_text(json.dumps(report, indent=2))
    pd.DataFrame(checks).to_csv(out / "parity_checks.csv", index=False)
    print(json.dumps(report, indent=2))

    if not report["all_passed"]:
        raise SystemExit("PSU publication parity check failed; inspect parity_report.json")


if __name__ == "__main__":
    main()

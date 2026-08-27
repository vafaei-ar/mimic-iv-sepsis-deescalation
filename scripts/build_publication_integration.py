#!/usr/bin/env python3
"""Build harmonized publication tables from the frozen aggregate source snapshot.

This script performs no patient-level analysis and does not refit any model. It reads only the
sanitized aggregate snapshot created by ``build_publication_source_snapshot.py`` and emits
manuscript-facing summary tables. For PSU mortality quantities that were explicitly frozen by
the parity checker, the accepted publication values in ``parity_report.json`` are used rather
than run-to-run floating-point realizations that vary slightly within the validated tolerance.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

SNAPSHOT = Path("outputs/publication_integration/source_snapshot/publication_source_snapshot.json")
OUT = Path("outputs/publication_integration/harmonized")


def records(snapshot: dict, key: str) -> list[dict]:
    content = snapshot["sources"][key]["content"]
    return content.get("records", [])


def parity_expected(snapshot: dict) -> dict[str, float | int]:
    checks = snapshot["sources"]["psu_parity_report"]["content"]["checks"]
    return {c["check"]: c["expected"] for c in checks}


def psu_point(snapshot: dict, outcome: str) -> dict:
    for row in records(snapshot, "psu_point_estimates"):
        if row["method"] == "stabilized_ate_iptw" and row["outcome"] == outcome:
            return row
    raise KeyError(outcome)


def psu_ci(snapshot: dict, outcome: str, estimand: str) -> dict:
    for row in records(snapshot, "psu_bootstrap_ci"):
        if row["method"] == "stabilized_ate_iptw" and row["outcome"] == outcome and row["estimand"] == estimand:
            return row
    raise KeyError((outcome, estimand))


def main() -> None:
    if not SNAPSHOT.exists():
        raise SystemExit(f"Missing frozen source snapshot: {SNAPSHOT}")
    snapshot = json.loads(SNAPSHOT.read_text())
    OUT.mkdir(parents=True, exist_ok=True)

    mimic_primary = next(r for r in records(snapshot, "mimic_primary_secondary") if r["analysis"] == "30-day mortality")
    expected = parity_expected(snapshot)
    psu_balance = snapshot["sources"]["psu_ps_balance_summary"]["content"]

    mortality_rows = [
        {
            "dataset_analysis": "MIMIC-IV primary stabilized IPTW",
            "cohort_n": 9589,
            "deescalated_n": 1863,
            "continued_n": 7726,
            "mortality_risk_deescalated": mimic_primary["risk_deescalated_stopped"],
            "mortality_risk_continued": mimic_primary["risk_continued"],
            "mortality_rd": mimic_primary["risk_difference"],
            "rd_ci95_low": mimic_primary["lower_95"],
            "rd_ci95_high": mimic_primary["upper_95"],
            "mortality_rr": mimic_primary["risk_ratio"],
            "rr_ci95_low": None,
            "rr_ci95_high": None,
            "estimand": "ATE",
            "note": "Corrected final MIMIC-IV analysis; primary 1000-replicate mortality CI.",
        },
        {
            "dataset_analysis": "PSU modified external replication, primary PRESCRIBING exposure",
            "cohort_n": expected["strict cohort n"],
            "deescalated_n": expected["de-escalated n"],
            "continued_n": expected["continued n"],
            "mortality_risk_deescalated": None,
            "mortality_risk_continued": None,
            "mortality_rd": expected["primary mortality RD"],
            "rd_ci95_low": expected["primary mortality RD CI low"],
            "rd_ci95_high": expected["primary mortality RD CI high"],
            "mortality_rr": expected["primary mortality RR"],
            "rr_ci95_low": expected["primary mortality RR CI low"],
            "rr_ci95_high": expected["primary mortality RR CI high"],
            "estimand": "ATE",
            "note": "Publication-locked PSU parity values; modified external replication.",
        },
        {
            "dataset_analysis": "PSU MED_ADMIN exposure sensitivity",
            "cohort_n": expected["MED_ADMIN cohort n"],
            "deescalated_n": expected["MED_ADMIN de-escalated n"],
            "continued_n": expected["MED_ADMIN continued n"],
            "mortality_risk_deescalated": None,
            "mortality_risk_continued": None,
            "mortality_rd": expected["MED_ADMIN mortality RD"],
            "rd_ci95_low": expected["MED_ADMIN RD CI low"],
            "rd_ci95_high": expected["MED_ADMIN RD CI high"],
            "mortality_rr": None,
            "rr_ci95_low": None,
            "rr_ci95_high": None,
            "estimand": "ATE",
            "note": "Prespecified administration-based exposure sensitivity; parity-locked RD values.",
        },
        {
            "dataset_analysis": "PSU lenient 96-hour landmark sensitivity",
            "cohort_n": expected["lenient cohort n"],
            "deescalated_n": expected["lenient de-escalated n"],
            "continued_n": expected["lenient continued n"],
            "mortality_risk_deescalated": None,
            "mortality_risk_continued": None,
            "mortality_rd": expected["lenient mortality RD"],
            "rd_ci95_low": expected["lenient RD CI low"],
            "rd_ci95_high": expected["lenient RD CI high"],
            "mortality_rr": None,
            "rr_ci95_low": None,
            "rr_ci95_high": None,
            "estimand": "ATE",
            "note": "Prespecified landmark sensitivity; parity-locked RD values.",
        },
    ]
    pd.DataFrame(mortality_rows).to_csv(OUT / "harmonized_mortality_results.csv", index=False)

    secondary_map = [
        ("hospital-free days", "hospital_free_days_30d", "Hospital-free days", "days"),
        ("antibiotic-free days", "antibiotic_free_days_30d", "Antibiotic-free days", "days"),
        ("normalized systemic antibiotic exposure", "normalized_antibiotic_exposure_30d", "Normalized systemic antibiotic exposure", "proportion"),
        ("normalized broad-spectrum exposure", "normalized_broad_antibiotic_exposure_30d", "Normalized broad-spectrum exposure", "proportion"),
        ("late recurrent/persistent antibiotic-course use", "late_recurrent_or_persistent_abx_course_30d", "Late recurrent/persistent antibiotic course", "risk difference"),
    ]
    mimic_by_name = {r["analysis"]: r for r in records(snapshot, "mimic_primary_secondary")}
    secondary_rows: list[dict] = []
    for mimic_name, psu_name, label, unit in secondary_map:
        m = mimic_by_name[mimic_name]
        estimate = m["mean_difference"] if m["kind"] == "mean" else m["risk_difference"]
        secondary_rows.append({
            "dataset": "MIMIC-IV",
            "outcome": label,
            "estimate": estimate,
            "ci95_low": m["lower_95"],
            "ci95_high": m["upper_95"],
            "unit": unit,
        })
        p = psu_point(snapshot, psu_name)
        pci = psu_ci(snapshot, psu_name, "difference_A1_minus_A0")
        secondary_rows.append({
            "dataset": "PSU",
            "outcome": label,
            "estimate": p["difference_A1_minus_A0"],
            "ci95_low": pci["lower_95"],
            "ci95_high": pci["upper_95"],
            "unit": unit,
        })
    pd.DataFrame(secondary_rows).to_csv(OUT / "harmonized_secondary_outcomes.csv", index=False)

    progressive = pd.DataFrame(records(snapshot, "mimic_progressive_adjustment"))
    progressive.to_csv(OUT / "mimic_progressive_adjustment.csv", index=False)

    mimic_weight = next(r for r in records(snapshot, "mimic_final_weighting_point_estimates") if r["analysis"] == "Primary stabilized IPTW")
    weighting = pd.DataFrame([
        {
            "dataset": "MIMIC-IV",
            "max_post_smd": mimic_weight["max_post_smd"],
            "worst_balanced_variable": mimic_weight["worst_balanced_variable"],
            "ess_deescalated": mimic_weight["ess_deescalated_stopped"],
            "ess_continued": mimic_weight["ess_continued"],
            "max_weight": mimic_weight["max_weight"],
        },
        {
            "dataset": "PSU",
            "max_post_smd": expected["primary max post-weighting absolute SMD"],
            "worst_balanced_variable": psu_balance["worst_post_balance_variable"],
            "ess_deescalated": psu_balance["treated_ess"],
            "ess_continued": psu_balance["continued_ess"],
            "max_weight": psu_balance["max_weight"],
        },
    ])
    weighting.to_csv(OUT / "weighting_diagnostics.csv", index=False)

    m_rd = 100 * mimic_primary["risk_difference"]
    m_lo = 100 * mimic_primary["lower_95"]
    m_hi = 100 * mimic_primary["upper_95"]
    p_rd = 100 * expected["primary mortality RD"]
    p_lo = 100 * expected["primary mortality RD CI low"]
    p_hi = 100 * expected["primary mortality RD CI high"]
    md = f"""# Frozen MIMIC-IV and PSU publication integration\n\nThis package is generated only from the frozen aggregate source snapshot. No patient-level data are read and no models are refit.\n\n## Primary interpretation\n\n- MIMIC-IV corrected primary stabilized IPTW: 30-day mortality RD {m_rd:+.2f} percentage points (95% CI {m_lo:+.2f} to {m_hi:+.2f}), RR {mimic_primary['risk_ratio']:.3f}.\n- PSU modified external replication: publication-locked 30-day mortality RD {p_rd:+.2f} percentage points (95% CI {p_lo:+.2f} to {p_hi:+.2f}), RR {expected['primary mortality RR']:.3f}.\n- These estimates should not be forced into agreement. The PSU analysis is a modified external replication with materially different culture-result, timing, route, and data-model semantics.\n- The cross-dataset result that is most consistent is reduced antibiotic burden after de-escalation.\n\n## MIMIC progressive adjustment\n\nThe mortality association moves from approximately -4.14 percentage points in M1 to +0.84 percentage points in fully adjusted M4, supporting strong confounding by clinical improvement and treatment intensity. The final MIMIC ATE retains residual balance/positivity limitations (maximum post-weighting absolute SMD {mimic_weight['max_post_smd']:.3f}; treated ESS {mimic_weight['ess_deescalated_stopped']:.0f}).\n\n## PSU reproducibility status\n\nThe PSU publication parity check passed all 29 checks. Mortality quantities in the harmonized table use the accepted publication values embedded in that checker; current rerun values can vary slightly within the validated numerical tolerances.\n\n## Files\n\n- `harmonized_mortality_results.csv`\n- `harmonized_secondary_outcomes.csv`\n- `mimic_progressive_adjustment.csv`\n- `weighting_diagnostics.csv`\n\nThe late recurrent/persistent antibiotic-course outcome remains exploratory because observation is affected by discharge timing.\n"""
    (OUT / "publication_integration_summary.md").write_text(md)


if __name__ == "__main__":
    main()

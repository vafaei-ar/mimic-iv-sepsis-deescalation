#!/usr/bin/env python3
from __future__ import annotations

import json
import pandas as pd

import build_reviewer_support_package as base
from sepsis_deescalation.fast_bootstrap import (
    OutcomeSpec,
    bootstrap_ci_long,
    bootstrap_multi_outcome_iptw,
)
from sepsis_deescalation.specification import CANDIDATE_PS_VARS
from sepsis_deescalation.stats import (
    balance_table,
    effective_sample_size,
    fit_stabilized_iptw,
    risks,
)


def subtype_sensitivity_fast(d: pd.DataFrame, reps: int = 1000) -> pd.DataFrame:
    contrasts = [
        ("narrowed_or_non_broad_only", "Narrowed/non-broad only vs continued broad"),
        ("stopped_all_observed_systemic_antibiotics", "Stopped all observed systemic antibiotics vs continued broad"),
    ]
    rows = []
    diag_rows = []
    for i, (typ, label) in enumerate(contrasts):
        s = d.loc[d["deescalation_type"].isin([typ, "continued_broad"])].copy()
        s["A"] = (s["deescalation_type"] == typ).astype(int)
        w, _, _ = fit_stabilized_iptw(s, CANDIDATE_PS_VARS)
        rt, rc, rd, rr = risks(w, "death_by_horizon", "SW_A")
        bal = balance_table(w, CANDIDATE_PS_VARS)

        boot, diag = bootstrap_multi_outcome_iptw(
            s,
            CANDIDATE_PS_VARS,
            [OutcomeSpec(label=label, column="death_by_horizon", kind="risk")],
            reps=reps,
            seed=20260827 + i,
            jobs="auto",
        )
        ci = bootstrap_ci_long(boot, requested=reps)
        rdrow = ci.loc[(ci["analysis"] == label) & (ci["estimand"] == "risk_difference")].iloc[0]
        rrrow = ci.loc[(ci["analysis"] == label) & (ci["estimand"] == "risk_ratio")].iloc[0]
        t = w["A"] == 1
        c = ~t
        rows.append({
            "analysis": label,
            "n_total": len(w),
            "n_subtype": int(t.sum()),
            "n_continued": int(c.sum()),
            "weighted_risk_subtype": rt,
            "weighted_risk_continued": rc,
            "risk_difference": rd,
            "rd_ci95_low": float(rdrow["lower_95"]),
            "rd_ci95_high": float(rdrow["upper_95"]),
            "risk_ratio": rr,
            "rr_ci95_low": float(rrrow["lower_95"]),
            "rr_ci95_high": float(rrrow["upper_95"]),
            "bootstrap_successes": int(rdrow["n_success"]),
            "max_post_smd": float(bal["after"].max()),
            "ess_subtype": effective_sample_size(w.loc[t, "SW_A"]),
            "ess_continued": effective_sample_size(w.loc[c, "SW_A"]),
            "max_weight": float(w["SW_A"].max()),
            "interpretation": "Exploratory subtype sensitivity; not a replacement for the frozen binary primary estimand.",
        })
        dd = diag.iloc[0].to_dict()
        dd["analysis"] = label
        diag_rows.append(dd)

    out = pd.DataFrame(rows)
    out.to_csv(base.OUT / "mimic_deescalation_subtype_adjusted_mortality.csv", index=False)
    pd.DataFrame(diag_rows).to_csv(base.OUT / "mimic_deescalation_subtype_bootstrap_diagnostics.csv", index=False)
    return out


def main() -> None:
    base.OUT.mkdir(parents=True, exist_ok=True)
    for p in [
        base.COHORT,
        base.FLOW,
        base.HARM / "harmonized_mortality_results.csv",
        base.HARM / "harmonized_secondary_outcomes.csv",
        base.HARM / "mimic_progressive_adjustment.csv",
    ]:
        if not p.exists():
            raise FileNotFoundError(p)

    d = pd.read_csv(base.COHORT, low_memory=False)
    w, bal = base.validate_primary(d)
    flow = base.build_flow()
    desc = base.subtype_descriptives(d)

    # Generate the inexpensive diagnostics and publication figures before the
    # exploratory subtype bootstrap so useful artifacts survive any later issue.
    base.diagnostics(w, bal)
    base.revised_main_figures()

    sub = subtype_sensitivity_fast(d, reps=1000)
    base.summary(flow, desc, sub, bal)
    metadata = {
        "corrected_cohort": str(base.COHORT),
        "n": len(d),
        "deescalated": int((d["A"] == 1).sum()),
        "continued": int((d["A"] == 0).sum()),
        "subtype_bootstrap_reps": 1000,
        "subtype_bootstrap_engine": "parallel matrix-form fast bootstrap validated against historical implementation",
        "primary_science_changed": False,
    }
    (base.OUT / "reviewer_support_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

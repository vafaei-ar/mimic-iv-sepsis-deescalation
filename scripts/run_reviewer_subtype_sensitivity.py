#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

import pandas as pd

import build_reviewer_support_package as base
from sepsis_deescalation.fast_bootstrap import OutcomeSpec, bootstrap_ci_long, bootstrap_multi_outcome_iptw
from sepsis_deescalation.specification import CANDIDATE_PS_VARS
from sepsis_deescalation.stats import balance_table, effective_sample_size, fit_stabilized_iptw, risks

CONTRASTS = {
    "narrowed": (
        "narrowed_or_non_broad_only",
        "Narrowed/non-broad only vs continued broad",
        20260827,
    ),
    "stopped": (
        "stopped_all_observed_systemic_antibiotics",
        "Stopped all observed systemic antibiotics vs continued broad",
        20260828,
    ),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contrast", choices=sorted(CONTRASTS), required=True)
    parser.add_argument("--reps", type=int, default=1000)
    args = parser.parse_args()

    base.OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(base.COHORT, low_memory=False)
    base.validate_primary(d)

    typ, label, seed = CONTRASTS[args.contrast]
    s = d.loc[d["deescalation_type"].isin([typ, "continued_broad"])].copy()
    s["A"] = (s["deescalation_type"] == typ).astype(int)

    w, _, _ = fit_stabilized_iptw(s, CANDIDATE_PS_VARS)
    rt, rc, rd, rr = risks(w, "death_by_horizon", "SW_A")
    bal = balance_table(w, CANDIDATE_PS_VARS)

    boot, diag = bootstrap_multi_outcome_iptw(
        s,
        CANDIDATE_PS_VARS,
        [OutcomeSpec(label=label, column="death_by_horizon", kind="risk")],
        reps=args.reps,
        seed=seed,
        jobs="auto",
    )
    ci = bootstrap_ci_long(boot, requested=args.reps)
    rdrow = ci.loc[(ci["analysis"] == label) & (ci["estimand"] == "risk_difference")].iloc[0]
    rrrow = ci.loc[(ci["analysis"] == label) & (ci["estimand"] == "risk_ratio")].iloc[0]

    t = w["A"] == 1
    c = ~t
    row = {
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
    }

    stem = f"mimic_deescalation_subtype_{args.contrast}"
    pd.DataFrame([row]).to_csv(base.OUT / f"{stem}_adjusted_mortality.csv", index=False)
    dd = diag.iloc[0].to_dict()
    dd["analysis"] = label
    pd.DataFrame([dd]).to_csv(base.OUT / f"{stem}_bootstrap_diagnostics.csv", index=False)
    metadata = {
        "contrast": args.contrast,
        "label": label,
        "bootstrap_reps": args.reps,
        "seed": seed,
        "primary_science_changed": False,
    }
    (base.OUT / f"{stem}_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

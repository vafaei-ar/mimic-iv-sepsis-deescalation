#!/usr/bin/env python3
"""Prespecified PSU robustness point-estimate bundle.

This script runs exactly two robustness analyses that were identified before PSU
outcome-effect inspection:

1. MED_ADMIN exposure sensitivity on the frozen strict cohort.
2. Lenient 96-hour landmark sensitivity using >= landmark calendar date rather than >.

Why these two sensitivities
---------------------------
PRESCRIBING is the primary PSU source because it is the closest available analogue to
the frozen MIMIC prescription/order phenotype. MED_ADMIN is valuable as a measurement
sensitivity because it asks whether the result depends on ordered versus recorded
administration data. The landmark sensitivity addresses the known limitation that PSU
only provides date-level discharge/death timing, making patients discharged or dying on
the landmark calendar date ambiguous relative to an exact 96-hour timestamp.

Everything else is intentionally frozen: covariates, propensity-score preprocessing and
model, outcome definitions, and weighting estimand. Only aggregate outputs are written.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

import run_psu_point_estimates as pe


# These literal markers are deliberate fail-closed guardrails. The robustness scripts
# modify only the exact frozen clauses shown below. If the parent PS script changes and
# a marker no longer appears exactly once, execution stops rather than silently changing
# more of the analysis than intended.
DF_MARKER = '    df = con.execute("select c.patid,c.encounterid,c.anchor_ts,c.admit_date,case when e.broad_72_96=1 then 0 else 1 end A from cohort c join exposure e using(patid,encounterid)").fetchdf()\n'
BAL_MARKER = '    max_pre=float(bdf.abs_pre_smd.max()); max_post=float(bdf.abs_post_smd.max()); worst=str(bdf.iloc[0].variable)\n'
STRICT_COHORT = "(b.discharge_date is null or b.discharge_date>cast(a.anchor_ts+interval 96 hour as date)) and (d.death_date is null or d.death_date>cast(a.anchor_ts+interval 96 hour as date))"
LENIENT_COHORT = "(b.discharge_date is null or b.discharge_date>=cast(a.anchor_ts+interval 96 hour as date)) and (d.death_date is null or d.death_date>=cast(a.anchor_ts+interval 96 hour as date))"


MEDADMIN_OVERRIDE = r'''
    # Prespecified MED_ADMIN exposure sensitivity.
    #
    # The PSU MED_ADMIN audit showed that the raw start/stop TIME fields are not usable
    # as independent clock times. Therefore, and *before looking at treatment effects*,
    # we froze a conservative date-span interpretation: start date begins at 00:00:00;
    # stop date ends at 23:59:59; a missing stop date is treated as the end of the start
    # calendar day. This sensitivity changes the exposure measurement source only.
    mcols=cols["med"]
    mp=first(mcols,["PATID"]); me=first(mcols,["ENCOUNTERID","ENCOUNTER_ID"])
    msd=first(mcols,["MEDADMIN_START_DATE","START_DATE"]); med=first(mcols,["MEDADMIN_STOP_DATE","STOP_DATE"])
    mn=first(mcols,["RAW_MEDADMIN_MED_NAME","MEDADMIN_MED_NAME","MEDICATION_NAME"])
    mc=first(mcols,["MEDADMIN_CODE","RXNORM_CUI"]); mr=first(mcols,["MEDADMIN_ROUTE","ROUTE"])
    if any(v is None for v in [mp,me,msd,mn]):
        raise RuntimeError("Required MED_ADMIN exposure fields missing")
    mname=f"lower(coalesce(cast(m.{qi(mn)} as varchar),''))"
    mcode=f"trim(coalesce(cast(m.{qi(mc)} as varchar),''))" if mc else "''"
    mroute=f"upper(trim(coalesce(cast(m.{qi(mr)} as varchar),'')))" if mr else "''"
    mbroad=f"(({mcode} in ({incsql}) or regexp_matches({mname},{q(BROAD_PATTERN)})) and {mcode} not in ({excsql}) and not regexp_matches({mname},{q(NON_SYSTEMIC_PATTERN)}) and {mroute} not in ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    mstart=f"cast(try_cast(m.{qi(msd)} as date) as timestamp)"
    if med:
        mstop=f"coalesce(cast(try_cast(m.{qi(med)} as date) as timestamp)+interval 23 hour+interval 59 minute+interval 59 second,{mstart}+interval 23 hour+interval 59 minute+interval 59 second)"
    else:
        mstop=f"{mstart}+interval 23 hour+interval 59 minute+interval 59 second"
    con.execute(f"""create temp table medadmin_exposure_sensitivity as
        select c.patid,c.encounterid,
        max(case when {mbroad} and {mstart}<c.anchor_ts+interval 96 hour and {mstop}>=c.anchor_ts+interval 72 hour then 1 else 0 end) broad_72_96
        from cohort c left join med m
          on cast(m.{qi(mp)} as varchar)=c.patid and cast(m.{qi(me)} as varchar)=c.encounterid
        group by 1,2""")
    ma=con.execute("select patid,encounterid,case when broad_72_96=1 then 0 else 1 end A_medadmin from medadmin_exposure_sensitivity").fetchdf()
    df=df.drop(columns=["A"]).merge(ma,on=["patid","encounterid"],how="left",validate="one_to_one")
    df["A"]=pd.to_numeric(df.pop("A_medadmin"),errors="raise").astype(int)
'''


def run_variant(source: str, variant: str, data_root: Path, outdir: Path) -> dict:
    s = source
    if variant == "medadmin_exposure":
        if s.count(DF_MARKER) != 1:
            raise RuntimeError("Frozen PS dataframe marker missing or ambiguous")
        s = s.replace(DF_MARKER, DF_MARKER + MEDADMIN_OVERRIDE + "\n", 1)
    elif variant == "lenient_landmark":
        if s.count(STRICT_COHORT) != 1:
            raise RuntimeError("Frozen strict landmark clause missing or ambiguous")
        s = s.replace(STRICT_COHORT, LENIENT_COHORT, 1)
    else:
        raise ValueError(variant)

    if s.count(BAL_MARKER) != 1:
        raise RuntimeError("Frozen PS balance marker missing or ambiguous")
    s = s.replace(BAL_MARKER, pe.INJECT + "\n" + BAL_MARKER, 1)

    outdir.mkdir(parents=True, exist_ok=True)
    old_argv = sys.argv[:]
    try:
        sys.argv = ["audit_psu_ps_balance.py", str(data_root), "--output-dir", str(outdir)]
        ns = {"__name__": "__main__", "__file__": str(Path(__file__).resolve().parent / "audit_psu_ps_balance.py")}
        exec(compile(s, ns["__file__"], "exec"), ns, ns)
    finally:
        sys.argv = old_argv

    summary = json.loads((outdir / "effect_summary.json").read_text())
    pts = pd.read_csv(outdir / "point_estimates.csv")
    death = pts[(pts["method"] == "stabilized_ate_iptw") & (pts["outcome"] == "death_30d")].iloc[0]
    w = pd.read_csv(outdir / "effect_weighting_diagnostics.csv")
    primary_w = w[w["method"] == "stabilized_ate_iptw"].iloc[0]
    return {
        "sensitivity": variant,
        "cohort_n": int(summary["strict_cohort_n"]),
        "deescalated_n": int(summary["deescalated_n"]),
        "continued_n": int(summary["continued_n"]),
        "death_risk_deescalated": float(death["deescalated_mean_or_risk"]),
        "death_risk_continued": float(death["continued_mean_or_risk"]),
        "death_risk_difference": float(death["difference_A1_minus_A0"]),
        "death_risk_ratio": float(death["risk_ratio_A1_over_A0"]),
        "treated_ess": float(primary_w["treated_ess"]),
        "continued_ess": float(primary_w["continued_ess"]),
        "max_weight": float(primary_w["max_weight"]),
        "max_abs_post_smd": float(primary_w["max_abs_post_smd"]),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    source_path = Path(__file__).resolve().parent / "audit_psu_ps_balance.py"
    source = source_path.read_text(encoding="utf-8")

    rows = [
        run_variant(source, "medadmin_exposure", args.data_root, args.output_dir / "medadmin_exposure"),
        run_variant(source, "lenient_landmark", args.data_root, args.output_dir / "lenient_landmark"),
    ]
    pd.DataFrame(rows).to_csv(args.output_dir / "robustness_summary.csv", index=False)

    meta = {
        "privacy_mode": "aggregate_only",
        "sensitivities": ["medadmin_exposure", "lenient_landmark"],
        "medadmin_rule": "strict primary PRESCRIBING-anchored cohort; binary day-3 exposure reclassified from MED_ADMIN using date-span intervals",
        "lenient_landmark_rule": "PRESCRIBING exposure; discharge/death on landmark calendar date allowed with >= date rule",
        "guardrail": "These sensitivities were identified before PSU outcome-effect inspection. Frozen covariate, PS, and outcome definitions are otherwise unchanged. Point estimates only; no bootstrap in this task.",
    }
    (args.output_dir / "robustness_metadata.json").write_text(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()

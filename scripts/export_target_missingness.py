#!/usr/bin/env python3
"""Export TARGET item 11 group-wise missingness from frozen, aggregate baseline files.

Does not access patient-level data, refit models, or modify frozen sources.
Only reports original missingness when preserved indicators or un-imputed
pre-propensity values establish it; otherwise reports 'NR' explicitly.
Small nonzero cells (1-10) are suppressed for safe sharing.
"""
from __future__ import annotations

import json
from pathlib import Path
import pandas as pd

IN = Path("outputs/publication_integration/baseline_characteristics")
OUT = Path("outputs/publication_integration/target_missingness")
EXPECTED = {"MIMIC-IV": (9589,1863,7726), "Penn State": (19841,5346,14495)}
MIN_CELL = 11

def num(x):
    try:
        v = int(float(x))
        return v
    except (ValueError,TypeError,OverflowError):
        return None

def safe_count(n, denominator):
    if n is None:
        return "NR", "NR"
    if not 0 <= n <= denominator:
        raise ValueError("Count violates frozen treatment-group denominator")
    if 0 < n < MIN_CELL:
        return "<11", "Suppressed"
    return str(n), f"{n/denominator*100:.1f}%"

def main():
    outrows = []
    tally = {}
    for source, basename in [
       ("MIMIC-IV", "mimic_baseline_characteristics"),
       ("Penn State","psu_baseline_characteristics")
    ]:
        dpath=IN / (basename+"_detailed.csv")
        mpath=IN / (basename+"_metadata.json")
        if not dpath.is_file() or not mpath.is_file():
            raise FileNotFoundError("Frozen aggregate baseline output missing: "+basename)
        df=pd.read_csv(dpath)
        meta=json.loads(mpath.read_text(encoding="utf-8"))
        n,n1,n0=EXPECTED[source]
        if int(meta["cohort_n"])!=n or int(meta["deescalation_or_stopping_n"])!=n1 or int(meta["continued_broad_n"])!=n0:
            raise RuntimeError("Frozen cohort parity failure for "+source)
        recovered=set(meta.get("descriptive_imputation_correction",{}).get("variables",[]))
        if source=="MIMIC-IV" and not meta.get("descriptive_imputation_correction",{}).get("applied",False):
            raise RuntimeError("Source not corrected for source-visible missingness")
        for row in df.to_dict("records"):
            var=str(row["variable"])
            typ=str(row["type"])
            if source=="MIMIC-IV":
                a=num(row.get("deescalated_or_stopped_missing_n"))
                b=num(row.get("continued_broad_missing_n"))
                if typ=="continuous" and var in recovered:
                    status="Original missingness recovered from stored flag"
                    complete=True
                elif typ=="continuous":
                    status="Original missingness not independently recoverable after source preprocessing"
                    complete=False
                else:
                    status="Original binary source missingness not separable from zero-coded category"
                    complete=False
            else:
                a=num(row.get("deescalated_missing_n"))
                b=num(row.get("continued_missing_n"))
                if typ=="continuous":
                    status="Missing in un-imputed PSU analysis covariate"
                else:
                    status="Missing in PSU analysis binary field; zero-coded source categories may not be distinguishable"
                complete=True
            # Never report an apparent zero as true source completeness when pre-imputation
            # missingness was not preserved in the final analytic MIMIC table.
            if not complete:
                a,b=None,None
            ca,pa=safe_count(a,n1)
            cb,pb=safe_count(b,n0)
            outrows.append({
                "dataset":source,
                "section":str(row["section"]),
                "variable":var,
                "characteristic":str(row["characteristic"]),
                "type":typ,
                "n_deescalation":n1,
                "deescalation_missing_n":ca,
                "deescalation_missing_percent":pa,
                "n_continued":n0,
                "continued_missing_n":cb,
                "continued_missing_percent":pb,
                "interpretation":status,
                "preimputation_status":"recoverable" if complete else "not_recoverable",
            })
        subset=[r for r in outrows if r["dataset"]==source]
        tally[source]={"variables":len(subset),
            "recoverable":sum(r["preimputation_status"]=="recoverable" for r in subset),
            "not_recoverable":sum(r["preimputation_status"]=="not_recoverable" for r in subset),
            "group_n_deescalation":n1,"group_n_continued":n0}
    OUT.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(outrows).to_csv(OUT/"target_groupwise_missingness.csv",index=False)
    report={
       "purpose":"TARGET statement item 11 descriptive missingness frequencies by observed treatment group",
       "sources":"Existing frozen baseline detailed aggregate outputs, no patient-level records read",
       "privacy":"Only aggregate counts; nonzero cells under 11 suppressed; unknown source missingness labeled NR",
       "model_or_outcome_reestimated":False,
       "denominators_fixed":True,
       "totals":tally,
       "note":"In MIMIC many source values were already zero-coded or imputed in the frozen cohort. Original measurement absence is only recoverable where pre-imputation flags were retained. For Penn State continuous covariates, missingness is before median imputation. Binary source-level missingness can be incompletely distinguishable from zero-coded clinical proxies.",
       "rows":outrows,
    }
    (OUT/"target_groupwise_missingness.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps({"totals":tally,"n_rows":len(outrows)},indent=2))

if __name__=="__main__":
    main()

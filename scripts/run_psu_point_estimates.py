#!/usr/bin/env python3
"""Run frozen PSU PS specification plus frozen outcomes and aggregate point estimates.

This wrapper deliberately reuses the already frozen PS diagnostic source verbatim and
injects outcome/effect aggregation only after PS fitting. No row-level or identifier
artifacts are exported. No bootstrap inference is performed in this task.
"""
from __future__ import annotations

import runpy
from pathlib import Path


INJECT = r'''
    # ---- frozen outcome construction and point-estimate block ----
    # This executes only after the frozen PS has been fit. It does not alter exposure,
    # covariates, imputation, model fitting, or the previously computed weights.
    con.execute("drop table if exists days_effects")
    con.execute("drop table if exists day_flags_effects")
    con.execute("drop table if exists day_triplets_effects")
    con.execute("drop table if exists abx_summary_effects")
    con.execute("drop table if exists outcomes_effects")

    any_systemic = f"(({broad}) or ({nonbroad}))"
    con.execute("""
      create temp table days_effects as
      select c.patid,c.encounterid,
             cast(c.anchor_ts+interval 96 hour as date) analysis_date,
             c.discharge_date,c.death_date,
             gs.day_idx::integer day_idx,
             cast(c.anchor_ts+interval 96 hour as date)+gs.day_idx::integer day_date
      from cohort c cross join generate_series(0,29) as gs(day_idx)
    """)
    con.execute(f"""
      create temp table day_flags_effects as
      select d.patid,d.encounterid,d.day_idx,d.day_date,
             max(case when {any_systemic} and {p_start_expr}<=d.day_date and {p_end_expr}>=d.day_date then 1 else 0 end) any_abx,
             max(case when {broad} and {p_start_expr}<=d.day_date and {p_end_expr}>=d.day_date then 1 else 0 end) broad_abx
      from days_effects d left join p
        on cast(p.{qi(pp)} as varchar)=d.patid
       and cast(p.{qi(pe)} as varchar)=d.encounterid
      group by 1,2,3,4
    """)
    con.execute("""
      create temp table day_triplets_effects as
      select *,
             lead(any_abx,1,0) over(partition by patid,encounterid order by day_idx) any_abx_d1,
             lead(any_abx,2,0) over(partition by patid,encounterid order by day_idx) any_abx_d2
      from day_flags_effects
    """)
    con.execute("""
      create temp table abx_summary_effects as
      select patid,encounterid,
             sum(any_abx)::integer antibiotic_days_30d,
             sum(broad_abx)::integer broad_antibiotic_days_30d,
             max(case when day_idx>=7 and day_idx<=27 and any_abx=1 and any_abx_d1=1 and any_abx_d2=1 then 1 else 0 end)::integer late_recurrent_or_persistent_abx_course_30d
      from day_triplets_effects group by 1,2
    """)
    con.execute("""
      create temp table outcomes_effects as
      select c.patid,c.encounterid,
        case when c.death_date is not null
               and c.death_date>=cast(c.anchor_ts+interval 96 hour as date)
               and c.death_date<=cast(c.anchor_ts+interval 96 hour as date)+30 then 1 else 0 end death_30d,
        case when c.death_date is not null
               and c.death_date>=cast(c.anchor_ts+interval 96 hour as date)
               and c.death_date<=cast(c.anchor_ts+interval 96 hour as date)+30 then 0.0
             else greatest(0.0,30.0-least(30.0,greatest(0.0,date_diff('day',cast(c.anchor_ts+interval 96 hour as date),c.discharge_date)::double))) end hospital_free_days_30d,
        a.antibiotic_days_30d,
        case when c.death_date is not null
               and c.death_date>=cast(c.anchor_ts+interval 96 hour as date)
               and c.death_date<=cast(c.anchor_ts+interval 96 hour as date)+30 then 0.0
             else greatest(0.0,30.0-a.antibiotic_days_30d::double) end antibiotic_free_days_30d,
        a.broad_antibiotic_days_30d,
        case when c.death_date is not null
               and c.death_date>=cast(c.anchor_ts+interval 96 hour as date)
               and c.death_date<=cast(c.anchor_ts+interval 96 hour as date)+30
             then greatest(0.1,least(30.0,date_diff('day',cast(c.anchor_ts+interval 96 hour as date),c.death_date)::double))
             else 30.0 end days_alive_30d,
        a.late_recurrent_or_persistent_abx_course_30d,
        least(30.0,greatest(0.0,date_diff('day',cast(c.anchor_ts+interval 96 hour as date),c.discharge_date)::double)) observable_hospital_days_after_landmark
      from cohort c join abx_summary_effects a using(patid,encounterid)
    """)
    con.execute("""
      alter table outcomes_effects add column normalized_antibiotic_exposure_30d double;
      update outcomes_effects set normalized_antibiotic_exposure_30d=least(1.0,greatest(0.0,antibiotic_days_30d/days_alive_30d));
      alter table outcomes_effects add column normalized_broad_antibiotic_exposure_30d double;
      update outcomes_effects set normalized_broad_antibiotic_exposure_30d=least(1.0,greatest(0.0,broad_antibiotic_days_30d/days_alive_30d));
    """)

    odf=con.execute("select * from outcomes_effects").fetchdf()
    df=df.merge(odf,on=["patid","encounterid"],how="left",validate="one_to_one")
    if len(df)!=n or int(df.death_30d.isna().sum())!=0:
        raise RuntimeError("Outcome merge parity failure")

    # Prespecified weighting estimands/sensitivities.
    ow=np.where(y==1,1.0-ps,ps)
    q01,q99=np.quantile(sw,[0.01,0.99]); tw_1_99=np.clip(sw,q01,q99)
    q025,q975=np.quantile(sw,[0.025,0.975]); tw_2p5_97p5=np.clip(sw,q025,q975)
    methods={
        "stabilized_ate_iptw":("ATE",sw),
        "overlap_weighting":("ATO",ow),
        "trunc_1_99":("ATE_truncated",tw_1_99),
        "trunc_2p5_97p5":("ATE_truncated",tw_2p5_97p5),
    }

    effect_outcomes=[
        ("death_30d","binary"),
        ("hospital_free_days_30d","continuous"),
        ("antibiotic_free_days_30d","continuous"),
        ("normalized_antibiotic_exposure_30d","continuous"),
        ("normalized_broad_antibiotic_exposure_30d","continuous"),
        ("late_recurrent_or_persistent_abx_course_30d","binary"),
    ]
    erows=[]
    for method,(estimand,w) in methods.items():
        for outcome,kind in effect_outcomes:
            z=pd.to_numeric(df[outcome],errors="coerce").to_numpy(float)
            m1=y==1; m0=y==0
            mu1=weighted_mean(z[m1],w[m1]); mu0=weighted_mean(z[m0],w[m0])
            erows.append({
                "method":method,"estimand":estimand,"outcome":outcome,"type":kind,
                "deescalated_mean_or_risk":float(mu1),"continued_mean_or_risk":float(mu0),
                "difference_A1_minus_A0":float(mu1-mu0),
                "risk_ratio_A1_over_A0":float(mu1/mu0) if kind=="binary" and mu0>0 else None,
            })
    effects=pd.DataFrame(erows)
    effects.to_csv(args.output_dir/"point_estimates.csv",index=False)

    wd=[]
    for method,(estimand,w) in methods.items():
        posts=[]
        for col in features:
            xv=X[col].to_numpy(float); posts.append(abs(smd(xv,y,w)))
        wd.append({
            "method":method,"estimand":estimand,
            "treated_ess":ess(w[y==1]),"continued_ess":ess(w[y==0]),
            "mean_weight":float(np.mean(w)),"p99_weight":float(np.quantile(w,.99)),"max_weight":float(np.max(w)),
            "max_abs_post_smd":float(np.nanmax(posts)),
        })
    wdiag=pd.DataFrame(wd)
    wdiag.to_csv(args.output_dir/"effect_weighting_diagnostics.csv",index=False)

    unadj=[]
    for outcome,kind in effect_outcomes:
        z=pd.to_numeric(df[outcome],errors="coerce").to_numpy(float)
        mu1=float(np.nanmean(z[y==1])); mu0=float(np.nanmean(z[y==0]))
        unadj.append({"outcome":outcome,"type":kind,"deescalated_mean_or_risk":mu1,"continued_mean_or_risk":mu0,"difference_A1_minus_A0":mu1-mu0,"risk_ratio_A1_over_A0":(mu1/mu0 if kind=="binary" and mu0>0 else None)})
    pd.DataFrame(unadj).to_csv(args.output_dir/"unadjusted_outcomes_by_exposure.csv",index=False)

    pdeath=effects[(effects.method=="stabilized_ate_iptw")&(effects.outcome=="death_30d")].iloc[0]
    summary_effects={
        "privacy_mode":"aggregate_only",
        "strict_cohort_n":safe(n),
        "deescalated_n":safe(int(y.sum())),
        "continued_n":safe(int((1-y).sum())),
        "frozen_ps_reproduced_max_abs_post_smd":float(wdiag.loc[wdiag.method=="stabilized_ate_iptw","max_abs_post_smd"].iloc[0]),
        "primary_method":"stabilized_ate_iptw",
        "primary_estimand":"ATE",
        "primary_outcome":"death_30d",
        "primary_deescalated_risk":float(pdeath.deescalated_mean_or_risk),
        "primary_continued_risk":float(pdeath.continued_mean_or_risk),
        "primary_risk_difference":float(pdeath.difference_A1_minus_A0),
        "primary_risk_ratio":float(pdeath.risk_ratio_A1_over_A0),
        "sensitivity_methods":["overlap_weighting","trunc_1_99","trunc_2p5_97p5"],
        "guardrail":"Point estimates only. No bootstrap confidence intervals. Exposure, covariates, PS, and outcomes were frozen before this task.",
    }
    (args.output_dir/"effect_summary.json").write_text(json.dumps(summary_effects,indent=2))
    # ---- end point-estimate block ----
'''


def main() -> None:
    here=Path(__file__).resolve().parent
    source_path=here/"audit_psu_ps_balance.py"
    source=source_path.read_text(encoding="utf-8")
    marker="    max_pre=float(bdf.abs_pre_smd.max()); max_post=float(bdf.abs_post_smd.max()); worst=str(bdf.iloc[0].variable)\n"
    if source.count(marker)!=1:
        raise RuntimeError("Frozen PS source marker missing or ambiguous")
    source=source.replace(marker,INJECT+"\n"+marker,1)
    ns={"__name__":"__main__","__file__":str(source_path)}
    exec(compile(source,str(source_path),"exec"),ns,ns)


if __name__=="__main__":
    main()

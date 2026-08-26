#!/usr/bin/env python3
"""Aggregate-only diagnostic for PSU lab/cohort linkage and clock alignment.

Explains why exact lab concepts are abundant globally yet have zero 0-72h coverage in the
strict modified PSU cohort. Exports only aggregate counts; no identifiers or row-level data.
"""
from __future__ import annotations
import argparse, ast, json, re
from pathlib import Path
import duckdb, pandas as pd

MIN_CELL=11
BROAD_PATTERN="vancomycin|vancocin|piperacillin|tazobactam|zosyn|cefepime|ceftazidime|fortaz|meropenem|merrem|imipenem|primaxin|aztreonam|azactam|linezolid|zyvox|daptomycin|cubicin|ceftolozane|zerbaxa|avibactam|avycaz"
NON_SYSTEMIC_PATTERN="cayston|inhal|nebul|tablet|capsule|oral solution|oral suspension|by mouth|\\bpo\\b"
LAB_CODES={"lactate":["2524-7","32693-4","19239-3"],"creatinine":["2160-0","38483-4"],"wbc":["6690-2"],"platelet":["777-3"],"bilirubin_total":["1975-2"]}

def q(s): return "'"+str(s).replace("'","''")+"'"
def qi(s): return '"'+s.replace('"','""')+'"'
def safe(n):
    if n is None:return None
    n=int(n); return n if n>=MIN_CELL else None

def find(root, stems):
    for stem in stems:
        c=sorted(root.glob(f"PCORnet/parquet/**/{stem}*.parquet")) or sorted(root.glob(f"PCORnet/parquet/**/{stem.replace('_','')}*.parquet"))
        if c:
            e=[p for p in c if p.stem.lower()==stem.lower()]
            return e[0] if e else max(c,key=lambda p:p.stat().st_size)
    raise FileNotFoundError(stems)

def first(cols,names):
    lut={c.upper():c for c in cols}
    for n in names:
        if n.upper() in lut:return lut[n.upper()]
    return None

def parse_codes(path):
    tree=ast.parse(path.read_text(errors='ignore')); inc=set(); exc=set()
    for node in tree.body:
        if not isinstance(node,ast.Assign):continue
        names=[t.id.lower() for t in node.targets if isinstance(t,ast.Name)]
        try:v=ast.literal_eval(node.value)
        except Exception:continue
        for name in names:
            if 'exclude' in name and isinstance(v,(set,list,tuple)):
                for x in v:
                    x=str(x).strip()
                    if re.fullmatch(r'\d{4,9}',x):exc.add(x)
            elif isinstance(v,(list,tuple)):
                for x in v:
                    if isinstance(x,(list,tuple)) and x:
                        z=str(x[0]).strip()
                        if re.fullmatch(r'\d{4,9}',z):inc.add(z)
    return inc-exc,exc

def ts(alias,date,time):
    if not date:return None
    if time:return f"try_cast(cast({alias}.{qi(date)} as varchar)||' '||coalesce(cast({alias}.{qi(time)} as varchar),'00:00:00') as timestamp)"
    return f"try_cast({alias}.{qi(date)} as timestamp)"

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('data_root',type=Path); ap.add_argument('--output-dir',type=Path,required=True); a=ap.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True)
    root=a.data_root
    paths={'s':find(root,['sepsis_encounter']),'p':find(root,['prescribing']),'d':find(root,['death']),'l':find(root,['lab_reduced','lab_result_cm'])}
    con=duckdb.connect(); con.execute('pragma threads=4')
    for k,p in paths.items():con.execute(f"create view {k} as select * from read_parquet({q(str(p))})")
    cols={k:set(con.execute(f'describe {k}').fetchdf()['column_name'].astype(str)) for k in paths}
    sp,se,ad,dc=first(cols['s'],['PATID']),first(cols['s'],['ENCOUNTERID']),first(cols['s'],['ADMIT_DATE']),first(cols['s'],['DISCHARGE_DATE'])
    pp,pe=first(cols['p'],['PATID']),first(cols['p'],['ENCOUNTERID']); od,ot=first(cols['p'],['RX_ORDER_DATE']),first(cols['p'],['RX_ORDER_TIME']); pn,pc,pr=first(cols['p'],['RAW_RX_MED_NAME']),first(cols['p'],['RXNORM_CUI']),first(cols['p'],['RX_ROUTE'])
    dp,dd=first(cols['d'],['PATID']),first(cols['d'],['DEATH_DATE'])
    inc,exc=parse_codes(root/'PCORnet/code/config/codes_antibiotics.py'); incsql=','.join(q(x) for x in sorted(inc)) or "''"; excsql=','.join(q(x) for x in sorted(exc)) or "''"
    pts=ts('p',od,ot); name=f"lower(coalesce(cast(p.{qi(pn)} as varchar),''))"; code=f"trim(coalesce(cast(p.{qi(pc)} as varchar),''))"; route=f"upper(trim(coalesce(cast(p.{qi(pr)} as varchar),'')))"
    broad=f"(({code} in ({incsql}) or regexp_matches({name},{q(BROAD_PATTERN)})) and {code} not in ({excsql}) and not regexp_matches({name},{q(NON_SYSTEMIC_PATTERN)}) and {route} not in ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    con.execute(f"create temp table base as select distinct cast({qi(sp)} as varchar) patid,cast({qi(se)} as varchar) encounterid,try_cast({qi(ad)} as date) admit_date,try_cast({qi(dc)} as date) discharge_date from s")
    con.execute(f"create temp table anchors as select cast(p.{qi(pp)} as varchar) patid,cast(p.{qi(pe)} as varchar) encounterid,min({pts}) anchor_ts from p join base b on cast(p.{qi(pp)} as varchar)=b.patid and cast(p.{qi(pe)} as varchar)=b.encounterid where {broad} and {pts}>=cast(b.admit_date as timestamp) and {pts}<cast(b.admit_date as timestamp)+interval 24 hour group by 1,2")
    con.execute(f"create temp table deaths as select cast({qi(dp)} as varchar) patid,min(try_cast({qi(dd)} as date)) death_date from d group by 1")
    con.execute("create temp table cohort as select a.*,b.admit_date,b.discharge_date,d.death_date from anchors a join base b using(patid,encounterid) left join deaths d using(patid) where (b.discharge_date is null or b.discharge_date>cast(a.anchor_ts+interval 96 hour as date)) and (d.death_date is null or d.death_date>cast(a.anchor_ts+interval 96 hour as date))")
    n=int(con.execute('select count(*) from cohort').fetchone()[0])

    lp,le=first(cols['l'],['PATID']),first(cols['l'],['ENCOUNTERID']); lc=first(cols['l'],['LAB_LOINC']); lv=first(cols['l'],['RESULT_NUM']);
    clocks=[]
    for label,ds,tsc in [('specimen',['SPECIMEN_DATE'],['SPECIMEN_TIME']),('result',['RESULT_DATE'],['RESULT_TIME']),('order',['LAB_ORDER_DATE'],['LAB_ORDER_TIME'])]:
        dcol=first(cols['l'],ds); tcol=first(cols['l'],tsc); expr=ts('l',dcol,tcol) if dcol else None
        if expr:clocks.append((label,expr,dcol,tcol))

    join_rows=[]
    for join_name,join_cond in [
        ('patid_encounter',f"cast(l.{qi(lp)} as varchar)=c.patid and cast(l.{qi(le)} as varchar)=c.encounterid"),
        ('encounter_only',f"cast(l.{qi(le)} as varchar)=c.encounterid")]:
        any_rows=int(con.execute(f"select count(distinct c.encounterid) from cohort c join l on {join_cond}").fetchone()[0])
        numeric=int(con.execute(f"select count(distinct c.encounterid) from cohort c join l on {join_cond} where try_cast(l.{qi(lv)} as double) is not null").fetchone()[0])
        exact_codes=','.join(q(x) for xs in LAB_CODES.values() for x in xs)
        mapped=int(con.execute(f"select count(distinct c.encounterid) from cohort c join l on {join_cond} where cast(l.{qi(lc)} as varchar) in ({exact_codes}) and try_cast(l.{qi(lv)} as double) is not null").fetchone()[0])
        join_rows.append({'join_rule':join_name,'any_lab':safe(any_rows),'numeric_lab':safe(numeric),'mapped_core_lab':safe(mapped),'cohort_n':safe(n)})
    pd.DataFrame(join_rows).to_csv(a.output_dir/'lab_join_diagnostic.csv',index=False)

    comp=[]
    for label,expr,dcol,tcol in clocks:
        present=int(con.execute(f"select count(distinct c.encounterid) from cohort c join l on cast(l.{qi(lp)} as varchar)=c.patid and cast(l.{qi(le)} as varchar)=c.encounterid where {expr} is not null").fetchone()[0])
        comp.append({'clock':label,'date_field':dcol,'time_field':tcol or '','cohort_encounters_with_timestamp':safe(present),'cohort_n':safe(n)})
    pd.DataFrame(comp).to_csv(a.output_dir/'lab_clock_completeness.csv',index=False)

    bins=[]
    for label,expr,_,_ in clocks:
        for concept,codes in LAB_CODES.items():
            csql=','.join(q(x) for x in codes)
            sql=f"""select
              sum(case when {expr}<c.anchor_ts-interval 30 day then 1 else 0 end) b_lt_m30,
              sum(case when {expr}>=c.anchor_ts-interval 30 day and {expr}<c.anchor_ts-interval 1 day then 1 else 0 end) b_m30_m1,
              sum(case when {expr}>=c.anchor_ts-interval 1 day and {expr}<c.anchor_ts then 1 else 0 end) b_m1_0,
              sum(case when {expr}>=c.anchor_ts and {expr}<c.anchor_ts+interval 1 day then 1 else 0 end) b_0_1,
              sum(case when {expr}>=c.anchor_ts+interval 1 day and {expr}<c.anchor_ts+interval 3 day then 1 else 0 end) b_1_3,
              sum(case when {expr}>=c.anchor_ts+interval 3 day and {expr}<c.anchor_ts+interval 7 day then 1 else 0 end) b_3_7,
              sum(case when {expr}>=c.anchor_ts+interval 7 day then 1 else 0 end) b_ge7
              from cohort c join l on cast(l.{qi(lp)} as varchar)=c.patid and cast(l.{qi(le)} as varchar)=c.encounterid
              where cast(l.{qi(lc)} as varchar) in ({csql}) and try_cast(l.{qi(lv)} as double) is not null and {expr} is not null"""
            r=con.execute(sql).fetchone()
            bins.append({'clock':label,'concept':concept,'lt_minus30d':safe(r[0]),'minus30_to_minus1d':safe(r[1]),'minus1_to_0d':safe(r[2]),'day0_to1':safe(r[3]),'day1_to3':safe(r[4]),'day3_to7':safe(r[5]),'ge7d':safe(r[6])})
    pd.DataFrame(bins).to_csv(a.output_dir/'lab_anchor_delta_bins.csv',index=False)

    overlap=[]
    for concept,codes in LAB_CODES.items():
        csql=','.join(q(x) for x in codes)
        anyc=int(con.execute(f"select count(distinct c.encounterid) from cohort c join l on cast(l.{qi(lp)} as varchar)=c.patid and cast(l.{qi(le)} as varchar)=c.encounterid where cast(l.{qi(lc)} as varchar) in ({csql}) and try_cast(l.{qi(lv)} as double) is not null").fetchone()[0])
        overlap.append({'concept':concept,'cohort_encounters_anytime':safe(anyc),'cohort_n':safe(n),'proportion_anytime':anyc/n if n else None})
    pd.DataFrame(overlap).to_csv(a.output_dir/'lab_concept_anytime_coverage.csv',index=False)
    summary={'privacy_mode':'aggregate_only','strict_cohort_n':safe(n),'purpose':'diagnose zero timed core-lab coverage before PS modeling','guardrail':'No propensity score or treatment effects fit.'}
    (a.output_dir/'summary.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__':main()

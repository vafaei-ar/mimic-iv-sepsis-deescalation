#!/usr/bin/env python3
"""Aggregate-only PSU MED_ADMIN clock audit for vasopressor overlap.

Rebuilds the established strict modified PSU cohort, inspects raw MEDADMIN start/stop
TIME encodings, compares standard timestamp parsing with numeric-seconds-since-midnight
parsing, and reports vasopressor overlap counts for 0-24 h and 48-72 h. No identifiers
or row-level data are exported.
"""
from __future__ import annotations
import argparse, ast, json, re
from pathlib import Path
import duckdb, pandas as pd

MIN_CELL=11
BROAD_PATTERN="vancomycin|vancocin|piperacillin|tazobactam|zosyn|cefepime|ceftazidime|fortaz|meropenem|merrem|imipenem|primaxin|aztreonam|azactam|linezolid|zyvox|daptomycin|cubicin|ceftolozane|zerbaxa|avibactam|avycaz"
NON_SYSTEMIC_PATTERN="cayston|inhal|nebul|tablet|capsule|oral solution|oral suspension|by mouth|\\bpo\\b"
VASO_TERMS=["norepinephrine","levophed","phenylephrine","vasopressin","epinephrine","dopamine"]
VASO_EXCLUDE=["racepinephrine","nasal","ophthalm","lidocaine","topical"]

def q(v): return "'"+str(v).replace("'","''")+"'"
def qi(v): return '"'+str(v).replace('"','""')+'"'
def safe(n):
    if n is None:return None
    n=int(n); return n if n>=MIN_CELL else None

def find(root,stems):
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
                exc|={str(x).strip() for x in v if re.fullmatch(r'\d{4,9}',str(x).strip())}
            elif isinstance(v,(list,tuple)):
                for x in v:
                    if isinstance(x,(list,tuple)) and x and re.fullmatch(r'\d{4,9}',str(x[0]).strip()):inc.add(str(x[0]).strip())
    return inc-exc,exc

def ts(alias,date,time,mode='auto'):
    d=f"try_cast({alias}.{qi(date)} as date)"
    if not time:return f"cast({d} as timestamp)"
    raw=f"trim(cast({alias}.{qi(time)} as varchar))"; num=f"try_cast({raw} as double)"
    if mode=='text': return f"try_cast(cast({alias}.{qi(date)} as varchar)||' '||coalesce(cast({alias}.{qi(time)} as varchar),'00:00:00') as timestamp)"
    return f"case when {d} is null then null when strpos({raw},':')>0 then try_cast(cast({alias}.{qi(date)} as varchar)||' '||{raw} as timestamp) when {num} between 0 and 86399 then cast({d} as timestamp)+({num}*interval 1 second) else null end"

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('data_root',type=Path); ap.add_argument('--output-dir',type=Path,required=True); a=ap.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True)
    root=a.data_root; paths={'s':find(root,['sepsis_encounter']),'p':find(root,['prescribing']),'d':find(root,['death']),'m':find(root,['med_admin'])}
    con=duckdb.connect(); con.execute('pragma threads=4')
    for k,p in paths.items():con.execute(f"create view {k} as select * from read_parquet({q(str(p))})")
    cols={k:set(con.execute(f'describe {k}').fetchdf()['column_name'].astype(str)) for k in paths}
    sp,se,ad,dc=first(cols['s'],['PATID']),first(cols['s'],['ENCOUNTERID']),first(cols['s'],['ADMIT_DATE']),first(cols['s'],['DISCHARGE_DATE'])
    pp,pe=first(cols['p'],['PATID']),first(cols['p'],['ENCOUNTERID']); od,ot=first(cols['p'],['RX_ORDER_DATE']),first(cols['p'],['RX_ORDER_TIME']); pn,pc,pr=first(cols['p'],['RAW_RX_MED_NAME']),first(cols['p'],['RXNORM_CUI']),first(cols['p'],['RX_ROUTE']); dp,dd=first(cols['d'],['PATID']),first(cols['d'],['DEATH_DATE'])
    inc,exc=parse_codes(root/'PCORnet/code/config/codes_antibiotics.py'); incsql=','.join(q(x) for x in sorted(inc)) or "''"; excsql=','.join(q(x) for x in sorted(exc)) or "''"
    pts=ts('p',od,ot); name=f"lower(coalesce(cast(p.{qi(pn)} as varchar),''))"; code=f"trim(coalesce(cast(p.{qi(pc)} as varchar),''))"; route=f"upper(trim(coalesce(cast(p.{qi(pr)} as varchar),'')))"; broad=f"(({code} in ({incsql}) or regexp_matches({name},{q(BROAD_PATTERN)})) and {code} not in ({excsql}) and not regexp_matches({name},{q(NON_SYSTEMIC_PATTERN)}) and {route} not in ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    con.execute(f"create temp table base as select distinct cast({qi(sp)} as varchar) patid,cast({qi(se)} as varchar) encounterid,try_cast({qi(ad)} as date) admit_date,try_cast({qi(dc)} as date) discharge_date from s")
    con.execute(f"create temp table anchors as select cast(p.{qi(pp)} as varchar) patid,cast(p.{qi(pe)} as varchar) encounterid,min({pts}) anchor_ts from p join base b on cast(p.{qi(pp)} as varchar)=b.patid and cast(p.{qi(pe)} as varchar)=b.encounterid where {broad} and {pts}>=cast(b.admit_date as timestamp) and {pts}<cast(b.admit_date as timestamp)+interval 24 hour group by 1,2")
    con.execute(f"create temp table deaths as select cast({qi(dp)} as varchar) patid,min(try_cast({qi(dd)} as date)) death_date from d group by 1")
    con.execute("create temp table cohort as select a.*,b.discharge_date,d.death_date from anchors a join base b using(patid,encounterid) left join deaths d using(patid) where (b.discharge_date is null or b.discharge_date>cast(a.anchor_ts+interval 96 hour as date)) and (d.death_date is null or d.death_date>cast(a.anchor_ts+interval 96 hour as date))")
    n=int(con.execute('select count(*) from cohort').fetchone()[0])
    mp,me=first(cols['m'],['PATID']),first(cols['m'],['ENCOUNTERID']); sd,st=first(cols['m'],['MEDADMIN_START_DATE']),first(cols['m'],['MEDADMIN_START_TIME']); ed,et=first(cols['m'],['MEDADMIN_STOP_DATE']),first(cols['m'],['MEDADMIN_STOP_TIME']); mn=first(cols['m'],['RAW_MEDADMIN_MED_NAME'])
    if not all([mp,me,sd,mn]): raise RuntimeError('MED_ADMIN fields missing')
    mname=f"lower(coalesce(cast(m.{qi(mn)} as varchar),''))"; incl=' or '.join(f"strpos({mname},{q(x)})>0" for x in VASO_TERMS); excl=' or '.join(f"strpos({mname},{q(x)})>0" for x in VASO_EXCLUDE); vaso=f"({incl}) and not ({excl})"
    enc=[]
    for label,tcol in [('start',st),('stop',et)]:
        if not tcol: continue
        raw=f"trim(cast(m.{qi(tcol)} as varchar))"; num=f"try_cast({raw} as double)"
        r=con.execute(f"select count(*) n,count(*) filter(where strpos({raw},':')>0) colon_rows,count(*) filter(where {num} is not null) numeric_rows,min({num}),quantile_cont({num},0.5),quantile_cont({num},0.95),max({num}) from m where {vaso} and m.{qi(tcol)} is not null").fetchone()
        enc.append({'clock':label,'rows':safe(r[0]),'colon_rows':safe(r[1]),'numeric_rows':safe(r[2]),'numeric_min':r[3],'numeric_median':r[4],'numeric_p95':r[5],'numeric_max':r[6]})
    pd.DataFrame(enc).to_csv(a.output_dir/'time_encoding.csv',index=False)
    rows=[]
    for mode in ['text','auto']:
        start=ts('m',sd,st,mode); stop=ts('m',ed,et,mode) if ed else start
        stop=f"coalesce({stop},{start}+interval 1 hour)"
        for win,lo,hi in [('0_24h',0,24),('48_72h',48,72),('0_72h',0,72)]:
            x=int(con.execute(f"select count(distinct c.encounterid) from cohort c join m on cast(m.{qi(mp)} as varchar)=c.patid and cast(m.{qi(me)} as varchar)=c.encounterid where {vaso} and {start}<c.anchor_ts+interval {hi} hour and {stop}>c.anchor_ts+interval {lo} hour").fetchone()[0])
            rows.append({'parser':mode,'window':win,'encounters':safe(x),'cohort_n':safe(n),'proportion':x/n})
    pd.DataFrame(rows).to_csv(a.output_dir/'vasopressor_overlap_by_parser.csv',index=False)
    summary={'privacy_mode':'aggregate_only','strict_cohort_n':safe(n),'purpose':'resolve zero vasopressor overlap before final covariate freeze','guardrail':'No propensity score or treatment effects fit.'}
    (a.output_dir/'summary.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__': main()

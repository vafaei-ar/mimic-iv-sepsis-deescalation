#!/usr/bin/env python3
"""Aggregate-only PSU laboratory clock semantics audit.

Compares LAB_ORDER_DATE/TIME, SPECIMEN_DATE/TIME, and RESULT_DATE/TIME for core
laboratory rows linked to the frozen strict PSU cohort. Raw specimen/result TIME fields
are decoded as seconds since midnight after the encoding audit demonstrated values from
0 to 86340 with all numeric rows valid under that representation. Exports aggregate
summaries only. No identifiers or patient-level rows are exported.
"""
from __future__ import annotations
import argparse, ast, json, re
from pathlib import Path
import duckdb, pandas as pd

MIN_CELL = 11
BROAD_PATTERN = "vancomycin|vancocin|piperacillin|tazobactam|zosyn|cefepime|ceftazidime|fortaz|meropenem|merrem|imipenem|primaxin|aztreonam|azactam|linezolid|zyvox|daptomycin|cubicin|ceftolozane|zerbaxa|avibactam|avycaz"
NON_SYSTEMIC_PATTERN = "cayston|inhal|nebul|tablet|capsule|oral solution|oral suspension|by mouth|\\bpo\\b"
LAB_CODES = {
    "lactate": ["2524-7", "32693-4", "19239-3"],
    "creatinine": ["2160-0", "38483-4"],
    "wbc": ["6690-2"],
    "platelet": ["777-3"],
    "bilirubin_total": ["1975-2"],
}

def q(v): return "'" + str(v).replace("'", "''") + "'"
def qi(v): return '"' + str(v).replace('"','""') + '"'
def safe(n):
    if n is None: return None
    n = int(n)
    return n if n >= MIN_CELL else None

def find(root, stems):
    for stem in stems:
        c = sorted(root.glob(f"PCORnet/parquet/**/{stem}*.parquet")) or sorted(root.glob(f"PCORnet/parquet/**/{stem.replace('_','')}*.parquet"))
        if c:
            e = [p for p in c if p.stem.lower() == stem.lower()]
            return e[0] if e else max(c, key=lambda p:p.stat().st_size)
    raise FileNotFoundError(stems)

def first(cols, names):
    lut = {c.upper(): c for c in cols}
    for n in names:
        if n.upper() in lut: return lut[n.upper()]
    return None

def parse_codes(path):
    tree = ast.parse(path.read_text(errors='ignore')); inc=set(); exc=set()
    for node in tree.body:
        if not isinstance(node, ast.Assign): continue
        names=[t.id.lower() for t in node.targets if isinstance(t,ast.Name)]
        try: v=ast.literal_eval(node.value)
        except Exception: continue
        for name in names:
            if 'exclude' in name and isinstance(v,(set,list,tuple)):
                exc |= {str(x).strip() for x in v if re.fullmatch(r'\d{4,9}',str(x).strip())}
            elif isinstance(v,(list,tuple)):
                for x in v:
                    if isinstance(x,(list,tuple)) and x and re.fullmatch(r'\d{4,9}',str(x[0]).strip()): inc.add(str(x[0]).strip())
    return inc-exc, exc

def ts_standard(alias,date,time):
    if not date: return None
    if time:
        return f"try_cast(cast({alias}.{qi(date)} as varchar)||' '||coalesce(cast({alias}.{qi(time)} as varchar),'00:00:00') as timestamp)"
    return f"try_cast({alias}.{qi(date)} as timestamp)"

def ts_seconds_since_midnight(alias,date,time):
    if not date: return None
    if not time: return f"try_cast({alias}.{qi(date)} as timestamp)"
    num=f"try_cast({alias}.{qi(time)} as double)"
    return f"case when try_cast({alias}.{qi(date)} as date) is not null and {num} between 0 and 86399 then cast(try_cast({alias}.{qi(date)} as date) as timestamp) + ({num} * interval '1 second') else NULL end"

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('data_root',type=Path); ap.add_argument('--output-dir',type=Path,required=True); a=ap.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True)
    root=a.data_root
    paths={'s':find(root,['sepsis_encounter']),'p':find(root,['prescribing']),'d':find(root,['death']),'l':find(root,['lab_reduced','lab_result_cm'])}
    con=duckdb.connect(); con.execute('pragma threads=4')
    for k,p in paths.items(): con.execute(f"create view {k} as select * from read_parquet({q(str(p))})")
    cols={k:set(con.execute(f'describe {k}').fetchdf()['column_name'].astype(str)) for k in paths}

    sp,se,ad,dc=first(cols['s'],['PATID']),first(cols['s'],['ENCOUNTERID']),first(cols['s'],['ADMIT_DATE']),first(cols['s'],['DISCHARGE_DATE'])
    pp,pe=first(cols['p'],['PATID']),first(cols['p'],['ENCOUNTERID']); od,ot=first(cols['p'],['RX_ORDER_DATE']),first(cols['p'],['RX_ORDER_TIME']); pn,pc,pr=first(cols['p'],['RAW_RX_MED_NAME']),first(cols['p'],['RXNORM_CUI']),first(cols['p'],['RX_ROUTE'])
    dp,dd=first(cols['d'],['PATID']),first(cols['d'],['DEATH_DATE'])
    inc,exc=parse_codes(root/'PCORnet/code/config/codes_antibiotics.py'); incsql=','.join(q(x) for x in sorted(inc)) or "''"; excsql=','.join(q(x) for x in sorted(exc)) or "''"
    pts=ts_standard('p',od,ot); name=f"lower(coalesce(cast(p.{qi(pn)} as varchar),''))"; code=f"trim(coalesce(cast(p.{qi(pc)} as varchar),''))"; route=f"upper(trim(coalesce(cast(p.{qi(pr)} as varchar),'')))"
    broad=f"(({code} in ({incsql}) or regexp_matches({name},{q(BROAD_PATTERN)})) and {code} not in ({excsql}) and not regexp_matches({name},{q(NON_SYSTEMIC_PATTERN)}) and {route} not in ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    con.execute(f"create temp table base as select distinct cast({qi(sp)} as varchar) patid,cast({qi(se)} as varchar) encounterid,try_cast({qi(ad)} as date) admit_date,try_cast({qi(dc)} as date) discharge_date from s")
    con.execute(f"create temp table anchors as select cast(p.{qi(pp)} as varchar) patid,cast(p.{qi(pe)} as varchar) encounterid,min({pts}) anchor_ts from p join base b on cast(p.{qi(pp)} as varchar)=b.patid and cast(p.{qi(pe)} as varchar)=b.encounterid where {broad} and {pts}>=cast(b.admit_date as timestamp) and {pts}<cast(b.admit_date as timestamp)+interval 24 hour group by 1,2")
    con.execute(f"create temp table deaths as select cast({qi(dp)} as varchar) patid,min(try_cast({qi(dd)} as date)) death_date from d group by 1")
    con.execute("create temp table cohort as select a.*,b.admit_date,b.discharge_date,d.death_date from anchors a join base b using(patid,encounterid) left join deaths d using(patid) where (b.discharge_date is null or b.discharge_date>cast(a.anchor_ts+interval 96 hour as date)) and (d.death_date is null or d.death_date>cast(a.anchor_ts+interval 96 hour as date))")
    n=int(con.execute('select count(*) from cohort').fetchone()[0])

    lp,le=first(cols['l'],['PATID']),first(cols['l'],['ENCOUNTERID']); lc=first(cols['l'],['LAB_LOINC']); lv=first(cols['l'],['RESULT_NUM'])
    clock_fields={}
    clocks={}
    for label,dnames,tnames in [('order',['LAB_ORDER_DATE'],['LAB_ORDER_TIME']),('specimen',['SPECIMEN_DATE'],['SPECIMEN_TIME']),('result',['RESULT_DATE'],['RESULT_TIME'])]:
        dcol=first(cols['l'],dnames); tcol=first(cols['l'],tnames)
        if dcol:
            if label in ('specimen','result'):
                clocks[label]=ts_seconds_since_midnight('l',dcol,tcol)
            else:
                clocks[label]=ts_standard('l',dcol,tcol)
            clock_fields[label]=(dcol,tcol)
    allcodes=','.join(q(x) for xs in LAB_CODES.values() for x in xs)
    basecond=f"cast(l.{qi(lp)} as varchar)=c.patid and cast(l.{qi(le)} as varchar)=c.encounterid and cast(l.{qi(lc)} as varchar) in ({allcodes}) and try_cast(l.{qi(lv)} as double) is not null"

    complete=[]
    encoding={}
    for label,expr in clocks.items():
        cnt=int(con.execute(f"select count(*) from cohort c join l on {basecond} where {expr} is not null").fetchone()[0])
        enc=int(con.execute(f"select count(distinct c.encounterid) from cohort c join l on {basecond} where {expr} is not null").fetchone()[0])
        dcol,tcol=clock_fields[label]
        raw={}
        if tcol:
            t=f"trim(cast(l.{qi(tcol)} as varchar))"; num=f"try_cast({t} as double)"; d=f"try_cast(l.{qi(dcol)} as date)"
            raw_nonnull=int(con.execute(f"select count(*) from cohort c join l on {basecond} where l.{qi(tcol)} is not null and {t}<>''").fetchone()[0])
            colon=int(con.execute(f"select count(*) from cohort c join l on {basecond} where l.{qi(tcol)} is not null and strpos({t},':')>0").fetchone()[0])
            numeric=int(con.execute(f"select count(*) from cohort c join l on {basecond} where {num} is not null").fetchone()[0])
            stats=con.execute(f"select min({num}),quantile_cont({num},0.05),quantile_cont({num},0.5),quantile_cont({num},0.95),max({num}) from cohort c join l on {basecond} where {num} is not null").fetchone()
            valid_min=f"{num} between 0 and 1439"
            valid_sec=f"{num} between 0 and 86399"
            hh=f"floor({num}/100)"; mm=f"mod({num},100)"; valid_hhmm=f"{num} between 0 and 2359 and {hh} between 0 and 23 and {mm} between 0 and 59"
            parse_min=int(con.execute(f"select count(*) from cohort c join l on {basecond} where {d} is not null and {valid_min}").fetchone()[0])
            parse_sec=int(con.execute(f"select count(*) from cohort c join l on {basecond} where {d} is not null and {valid_sec}").fetchone()[0])
            parse_hhmm=int(con.execute(f"select count(*) from cohort c join l on {basecond} where {d} is not null and {valid_hhmm}").fetchone()[0])
            raw={'raw_time_nonnull_rows':safe(raw_nonnull),'colon_string_rows':safe(colon),'numeric_rows':safe(numeric),'numeric_min':stats[0],'numeric_p05':stats[1],'numeric_median':stats[2],'numeric_p95':stats[3],'numeric_max':stats[4],'rows_valid_minutes_since_midnight':safe(parse_min),'rows_valid_seconds_since_midnight':safe(parse_sec),'rows_valid_hhmm_numeric':safe(parse_hhmm)}
            encoding[label]=raw
        complete.append({'clock':label,'rows_with_clock':safe(cnt),'encounters_with_clock':safe(enc),'cohort_n':safe(n),**raw})
    pd.DataFrame(complete).to_csv(a.output_dir/'clock_completeness.csv',index=False)

    pairrows=[]
    labels=list(clocks)
    for i,a1 in enumerate(labels):
        for a2 in labels[i+1:]:
            e1,e2=clocks[a1],clocks[a2]
            sql=f"select count(*) n, count(distinct c.encounterid) enc, avg(abs(date_diff('minute',{e1},{e2}))) mean_abs_min, median(abs(date_diff('minute',{e1},{e2}))) median_abs_min, quantile_cont(abs(date_diff('minute',{e1},{e2})),0.95) p95_abs_min, sum(case when cast({e1} as date)=cast({e2} as date) then 1 else 0 end) same_date from cohort c join l on {basecond} where {e1} is not null and {e2} is not null"
            r=con.execute(sql).fetchone(); nrows=int(r[0]) if r[0] is not None else 0
            pairrows.append({'clock_a':a1,'clock_b':a2,'paired_rows':safe(nrows),'paired_encounters':safe(r[1]),'mean_abs_minutes':float(r[2]) if r[2] is not None else None,'median_abs_minutes':float(r[3]) if r[3] is not None else None,'p95_abs_minutes':float(r[4]) if r[4] is not None else None,'same_date_rows':safe(r[5]),'same_date_proportion':(float(r[5])/nrows) if nrows else None})
    pd.DataFrame(pairrows).to_csv(a.output_dir/'clock_pair_offsets.csv',index=False)

    anch=[]
    for label,expr in clocks.items():
        r=con.execute(f"select count(*) n, median(date_diff('hour',c.anchor_ts,{expr})) med_h, quantile_cont(date_diff('hour',c.anchor_ts,{expr}),0.05) p05_h, quantile_cont(date_diff('hour',c.anchor_ts,{expr}),0.95) p95_h from cohort c join l on {basecond} where {expr} is not null").fetchone()
        anch.append({'clock':label,'rows':safe(r[0]),'median_hours_from_anchor':float(r[1]) if r[1] is not None else None,'p05_hours_from_anchor':float(r[2]) if r[2] is not None else None,'p95_hours_from_anchor':float(r[3]) if r[3] is not None else None})
    pd.DataFrame(anch).to_csv(a.output_dir/'clock_anchor_offsets.csv',index=False)

    concept=[]
    for clock_name,clock_expr in clocks.items():
        for concept_name,codes in LAB_CODES.items():
            csql=','.join(q(x) for x in codes)
            for win,lo,hi in [('0_24h',0,24),('48_72h',48,72),('pre72',0,72)]:
                x=int(con.execute(f"select count(distinct c.encounterid) from cohort c join l on cast(l.{qi(lp)} as varchar)=c.patid and cast(l.{qi(le)} as varchar)=c.encounterid where cast(l.{qi(lc)} as varchar) in ({csql}) and try_cast(l.{qi(lv)} as double) is not null and {clock_expr}>=c.anchor_ts+interval {lo} hour and {clock_expr}<c.anchor_ts+interval {hi} hour").fetchone()[0])
                concept.append({'clock':clock_name,'concept':concept_name,'window':win,'count':safe(x),'cohort_n':safe(n),'proportion':x/n if n else None})
    pd.DataFrame(concept).to_csv(a.output_dir/'order_clock_core_lab_coverage.csv',index=False)

    summary={'privacy_mode':'aggregate_only','strict_cohort_n':safe(n),'available_clocks':list(clocks),'specimen_result_time_decoding':'numeric seconds since midnight','candidate_primary_lab_clock':'pending comparison of corrected specimen/result clocks with order clock','guardrail':'No propensity score or treatment effects fit.'}
    (a.output_dir/'summary.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__': main()

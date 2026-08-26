#!/usr/bin/env python3
"""Aggregate-only diagnostic for PSU core-lab timing and MAP-code resolution.

Explains why prespecified core-lab LOINCs had zero strict-cohort time-window coverage in
covariate_mapping_freeze_v2, and searches for an aggregate MAP candidate in OBS_CLIN.
No identifiers, row-level patient data, or result free text are exported.
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

def q(s): return "'" + str(s).replace("'", "''") + "'"
def qi(s): return '"' + s.replace('"', '""') + '"'
def safe(n):
    if n is None: return None
    n = int(n)
    return n if n >= MIN_CELL else None

def find(root, stems):
    for stem in stems:
        c = sorted(root.glob(f"PCORnet/parquet/**/{stem}*.parquet")) or sorted(root.glob(f"PCORnet/parquet/**/{stem.replace('_','')}*.parquet"))
        if c:
            e = [p for p in c if p.stem.lower() == stem.lower()]
            return e[0] if e else max(c, key=lambda p: p.stat().st_size)
    raise FileNotFoundError(stems)

def first(cols, names):
    lut = {c.upper(): c for c in cols}
    for n in names:
        if n.upper() in lut: return lut[n.upper()]
    return None

def parse_codes(path):
    tree = ast.parse(path.read_text(errors="ignore")); inc=set(); exc=set()
    for node in tree.body:
        if not isinstance(node, ast.Assign): continue
        names=[t.id.lower() for t in node.targets if isinstance(t, ast.Name)]
        try: v=ast.literal_eval(node.value)
        except Exception: continue
        for name in names:
            if "exclude" in name and isinstance(v,(set,list,tuple)):
                for x in v:
                    x=str(x).strip()
                    if re.fullmatch(r"\d{4,9}",x): exc.add(x)
            elif isinstance(v,(list,tuple)):
                for x in v:
                    if isinstance(x,(list,tuple)) and x:
                        z=str(x[0]).strip()
                        if re.fullmatch(r"\d{4,9}",z): inc.add(z)
    return inc-exc, exc

def ts(alias, date_col, time_col):
    if not date_col: return "NULL::TIMESTAMP"
    if time_col:
        return f"try_cast(cast({alias}.{qi(date_col)} as varchar)||' '||coalesce(cast({alias}.{qi(time_col)} as varchar),'00:00:00') as timestamp)"
    return f"try_cast({alias}.{qi(date_col)} as timestamp)"

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("data_root",type=Path); ap.add_argument("--output-dir",type=Path,required=True)
    a=ap.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True); root=a.data_root
    paths={"s":find(root,["sepsis_encounter"]),"p":find(root,["prescribing"]),"d":find(root,["death"]),"l":find(root,["lab_reduced","lab_result_cm"]),"o":find(root,["obs_clin"])}
    con=duckdb.connect(); con.execute("PRAGMA threads=4")
    for k,p in paths.items(): con.execute(f"create view {k} as select * from read_parquet({q(str(p))})")
    cols={k:set(con.execute(f"describe {k}").fetchdf()["column_name"].astype(str)) for k in paths}
    sp,se,ad,dc=first(cols['s'],['PATID']),first(cols['s'],['ENCOUNTERID']),first(cols['s'],['ADMIT_DATE']),first(cols['s'],['DISCHARGE_DATE'])
    pp,pe=first(cols['p'],['PATID']),first(cols['p'],['ENCOUNTERID']); od,ot=first(cols['p'],['RX_ORDER_DATE']),first(cols['p'],['RX_ORDER_TIME']); pn,pc,pr=first(cols['p'],['RAW_RX_MED_NAME']),first(cols['p'],['RXNORM_CUI']),first(cols['p'],['RX_ROUTE'])
    dp,dd=first(cols['d'],['PATID']),first(cols['d'],['DEATH_DATE'])
    inc,exc=parse_codes(root/'PCORnet/code/config/codes_antibiotics.py'); incsql=','.join(q(x) for x in sorted(inc)) or "''"; excsql=','.join(q(x) for x in sorted(exc)) or "''"
    pts=ts('p',od,ot); name=f"lower(coalesce(cast(p.{qi(pn)} as varchar),''))"; code=f"trim(coalesce(cast(p.{qi(pc)} as varchar),''))"; route=f"upper(trim(coalesce(cast(p.{qi(pr)} as varchar),'')))"
    broad=f"(({code} in ({incsql}) or regexp_matches({name},{q(BROAD_PATTERN)})) and {code} not in ({excsql}) and not regexp_matches({name},{q(NON_SYSTEMIC_PATTERN)}) and {route} not in ('ORAL','RESPIRATORY_TRACT','INHALATION'))"
    con.execute(f"create temp table base as select distinct cast({qi(sp)} as varchar) patid,cast({qi(se)} as varchar) encounterid,try_cast({qi(ad)} as date) admit_date,try_cast({qi(dc)} as date) discharge_date from s")
    con.execute(f"create temp table anchors as select cast(p.{qi(pp)} as varchar) patid,cast(p.{qi(pe)} as varchar) encounterid,min({pts}) anchor_ts from p join base b on cast(p.{qi(pp)} as varchar)=b.patid and cast(p.{qi(pe)} as varchar)=b.encounterid where {broad} and {pts}>=cast(b.admit_date as timestamp) and {pts}<cast(b.admit_date as timestamp)+interval 24 hour group by 1,2")
    con.execute(f"create temp table deaths as select cast({qi(dp)} as varchar) patid,min(try_cast({qi(dd)} as date)) death_date from d group by 1")
    con.execute("create temp table cohort as select a.*,b.discharge_date,d.death_date from anchors a join base b using(patid,encounterid) left join deaths d using(patid) where (b.discharge_date is null or b.discharge_date>cast(a.anchor_ts+interval 96 hour as date)) and (d.death_date is null or d.death_date>cast(a.anchor_ts+interval 96 hour as date))")
    n=int(con.execute("select count(*) from cohort").fetchone()[0])

    lp,le=first(cols['l'],['PATID']),first(cols['l'],['ENCOUNTERID']); lc=first(cols['l'],['LAB_LOINC']); lv=first(cols['l'],['RESULT_NUM'])
    sd,st=first(cols['l'],['SPECIMEN_DATE']),first(cols['l'],['SPECIMEN_TIME']); rd,rt=first(cols['l'],['RESULT_DATE']),first(cols['l'],['RESULT_TIME']); od2=first(cols['l'],['LAB_ORDER_DATE'])
    clock_exprs={"specimen":ts('l',sd,st),"result":ts('l',rd,rt),"order_date":ts('l',od2,None)}
    clock_rows=[]
    for concept,codes in LAB_CODES.items():
        csql=','.join(q(x) for x in codes)
        base_where=f"cast(l.{qi(lc)} as varchar) in ({csql}) and try_cast(l.{qi(lv)} as double) is not null"
        any_join=int(con.execute(f"select count(distinct c.encounterid) from cohort c join l on cast(l.{qi(lp)} as varchar)=c.patid and cast(l.{qi(le)} as varchar)=c.encounterid where {base_where}").fetchone()[0])
        for clock,expr in clock_exprs.items():
            valid=int(con.execute(f"select count(distinct c.encounterid) from cohort c join l on cast(l.{qi(lp)} as varchar)=c.patid and cast(l.{qi(le)} as varchar)=c.encounterid where {base_where} and {expr} is not null").fetchone()[0])
            for win,lo,hi in [('0_24h',0,24),('48_72h',48,72),('pre72',0,72)]:
                within=int(con.execute(f"select count(distinct c.encounterid) from cohort c join l on cast(l.{qi(lp)} as varchar)=c.patid and cast(l.{qi(le)} as varchar)=c.encounterid where {base_where} and {expr}>=c.anchor_ts+interval {lo} hour and {expr}<=c.anchor_ts+interval {hi} hour").fetchone()[0])
                clock_rows.append({'concept':concept,'clock':clock,'window':win,'cohort_any_matching_lab':safe(any_join),'cohort_valid_clock':safe(valid),'cohort_within_window':safe(within),'cohort_n':safe(n)})
    pd.DataFrame(clock_rows).to_csv(a.output_dir/'lab_clock_coverage.csv',index=False)

    op,oe=first(cols['o'],['PATID']),first(cols['o'],['ENCOUNTERID']); oc=first(cols['o'],['OBSCLIN_CODE']); on=first(cols['o'],['RAW_OBSCLIN_NAME']); ou=first(cols['o'],['OBSCLIN_RESULT_UNIT','RAW_OBSCLIN_UNIT']); ov=first(cols['o'],['OBSCLIN_RESULT_NUM'])
    map_rows=[]
    if on and oc and ov:
        nameexpr=f"lower(coalesce(cast(o.{qi(on)} as varchar),''))"
        mdf=con.execute(f"select cast(o.{qi(oc)} as varchar) code,coalesce(cast(o.{qi(on)} as varchar),'') raw_name,coalesce(cast(o.{qi(ou)} as varchar),'') unit,count(*) as row_count,count(distinct cast(o.{qi(oe)} as varchar)) as encounters from o where try_cast(o.{qi(ov)} as double) is not null and (strpos({nameexpr},'mean arterial')>0 or regexp_matches({nameexpr},'(^|[^a-z])map([^a-z]|$)')) group by 1,2,3 order by row_count desc limit 30").fetchdf()
        for _,r in mdf.iterrows():
            if int(r.row_count)>=MIN_CELL:
                map_rows.append({'code':r.code,'raw_name':r.raw_name,'unit':r.unit,'rows':safe(r.row_count),'encounters':safe(r.encounters)})
    pd.DataFrame(map_rows).to_csv(a.output_dir/'map_candidate_codes.csv',index=False)
    summary={'privacy_mode':'aggregate_only','strict_cohort_n':safe(n),'purpose':'resolve zero lab-window coverage and MAP mapping before covariate freeze','guardrail':'No propensity score or treatment effects fit.'}
    (a.output_dir/'summary.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__': main()

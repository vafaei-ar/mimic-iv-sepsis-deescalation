#!/usr/bin/env python3
"""Corrected aggregate-only PSU covariate mapping audit.

Uses exact clinical codes where they work, and conservative raw-name + unit fallbacks only
inside the strict PSU cohort when local LOINC population is absent. It also resolves a MAP
candidate from OBS_CLIN and excludes obvious non-systemic vasopressor false positives.
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
# Conservative local-name fallback rules, always paired with expected units and exclusions.
LAB_NAME_RULES = {
    "lactate": {
        "include": ["lactate", "lactic acid"],
        "exclude": ["dehydrogenase", "ldh", "csf", "fluid", "urine", "organic acid", "isoenzyme"],
        "units": ["mmol/l"],
    },
    "creatinine": {
        "include": ["creatinine level", "creatinine, poc", "| cret", "creat w/gfr"],
        "exclude": ["urine", "fluid", "clearance", "ratio", "protein/creatinine", "estimated gfr"],
        "units": ["mg/dl"],
    },
    "wbc": {
        "include": ["white blood cell count", "| wbc"],
        "exclude": ["urine", "fluid", "csf", "body fluid"],
        "units": ["10*3/ul", "10^3/ul", "k/ul"],
    },
    "platelet": {
        "include": ["platelet count", "| plts"],
        "exclude": ["morphology", "prepare platelets", "transfuse", "function", "antibody", "immature"],
        "units": ["10*3/ul", "10^3/ul", "k/ul"],
    },
    "bilirubin_total": {
        "include": ["total bilirubin", "bilirubin total", "bilirubin level"],
        "exclude": ["direct", "indirect", "conjugated", "neonatal", "fluid", "urine"],
        "units": ["mg/dl"],
    },
}
OBS_CODES = {
    "heart_rate": ["8867-4"],
    "resp_rate": ["9279-1"],
    "spo2": ["59408-5"],
    "temperature": ["8310-5"],
}
VASO_TERMS = ["norepinephrine", "levophed", "phenylephrine", "vasopressin", "epinephrine", "dopamine"]
VASO_EXCLUDE = ["racepinephrine", "nasal", "ophthalm", "lidocaine", "topical"]


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
                    if re.fullmatch(r"\d{4,9}", x): exc.add(x)
            elif isinstance(v,(list,tuple)):
                for x in v:
                    if isinstance(x,(list,tuple)) and x:
                        z=str(x[0]).strip()
                        if re.fullmatch(r"\d{4,9}", z): inc.add(z)
    return inc-exc, exc

def ts(alias, date_col, time_col):
    if not date_col: return "NULL::TIMESTAMP"
    if time_col:
        return f"try_cast(cast({alias}.{qi(date_col)} as varchar)||' '||coalesce(cast({alias}.{qi(time_col)} as varchar),'00:00:00') as timestamp)"
    return f"try_cast({alias}.{qi(date_col)} as timestamp)"

def str_contains(expr, terms):
    return "(" + " OR ".join(f"strpos({expr},{q(t)})>0" for t in terms) + ")"

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("data_root",type=Path); ap.add_argument("--output-dir",type=Path,required=True)
    a=ap.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True); root=a.data_root
    paths={"s":find(root,["sepsis_encounter"]),"p":find(root,["prescribing"]),"d":find(root,["death"]),"l":find(root,["lab_reduced","lab_result_cm"]),"o":find(root,["obs_clin"]),"m":find(root,["med_admin"])}
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

    # Core labs: audit exact LOINC and conservative cohort-local name/unit fallbacks.
    lp,le=first(cols['l'],['PATID']),first(cols['l'],['ENCOUNTERID']); lc=first(cols['l'],['LAB_LOINC']); ln=first(cols['l'],['RAW_LAB_NAME']); lu=first(cols['l'],['RESULT_UNIT','RAW_UNIT']); lv=first(cols['l'],['RESULT_NUM'])
    sd,st=first(cols['l'],['SPECIMEN_DATE']),first(cols['l'],['SPECIMEN_TIME']); rd,rt=first(cols['l'],['RESULT_DATE']),first(cols['l'],['RESULT_TIME'])
    specimen_ts=ts('l',sd,st); result_ts=ts('l',rd,rt)
    lname=f"lower(coalesce(cast(l.{qi(ln)} as varchar),''))"; lunit=f"lower(trim(coalesce(cast(l.{qi(lu)} as varchar),'')))"
    lab_map=[]; cov=[]
    for concept,codes in LAB_CODES.items():
        csql=','.join(q(x) for x in codes)
        rule=LAB_NAME_RULES[concept]
        incl=str_contains(lname,rule['include']); excl=str_contains(lname,rule['exclude']); units=','.join(q(x) for x in rule['units'])
        exact=f"cast(l.{qi(lc)} as varchar) in ({csql})"
        fallback=f"({incl} and not {excl} and {lunit} in ({units}))"
        # Map only aggregate cohort-local name/code/unit combinations.
        mdf=con.execute(f"select coalesce(cast(l.{qi(lc)} as varchar),'') code,coalesce(cast(l.{qi(ln)} as varchar),'') raw_name,coalesce(cast(l.{qi(lu)} as varchar),'') unit,case when {exact} then 'exact_loinc' else 'name_unit_fallback' end mapping_source,count(*) as row_count,count(distinct c.encounterid) encounters from cohort c join l on cast(l.{qi(lp)} as varchar)=c.patid and cast(l.{qi(le)} as varchar)=c.encounterid where try_cast(l.{qi(lv)} as double) is not null and ({exact} or {fallback}) group by 1,2,3,4 order by row_count desc limit 40").fetchdf()
        for _,r in mdf.iterrows():
            if int(r.row_count)>=MIN_CELL:
                lab_map.append({'concept':concept,'code':r.code,'raw_name':r.raw_name,'unit':r.unit,'mapping_source':r.mapping_source,'rows':safe(r.row_count),'encounters':safe(r.encounters)})
        # Coverage for exact and fallback, testing specimen and result clocks separately.
        for source,cond in [('exact_loinc',exact),('name_unit_fallback',fallback)]:
            for clock,clock_expr in [('specimen',specimen_ts),('result',result_ts)]:
                for win,lo,hi in [('0_24h',0,24),('48_72h',48,72),('pre72',0,72)]:
                    x=con.execute(f"select count(distinct c.encounterid) from cohort c join l on cast(l.{qi(lp)} as varchar)=c.patid and cast(l.{qi(le)} as varchar)=c.encounterid where try_cast(l.{qi(lv)} as double) is not null and {cond} and {clock_expr}>=c.anchor_ts+interval {lo} hour and {clock_expr}<=c.anchor_ts+interval {hi} hour").fetchone()[0]
                    cov.append({'domain':'lab','construct':concept+'_'+win,'mapping_source':source,'clock':clock,'count':safe(x),'cohort_n':safe(n),'proportion':x/n if n else None})
    pd.DataFrame(lab_map).to_csv(a.output_dir/'lab_code_unit_map.csv',index=False)

    # OBS_CLIN: exact codes for HR/RR/SpO2/temp; conservative name candidate for MAP.
    op,oe=first(cols['o'],['PATID']),first(cols['o'],['ENCOUNTERID']); oc=first(cols['o'],['OBSCLIN_CODE']); on=first(cols['o'],['RAW_OBSCLIN_NAME']); ou=first(cols['o'],['OBSCLIN_RESULT_UNIT','RAW_OBSCLIN_UNIT']); ov=first(cols['o'],['OBSCLIN_RESULT_NUM']); odt=first(cols['o'],['OBSCLIN_START_DATE']); ott=first(cols['o'],['OBSCLIN_START_TIME']); ots=ts('o',odt,ott)
    oname=f"lower(coalesce(cast(o.{qi(on)} as varchar),''))"; orows=[]
    for concept,codes in OBS_CODES.items():
        csql=','.join(q(x) for x in codes); cond=f"cast(o.{qi(oc)} as varchar) in ({csql})"
        udf=con.execute(f"select cast(o.{qi(oc)} as varchar) code,coalesce(cast(o.{qi(on)} as varchar),'') raw_name,coalesce(cast(o.{qi(ou)} as varchar),'') unit,'exact_code' mapping_source,count(*) as row_count,count(distinct c.encounterid) encounters from cohort c join o on cast(o.{qi(op)} as varchar)=c.patid and cast(o.{qi(oe)} as varchar)=c.encounterid where {cond} and try_cast(o.{qi(ov)} as double) is not null group by 1,2,3,4 order by row_count desc limit 20").fetchdf()
        for _,r in udf.iterrows():
            if int(r.row_count)>=MIN_CELL: orows.append({'concept':concept,'code':r.code,'raw_name':r.raw_name,'unit':r.unit,'mapping_source':'exact_code','rows':safe(r.row_count),'encounters':safe(r.encounters)})
        for win,lo,hi in [('0_24h',0,24),('48_72h',48,72),('pre72',0,72)]:
            x=con.execute(f"select count(distinct c.encounterid) from cohort c join o on cast(o.{qi(op)} as varchar)=c.patid and cast(o.{qi(oe)} as varchar)=c.encounterid where {cond} and try_cast(o.{qi(ov)} as double) is not null and {ots}>=c.anchor_ts+interval {lo} hour and {ots}<=c.anchor_ts+interval {hi} hour").fetchone()[0]
            cov.append({'domain':'obsclin','construct':concept+'_'+win,'mapping_source':'exact_code','clock':'obsclin_start','count':safe(x),'cohort_n':safe(n),'proportion':x/n if n else None})
    map_cond=f"(strpos({oname},'mean arterial')>0 or regexp_matches({oname},'(^|[^a-z])map([^a-z]|$)'))"
    mdf=con.execute(f"select coalesce(cast(o.{qi(oc)} as varchar),'') code,coalesce(cast(o.{qi(on)} as varchar),'') raw_name,coalesce(cast(o.{qi(ou)} as varchar),'') unit,count(*) as row_count,count(distinct c.encounterid) encounters from cohort c join o on cast(o.{qi(op)} as varchar)=c.patid and cast(o.{qi(oe)} as varchar)=c.encounterid where {map_cond} and try_cast(o.{qi(ov)} as double) is not null group by 1,2,3 order by row_count desc limit 30").fetchdf()
    for _,r in mdf.iterrows():
        if int(r.row_count)>=MIN_CELL: orows.append({'concept':'map','code':r.code,'raw_name':r.raw_name,'unit':r.unit,'mapping_source':'name_candidate','rows':safe(r.row_count),'encounters':safe(r.encounters)})
    for win,lo,hi in [('0_24h',0,24),('48_72h',48,72),('pre72',0,72)]:
        x=con.execute(f"select count(distinct c.encounterid) from cohort c join o on cast(o.{qi(op)} as varchar)=c.patid and cast(o.{qi(oe)} as varchar)=c.encounterid where {map_cond} and try_cast(o.{qi(ov)} as double) is not null and {ots}>=c.anchor_ts+interval {lo} hour and {ots}<=c.anchor_ts+interval {hi} hour").fetchone()[0]
        cov.append({'domain':'obsclin','construct':'map_'+win,'mapping_source':'name_candidate','clock':'obsclin_start','count':safe(x),'cohort_n':safe(n),'proportion':x/n if n else None})
    pd.DataFrame(orows).to_csv(a.output_dir/'obsclin_code_unit_map.csv',index=False)
    pd.DataFrame(cov).to_csv(a.output_dir/'strict_cohort_coverage.csv',index=False)

    # Vasopressor aggregate map with obvious non-systemic formulations removed.
    me=first(cols['m'],['ENCOUNTERID']); mn=first(cols['m'],['RAW_MEDADMIN_MED_NAME']); mr=first(cols['m'],['MEDADMIN_ROUTE'])
    mname=f"lower(coalesce(cast(m.{qi(mn)} as varchar),''))"; mroute=f"upper(trim(coalesce(cast(m.{qi(mr)} as varchar),'')))"
    incl=' or '.join(f"strpos({mname},{q(x)})>0" for x in VASO_TERMS); excl=' or '.join(f"strpos({mname},{q(x)})>0" for x in VASO_EXCLUDE)
    vdf=con.execute(f"select case when strpos({mname},'norepinephrine')>0 or strpos({mname},'levophed')>0 then 'norepinephrine' when strpos({mname},'phenylephrine')>0 then 'phenylephrine' when strpos({mname},'vasopressin')>0 then 'vasopressin' when strpos({mname},'epinephrine')>0 then 'epinephrine' when strpos({mname},'dopamine')>0 then 'dopamine' else 'other' end agent,{mroute} route,count(*) as row_count,count(distinct cast(m.{qi(me)} as varchar)) encounters from m where ({incl}) and not ({excl}) group by 1,2 order by row_count desc").fetchdf(); vdf=vdf.rename(columns={'row_count':'rows'}); vdf['rows']=vdf['rows'].map(safe); vdf['encounters']=vdf['encounters'].map(safe); vdf.to_csv(a.output_dir/'vasopressor_validated_map.csv',index=False)

    gcs_rows=0
    if on:
        gcs_rows=int(con.execute(f"select count(*) from cohort c join o on cast(o.{qi(op)} as varchar)=c.patid and cast(o.{qi(oe)} as varchar)=c.encounterid where strpos({oname},'glasgow')>0 or regexp_matches({oname},'(^|[^a-z])gcs([^a-z]|$)')").fetchone()[0])
    summary={'privacy_mode':'aggregate_only','strict_cohort_n':safe(n),'lab_mapping':'exact LOINC plus conservative cohort-local name/unit fallback audit','obsclin_mapping':'exact codes for HR/RR/SpO2/temperature plus MAP name-candidate audit','gcs_candidate_rows_strict_cohort':safe(gcs_rows),'fio2_primary':'exclude','ventilation_primary':'unavailable_current_extract','guardrail':'No propensity score or treatment effects fit; freeze only after reviewing these aggregate maps.'}
    (a.output_dir/'summary.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__': main()

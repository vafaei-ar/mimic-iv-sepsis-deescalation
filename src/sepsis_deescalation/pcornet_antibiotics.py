from __future__ import annotations

import re
from typing import Iterable

import pandas as pd

# RxNorm mappings inherited from prior PSU work. These must be audited against the actual
# PSU datamart before the external validation is frozen.
RXCUI_TO_DRUG_GROUP = {
    "11124": "vancomycin", "202368": "vancomycin", "239209": "vancomycin",
    "74170": "piperacillin_tazobactam", "8339": "piperacillin_tazobactam", "203134": "piperacillin_tazobactam",
    "20481": "cefepime",
    "2191": "ceftazidime", "203421": "ceftazidime",
    "29561": "meropenem", "83175": "meropenem",
    "5690": "imipenem_cilastatin", "33533": "imipenem_cilastatin",
    "1272": "aztreonam", "202561": "aztreonam",
    "190376": "linezolid", "261710": "linezolid",
    "22299": "daptomycin", "404965": "daptomycin",
    "1597616": "ceftolozane_tazobactam",
    "1603841": "ceftazidime_avibactam",
}

NAME_PATTERNS = {
    "vancomycin": ["vancomycin", "vancocin"],
    "piperacillin_tazobactam": ["piperacillin", "zosyn", "tazobactam"],
    "cefepime": ["cefepime", "maxipime"],
    "ceftazidime": ["ceftazidime", "fortaz", "tazicef"],
    "meropenem": ["meropenem", "merrem"],
    "imipenem_cilastatin": ["imipenem", "primaxin"],
    "aztreonam": ["aztreonam", "azactam"],
    "linezolid": ["linezolid", "zyvox"],
    "daptomycin": ["daptomycin", "cubicin"],
    "ceftolozane_tazobactam": ["ceftolozane", "zerbaxa"],
    "ceftazidime_avibactam": ["avibactam", "avycaz"],
}

NON_IV_NAME_PATTERNS = [
    "oral", "capsule", "tablet", "suspension", "powder for oral", "enema", "ophthalmic",
    "topical", "inhalation", "nebulizer", "cream", "ointment", "drops",
]
NON_IV_RXCUI = {"313570", "313571", "2000134"}
NON_IV_ROUTES = {
    "ORAL", "TOPICAL", "OPHTHALMIC", "RECTAL", "NASAL", "OTIC", "SUBLINGUAL", "BUCCAL",
    "TRANSDERMAL", "OROMUCOSAL", "VAGINAL", "INTRAVESICAL", "INTRAUTERINE", "RESPIRATORY_TRACT",
    "INTRATHECAL", "INTRADERMAL", "INTRAOCULAR", "INTRAVITREAL", "INTRA_ARTICULAR",
}

NON_BROAD_NAMES = [
    "ceftriaxone", "cefazolin", "ampicillin", "amoxicillin", "doxycycline", "azithromycin",
    "metronidazole", "clindamycin", "cephalexin", "ciprofloxacin", "levofloxacin", "gentamicin", "tobramycin",
]


def classify_broad_drug(code: object, raw_name: object, route: object = None) -> str | None:
    code_str = str(code).strip() if code is not None else ""
    raw = str(raw_name).lower().strip() if raw_name is not None else ""
    route_str = str(route).upper().strip() if route is not None else ""
    if code_str in NON_IV_RXCUI or any(p in raw for p in NON_IV_NAME_PATTERNS) or route_str in NON_IV_ROUTES:
        return None
    if code_str in RXCUI_TO_DRUG_GROUP:
        return RXCUI_TO_DRUG_GROUP[code_str]
    for group, patterns in NAME_PATTERNS.items():
        if any(pattern in raw for pattern in patterns):
            return group
    return None


def broad_mask(df: pd.DataFrame, code_col: str, name_col: str, route_col: str | None = None) -> pd.Series:
    route = df[route_col] if route_col and route_col in df else pd.Series(None, index=df.index)
    values = [classify_broad_drug(c, n, r) for c, n, r in zip(df[code_col], df[name_col], route)]
    return pd.Series([v is not None for v in values], index=df.index)


def add_broad_group(df: pd.DataFrame, code_col: str, name_col: str, route_col: str | None = None) -> pd.DataFrame:
    d = df.copy()
    route = d[route_col] if route_col and route_col in d else pd.Series(None, index=d.index)
    d["broad_drug_group"] = [classify_broad_drug(c, n, r) for c, n, r in zip(d[code_col], d[name_col], route)]
    d["is_broad_spectrum"] = d["broad_drug_group"].notna().astype(int)
    return d


def antibiotic_mapping_audit(df: pd.DataFrame, code_col: str, name_col: str, route_col: str | None = None) -> pd.DataFrame:
    cols = [code_col, name_col] + ([route_col] if route_col and route_col in df else [])
    d = add_broad_group(df[cols].copy(), code_col, name_col, route_col)
    return d.groupby(cols + ["broad_drug_group"], dropna=False).size().reset_index(name="n").sort_values("n", ascending=False)

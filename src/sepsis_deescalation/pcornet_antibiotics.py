from __future__ import annotations

import pandas as pd

# RxNorm mappings inherited from prior PSU work. Audit against the current PSU datamart before freezing.
RXCUI_TO_DRUG_GROUP = {
    "11124": "vancomycin", "202368": "vancomycin", "239209": "vancomycin",
    "74170": "piperacillin_tazobactam", "8339": "piperacillin_tazobactam", "203134": "piperacillin_tazobactam",
    "20481": "cefepime", "2191": "ceftazidime", "203421": "ceftazidime",
    "29561": "meropenem", "83175": "meropenem", "5690": "imipenem_cilastatin", "33533": "imipenem_cilastatin",
    "1272": "aztreonam", "202561": "aztreonam", "190376": "linezolid", "261710": "linezolid",
    "22299": "daptomycin", "404965": "daptomycin", "1597616": "ceftolozane_tazobactam", "1603841": "ceftazidime_avibactam",
}
NAME_PATTERNS = {
    "vancomycin": ["vancomycin", "vancocin"], "piperacillin_tazobactam": ["piperacillin", "zosyn", "tazobactam"],
    "cefepime": ["cefepime", "maxipime"], "ceftazidime": ["ceftazidime", "fortaz", "tazicef"],
    "meropenem": ["meropenem", "merrem"], "imipenem_cilastatin": ["imipenem", "primaxin"],
    "aztreonam": ["aztreonam", "azactam"], "linezolid": ["linezolid", "zyvox"],
    "daptomycin": ["daptomycin", "cubicin"], "ceftolozane_tazobactam": ["ceftolozane", "zerbaxa"],
    "ceftazidime_avibactam": ["avibactam", "avycaz"],
}
NON_IV_NAME_PATTERNS = ["oral", "capsule", "tablet", "suspension", "powder for oral", "enema", "ophthalmic", "topical", "inhalation", "nebulizer", "cream", "ointment", "drops"]
NON_IV_RXCUI = {"313570", "313571", "2000134"}
NON_IV_ROUTES = {"ORAL", "TOPICAL", "OPHTHALMIC", "RECTAL", "NASAL", "OTIC", "SUBLINGUAL", "BUCCAL", "TRANSDERMAL", "OROMUCOSAL", "VAGINAL", "INTRAVESICAL", "INTRAUTERINE", "RESPIRATORY_TRACT", "INTRATHECAL", "INTRADERMAL", "INTRAOCULAR", "INTRAVITREAL", "INTRA_ARTICULAR"}
NON_BROAD_NAMES = ["ceftriaxone", "cefazolin", "ampicillin", "amoxicillin", "doxycycline", "azithromycin", "metronidazole", "clindamycin", "cephalexin", "ciprofloxacin", "levofloxacin", "gentamicin", "tobramycin"]


def _clean(value: object) -> str:
    return "" if value is None or pd.isna(value) else str(value).strip()


def is_non_iv(code: object, raw_name: object, route: object = None) -> bool:
    code_str = _clean(code)
    raw = _clean(raw_name).lower()
    route_str = _clean(route).upper()
    return code_str in NON_IV_RXCUI or any(p in raw for p in NON_IV_NAME_PATTERNS) or route_str in NON_IV_ROUTES


def classify_broad_drug(code: object, raw_name: object, route: object = None) -> str | None:
    if is_non_iv(code, raw_name, route):
        return None
    code_str = _clean(code)
    raw = _clean(raw_name).lower()
    if code_str in RXCUI_TO_DRUG_GROUP:
        return RXCUI_TO_DRUG_GROUP[code_str]
    for group, patterns in NAME_PATTERNS.items():
        if any(pattern in raw for pattern in patterns):
            return group
    return None


def is_observed_systemic_antibiotic(code: object, raw_name: object, route: object = None) -> bool:
    if is_non_iv(code, raw_name, route):
        return False
    if classify_broad_drug(code, raw_name, route) is not None:
        return True
    raw = _clean(raw_name).lower()
    return any(name in raw for name in NON_BROAD_NAMES)


def add_antibiotic_classification(df: pd.DataFrame, code_col: str, name_col: str, route_col: str | None = None) -> pd.DataFrame:
    d = df.copy()
    route = d[route_col] if route_col and route_col in d else pd.Series(None, index=d.index)
    d["broad_drug_group"] = [classify_broad_drug(c, n, r) for c, n, r in zip(d[code_col], d[name_col], route)]
    d["is_broad_spectrum"] = d["broad_drug_group"].notna().astype(int)
    d["is_systemic_antibiotic"] = [int(is_observed_systemic_antibiotic(c, n, r)) for c, n, r in zip(d[code_col], d[name_col], route)]
    return d


def broad_mask(df: pd.DataFrame, code_col: str, name_col: str, route_col: str | None = None) -> pd.Series:
    return add_antibiotic_classification(df, code_col, name_col, route_col)["is_broad_spectrum"].astype(bool)


def antibiotic_mapping_audit(df: pd.DataFrame, code_col: str, name_col: str, route_col: str | None = None) -> pd.DataFrame:
    cols = [code_col, name_col] + ([route_col] if route_col and route_col in df else [])
    d = add_antibiotic_classification(df[cols].copy(), code_col, name_col, route_col)
    return d.groupby(cols + ["broad_drug_group", "is_broad_spectrum", "is_systemic_antibiotic"], dropna=False).size().reset_index(name="n").sort_values("n", ascending=False)

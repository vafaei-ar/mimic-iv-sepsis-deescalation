#!/usr/bin/env python3
"""Compatibility wrapper for the final PSU covariate freeze using date-span MED_ADMIN timing.

MED_ADMIN vasopressor rows in the audited PSU extract have date-level start/stop information
without usable administration TIME values. For interval-overlap construction, start dates are
interpreted at 00:00:00 and stop dates at 23:59:59 on the recorded calendar date. All other
clock handling and covariate definitions are delegated unchanged to the frozen base audit.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import audit_psu_final_covariate_freeze_base as base

_original_ts = base.ts


def _ts_with_medadmin_date_span(alias, date_col, time_col, numeric_seconds=False):
    if alias == "m" and date_col:
        d = f"try_cast({alias}.{base.qi(date_col)} as date)"
        if "STOP" in str(date_col).upper():
            return f"case when {d} is null then null else cast({d} as timestamp)+interval 1 day-interval 1 second end"
        return f"case when {d} is null then null else cast({d} as timestamp) end"
    return _original_ts(alias, date_col, time_col, numeric_seconds)


base.ts = _ts_with_medadmin_date_span
base.main()

# Record the externally validated MED_ADMIN clock choice in aggregate outputs.
try:
    out = Path(sys.argv[sys.argv.index("--output-dir") + 1])
    summary_path = out / "summary.json"
    if summary_path.exists():
        s = json.loads(summary_path.read_text())
        s["medadmin_clock"] = "date-span: start date 00:00:00; stop date 23:59:59 because audited MEDADMIN TIME fields are unusable"
        s["medadmin_clock_validation_job"] = "MIMICIV-SEPSIS-DEESCALATION-PSU-MEDADMIN-CLOCK-SEMANTICS-0003"
        summary_path.write_text(json.dumps(s, indent=2))
except Exception:
    pass

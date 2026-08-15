#!/usr/bin/env python
from __future__ import annotations

import argparse
import logging

from sepsis_deescalation.pcornet_pipeline import run_pcornet


def main() -> None:
    parser = argparse.ArgumentParser(description="Run harmonized PSU/PCORnet cohort and exposure pipeline.")
    parser.add_argument("--config", default="config/pcornet_psu.local.yaml")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    run_dir = run_pcornet(args.config)
    print(f"Completed PSU/PCORnet cohort/exposure stage: {run_dir}")
    print("This stage intentionally does not emit the final adjusted external effect until PSU covariate mappings are audited and frozen.")


if __name__ == "__main__":
    main()

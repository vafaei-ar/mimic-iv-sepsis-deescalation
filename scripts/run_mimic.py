#!/usr/bin/env python
from __future__ import annotations

import argparse
import logging

from sepsis_deescalation.mimic_pipeline import run_mimic


def main() -> None:
    parser = argparse.ArgumentParser(description="Run MIMIC-IV day-3 antibiotic de-escalation analysis.")
    parser.add_argument("--config", default="config/mimic.yaml")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    run_dir = run_mimic(args.config)
    print(f"Completed MIMIC analysis: {run_dir}")
    print(f"ZIP: {run_dir}.zip")


if __name__ == "__main__":
    main()

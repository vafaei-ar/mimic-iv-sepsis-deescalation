#!/usr/bin/env python3
"""Run repository-level public validation without accessing restricted data."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

LOG_PATH = Path("outputs/validation/repository_validation.log")


def run_and_log(log, *args: str) -> int:
    proc = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log.write("$ " + " ".join(args) + "\n")
    log.write(proc.stdout or "")
    if proc.stdout and not proc.stdout.endswith("\n"):
        log.write("\n")
    log.write(f"[exit {proc.returncode}]\n\n")
    log.flush()
    return int(proc.returncode)


def main() -> None:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LOG_PATH.open("w", encoding="utf-8") as log:
        rc_pytest = run_and_log(log, sys.executable, "-m", "pytest", "-q")
        rc_ruff = run_and_log(log, sys.executable, "-m", "ruff", "check", ".")
    raise SystemExit(rc_pytest or rc_ruff)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run repository-level public validation without accessing restricted data."""
from __future__ import annotations

import subprocess
import sys


def run(*args: str) -> None:
    subprocess.run(args, check=True)


def main() -> None:
    run(sys.executable, "-m", "pytest", "-q")
    run(sys.executable, "-m", "ruff", "check", ".")


if __name__ == "__main__":
    main()

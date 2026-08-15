from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml


def load_config(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    with path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    cfg["_config_path"] = str(path.resolve())
    return cfg


def canonical_config_json(cfg: dict[str, Any]) -> str:
    clean = {k: v for k, v in cfg.items() if not str(k).startswith("_")}
    return json.dumps(clean, sort_keys=True, separators=(",", ":"), default=str)


def config_sha256(cfg: dict[str, Any]) -> str:
    return hashlib.sha256(canonical_config_json(cfg).encode("utf-8")).hexdigest()


def resolve_mimic_source(cfg: dict[str, Any]) -> Path:
    env_name = cfg.get("source_env", "MIMIC_SOURCE")
    candidates: list[Path] = []
    if os.environ.get(env_name):
        candidates.append(Path(os.environ[env_name]).expanduser())
    candidates.extend(Path(p).expanduser() for p in cfg.get("source_candidates", []))
    for path in candidates:
        if path.exists():
            return path.resolve()
    attempted = "\n".join(f"  - {p}" for p in candidates) or "  (none)"
    raise FileNotFoundError(
        f"Could not locate MIMIC-IV. Set {env_name} or edit source_candidates. Tried:\n{attempted}"
    )


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


@dataclass(frozen=True)
class RunPaths:
    run_dir: Path
    tables: Path
    diagnostics: Path
    figures: Path
    logs: Path
    audits: Path


def make_run_paths(output_root: str | Path, site: str, version: str) -> RunPaths:
    root = Path(output_root)
    run_dir = root / f"{site}_{version}_{utc_timestamp()}"
    tables = run_dir / "tables"
    diagnostics = run_dir / "diagnostics"
    figures = run_dir / "figure_data"
    logs = run_dir / "logs"
    audits = run_dir / "audits"
    for p in (run_dir, tables, diagnostics, figures, logs, audits):
        p.mkdir(parents=True, exist_ok=True)
    return RunPaths(run_dir, tables, diagnostics, figures, logs, audits)

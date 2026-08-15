from __future__ import annotations

import json
import platform
import shutil
import subprocess
import sys
from importlib import metadata
from pathlib import Path
from typing import Any

from .config import config_sha256


def git_commit() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return None


def package_versions(names: list[str]) -> dict[str, str | None]:
    out: dict[str, str | None] = {}
    for name in names:
        try:
            out[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            out[name] = None
    return out


def write_run_manifest(
    run_dir: str | Path,
    cfg: dict[str, Any],
    source: str | Path | None,
    extra: dict[str, Any] | None = None,
) -> Path:
    run_dir = Path(run_dir)
    manifest = {
        "git_commit": git_commit(),
        "config_sha256": config_sha256(cfg),
        "config_path": cfg.get("_config_path"),
        "source": str(source) if source is not None else None,
        "python": sys.version,
        "platform": platform.platform(),
        "packages": package_versions(
            ["numpy", "pandas", "scipy", "statsmodels", "matplotlib", "PyYAML", "duckdb", "pyarrow"]
        ),
    }
    if extra:
        manifest.update(extra)
    path = run_dir / "run_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return path


def copy_config(cfg: dict[str, Any], run_dir: str | Path) -> Path | None:
    source = cfg.get("_config_path")
    if not source:
        return None
    source_path = Path(source)
    if not source_path.exists():
        return None
    target = Path(run_dir) / "analysis_config.yaml"
    shutil.copy2(source_path, target)
    return target


def zip_run(run_dir: str | Path) -> Path:
    run_dir = Path(run_dir)
    archive_base = run_dir.parent / run_dir.name
    path = shutil.make_archive(str(archive_base), "zip", root_dir=run_dir)
    return Path(path)

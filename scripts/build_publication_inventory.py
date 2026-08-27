#!/usr/bin/env python3
"""Create a safe inventory of project-local aggregate publication outputs.

This helper does not read raw MIMIC-IV or PSU patient-level data. It scans only the
project's ``outputs`` directory, records file metadata, and captures small text previews
for aggregate JSON/CSV/TSV/MD/TXT files so the frozen result packages can be identified
before manuscript integration.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"
DEST = OUTPUTS / "publication_integration" / "inventory"
TEXT_EXTENSIONS = {".json", ".csv", ".tsv", ".md", ".txt"}
MAX_PREVIEW_BYTES = 12000


def safe_preview(path: Path) -> str:
    if path.suffix.lower() not in TEXT_EXTENSIONS:
        return ""
    if path.stat().st_size > MAX_PREVIEW_BYTES:
        return ""
    try:
        return path.read_text(errors="replace")
    except OSError:
        return ""


def main() -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in sorted(OUTPUTS.rglob("*")):
        if not path.is_file():
            continue
        try:
            rel = path.relative_to(ROOT)
        except ValueError:
            continue
        if str(rel).startswith("outputs/publication_integration/inventory/"):
            continue
        stat = path.stat()
        rows.append(
            {
                "path": str(rel),
                "suffix": path.suffix.lower(),
                "size_bytes": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "preview": safe_preview(path),
            }
        )

    json_path = DEST / "publication_output_inventory.json"
    csv_path = DEST / "publication_output_inventory.csv"
    json_path.write_text(json.dumps({"n_files": len(rows), "files": rows}, indent=2))
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["path", "suffix", "size_bytes", "mtime_ns"])
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in writer.fieldnames})

    print(json.dumps({"n_files": len(rows), "json": str(json_path), "csv": str(csv_path)}, indent=2))


if __name__ == "__main__":
    main()

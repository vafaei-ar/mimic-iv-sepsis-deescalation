#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

EXPECTED_ROOT = Path('/home/asadr/Depts/PHS/PATH_CDM/Hwang_Bonavia').resolve()

TARGETS = [
    'PCORnet/parquet/sepsis_encounter.parquet',
    'PCORnet/parquet/sepsis_demographic.parquet',
    'PCORnet/parquet/sepsis_diagnosis.parquet',
    'PCORnet/parquet/sepsis_vital.parquet',
    'PCORnet/parquet/condition.parquet',
    'PCORnet/parquet/death.parquet',
    'PCORnet/parquet/prescribing.parquet',
    'PCORnet/parquet/med_admin.parquet',
    'PCORnet/parquet/obs_clin.parquet',
    'PCORnet/parquet/lab/lab_reduced.parquet',
    'PCORnet/Full/lab_result_cm.sas7bdat',
    'PCORnet/Full/prescribing.sas7bdat',
    'PCORnet/Full/med_admin.sas7bdat',
]


def safe_root(root: Path) -> Path:
    root = root.expanduser().resolve()
    if root != EXPECTED_ROOT:
        raise ValueError(f'Refusing unexpected data root: {root}')
    return root


def mount_info_for(path: Path) -> dict | None:
    best = None
    best_len = -1
    try:
        lines = Path('/proc/self/mountinfo').read_text(encoding='utf-8', errors='replace').splitlines()
    except Exception as exc:
        return {'error': f'{type(exc).__name__}: {exc}'}
    for line in lines:
        parts = line.split(' - ', 1)
        if len(parts) != 2:
            continue
        left, right = parts
        fields = left.split()
        rfields = right.split()
        if len(fields) < 6 or len(rfields) < 3:
            continue
        mount_point = fields[4].replace('\\040', ' ')
        try:
            path.relative_to(Path(mount_point))
        except ValueError:
            continue
        if len(mount_point) > best_len:
            best_len = len(mount_point)
            best = {
                'mount_point': mount_point,
                'filesystem_type': rfields[0],
                'mount_source': rfields[1],
                'mount_options': fields[5],
            }
    return best


def probe(path: Path, root: Path) -> dict:
    row = {'relative_path': str(path.relative_to(root))}
    for name, fn in [
        ('lexists', lambda: os.path.lexists(path)),
        ('lstat', lambda: os.lstat(path)),
        ('stat', lambda: os.stat(path)),
    ]:
        try:
            value = fn()
            row[name] = True if name == 'lexists' else 'ok'
            if name == 'stat':
                row['size_bytes'] = int(value.st_size)
        except Exception as exc:
            row[name] = False if name == 'lexists' else 'error'
            row[f'{name}_error'] = f'{type(exc).__name__}: {exc}'
    try:
        fd = os.open(path, os.O_RDONLY | getattr(os, 'O_NONBLOCK', 0))
        os.close(fd)
        row['open_readonly'] = 'ok'
    except Exception as exc:
        row['open_readonly'] = 'error'
        row['open_error'] = f'{type(exc).__name__}: {exc}'
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description='Metadata-only PSU filesystem mount diagnostic; no file contents are read or exported.')
    parser.add_argument('data_root', type=Path)
    parser.add_argument('--output', type=Path, default=Path('outputs/psu_mount_diagnostic/latest/summary.json'))
    args = parser.parse_args()

    root = safe_root(args.data_root)
    out = args.output.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    parent_checks = []
    for rel in ['', 'PCORnet', 'PCORnet/parquet', 'PCORnet/Full']:
        p = root / rel
        try:
            os.stat(p)
            status = 'ok'
            error = None
        except Exception as exc:
            status = 'error'
            error = f'{type(exc).__name__}: {exc}'
        parent_checks.append({'relative_path': rel or '.', 'stat': status, 'error': error})

    probes = [probe(root / rel, root) for rel in TARGETS]
    error_counts: dict[str, int] = {}
    for row in probes:
        for key in ('lstat_error', 'stat_error', 'open_error'):
            msg = row.get(key)
            if not msg:
                continue
            bucket = 'stale_file_handle' if 'Stale file handle' in msg or '[Errno 116]' in msg else msg.split(':', 1)[0]
            error_counts[bucket] = error_counts.get(bucket, 0) + 1

    summary = {
        'privacy_mode': 'filesystem_metadata_only_no_patient_rows_no_file_contents',
        'data_root': str(root),
        'mount': mount_info_for(root),
        'parent_checks': parent_checks,
        'targets_tested': len(probes),
        'targets_openable': sum(1 for x in probes if x.get('open_readonly') == 'ok'),
        'error_counts': error_counts,
        'probes': probes,
        'interpretation_hint': 'If parent directories stat successfully but target files repeatedly return Errno 116, refresh/remount the network filesystem before cohort extraction.',
    }
    out.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: summary[k] for k in ['privacy_mode','targets_tested','targets_openable','error_counts']}, indent=2))


if __name__ == '__main__':
    main()

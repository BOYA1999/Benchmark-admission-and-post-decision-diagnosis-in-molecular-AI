from pathlib import Path
from collections import Counter
import ast
import csv
import gzip
import hashlib
import os
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / 'MANIFEST_SHA256.txt'

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

failures = []
files = sorted(p for p in ROOT.rglob('*') if p.is_file())
listed = {}
for line in MANIFEST.read_text(encoding='utf-8').splitlines():
    try:
        expected, relative = line.split('  ', 1)
    except ValueError:
        failures.append('malformed manifest entry')
        continue
    target = (ROOT / relative).resolve()
    if ROOT not in target.parents or relative in listed:
        failures.append('unsafe or duplicate manifest entry')
        continue
    listed[relative] = expected
    if not target.is_file() or digest(target) != expected:
        failures.append('manifest content mismatch')
actual = {p.relative_to(ROOT).as_posix() for p in files if p != MANIFEST}
if set(listed) != actual:
    failures.append('manifest file set mismatch')

forbidden_parts = {'.git', '__pycache__', '.venv', 'venv', 'cache', 'caches', 'logs', 'checkpoints', 'score_cache'}
workflow_name = re.compile(r'(^|[_.-])(manuscript|dingzhen|review|reviewer|rebuttal|backups?)([_.-]|$)', re.I)
forbidden_suffixes = {'.doc', '.docx', '.tex', '.bib', '.log', '.pt', '.pth', '.ckpt', '.npz', '.npy', '.zip', '.7z', '.rar'}
forbidden_names = {'split_manifest.csv', 'predictions_long.csv', 'overlap_examples.csv'}
structure_columns = {'smiles', 'source_smiles', 'canonical_smiles', 'reactant', 'reagent', 'product'}
text_suffixes = {'', '.md', '.txt', '.csv', '.json', '.py', '.yml', '.yaml', '.toml', '.svg', '.mplstyle'}
patterns = [
    re.compile(r'(?<![A-Za-z0-9])[A-Za-z]:[\\/]'),
    re.compile(r'(?:^|["\'])\\\\{2,}[A-Za-z0-9_.-]+\\'),
    re.compile(r'/(?:home|Users)/[^/\s]+/', re.I),
    re.compile(r'(?<![A-Za-z0-9_.+-])[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}(?![A-Za-z0-9_.-])'),
    re.compile('BEGIN ' + 'PRIVATE KEY'),
    re.compile('gh' + r'[pousr]_[A-Za-z0-9]{20,}'),
    re.compile('AKIA' + r'[A-Z0-9]{16}'),
    re.compile('sk' + r'-[A-Za-z0-9_-]{20,}'),
]
if any(count > 100 for count in Counter(path.parent for path in files).values()):
    failures.append('directory contains more than 100 files')
for path in files:
    relative = path.relative_to(ROOT)
    parts = {part.lower() for part in relative.parts}
    if parts & forbidden_parts or path.suffix.lower() in forbidden_suffixes or any(workflow_name.search(part) for part in relative.parts) or path.name == '.env' or path.name.lower() in forbidden_names:
        failures.append('disallowed repository artifact')
    if path.suffix.lower() in text_suffixes or path.name == '.gitignore' or path.suffix.lower() == '.gz':
        if path.suffix.lower() == '.gz':
            try:
                text = gzip.decompress(path.read_bytes()).decode('utf-8-sig', errors='replace')
            except (OSError, EOFError):
                failures.append('unreadable compressed content')
                continue
        else:
            text = path.read_text(encoding='utf-8-sig', errors='replace')
        if any(pattern.search(text) for pattern in patterns):
            failures.append('privacy or credential pattern')
        if path.name.lower().endswith(('.csv', '.csv.gz')):
            header = {column.strip().lower() for column in next(csv.reader(text.splitlines()), [])}
            if header & structure_columns:
                failures.append('structure-bearing data table')
        if path.suffix.lower() == '.py':
            try:
                ast.parse(text, filename=relative.as_posix())
            except SyntaxError:
                failures.append('Python syntax error')
    else:
        payload = path.read_bytes()
        if re.search(rb'[A-Za-z]:[\\/](?:Users|Documents and Settings)[\\/]|/(?:Users|home)/[^/\s]+/', payload):
            failures.append('binary local-path marker')

if failures:
    raise SystemExit('FAIL: surface verification; categories=' + ','.join(sorted(set(failures))))

environment = os.environ.copy()
environment['PYTHONDONTWRITEBYTECODE'] = '1'
commands = [
    [sys.executable, '-B', str(ROOT / 'audit_cases/verify_release.py')],
    [sys.executable, '-B', str(ROOT / 'molxai_crc/scripts/verify_release.py')],
    [sys.executable, '-B', '-m', 'unittest', 'discover', '-s', str(ROOT / 'molxai_crc/src'), '-p', 'test_*.py'],
]
for command in commands:
    completed = subprocess.run(command, cwd=ROOT, env=environment)
    if completed.returncode:
        raise SystemExit('FAIL: nested scientific verification')

if any(p.is_file() and p.relative_to(ROOT).as_posix() not in actual | {'MANIFEST_SHA256.txt'} for p in ROOT.rglob('*')):
    raise SystemExit('FAIL: verification generated an unmanifested file')
print(f'PASS: {len(files)} files; manifest, privacy, archive policy, directory file limits, nested evidence checks, and unit tests')

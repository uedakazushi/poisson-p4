#!/usr/bin/env python3
"""Validate the current classification plan; never certify mathematics or start jobs.

Python >=3.11, standard library only. No network, no writes, no scheduler.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import sys
from typing import Any

PLAN_PATH = 'research/classification/PLAN.json'
STATES = {'ready', 'blocked', 'conditional', 'suspended', 'done'}
KINDS = {'primary', 'optional_primary', 'secondary'}

class PlanError(ValueError):
    pass

def relative_path(value: Any) -> str:
    if not isinstance(value, str) or not value or '\\' in value:
        raise PlanError(f'Invalid relative path: {value!r}')
    p = PurePosixPath(value)
    if p.is_absolute() or '..' in p.parts or str(p) != value or value == '.':
        raise PlanError(f'Unsafe/noncanonical path: {value!r}')
    return value

def local_file(root: Path, value: str) -> Path:
    p = root / relative_path(value)
    try:
        p.resolve().relative_to(root.resolve())
    except ValueError as exc:
        raise PlanError(f'Path escapes repository: {value}') from exc
    if not p.is_file():
        raise PlanError(f'Missing file: {value}')
    return p

def read_json(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, ValueError) as exc:
        raise PlanError(f'Cannot read JSON: {path}: {exc}') from exc
    if not isinstance(data, dict):
        raise PlanError(f'Expected JSON object: {path}')
    return data

def digest(path: Path, *, git_blob: bool = False) -> str:
    h = hashlib.sha1() if git_blob else hashlib.sha256()
    if git_blob:
        h.update(f'blob {path.stat().st_size}\0'.encode())
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def nonempty(value: Any, field: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise PlanError(f'Required nonempty string: {field}')

def validate(plan: dict, sources: dict, root: Path, *, check_inputs: bool = True) -> dict:
    if plan.get('schema_version') != 1:
        raise PlanError('Unsupported schema_version')
    for key in ('version', 'date', 'baseline_commit', 'target'):
        nonempty(plan.get(key), key)
    if sources.get('baseline_commit') != plan['baseline_commit']:
        raise PlanError('Source/plan baseline mismatch')
    registry = sources.get('repository_sources', [])
    ids = [s.get('id') for s in registry]
    if not ids or len(set(ids)) != len(ids) or any(not isinstance(i, str) for i in ids):
        raise PlanError('Source ids must be unique nonempty strings')
    source_ids = set(ids)
    for s in registry:
        relative_path(s['path'])
        for field in ('purpose', 'review_scope'):
            nonempty(s.get(field), field)
        if check_inputs:
            p = local_file(root, s['path'])
            expected = s.get('baseline_git_blob')
            if expected and digest(p, git_blob=True) != expected:
                raise PlanError(f'Source changed; review and update binding explicitly: {s["path"]}')
    tasks = plan.get('tasks')
    if not isinstance(tasks, list) or not tasks:
        raise PlanError('Nonempty tasks required')
    by_id: dict[str, dict] = {}
    for t in tasks:
        tid = t.get('id')
        nonempty(tid, 'task id')
        if tid in by_id:
            raise PlanError(f'Duplicate task: {tid}')
        by_id[tid] = t
        if t.get('kind') not in KINDS or t.get('status') not in STATES:
            raise PlanError(f'Invalid kind/status: {tid}')
        if type(t.get('priority')) is not int or t['priority'] < 0:
            raise PlanError(f'Invalid priority: {tid}')
        if type(t.get('nonclosing_attempts')) is not int or t['nonclosing_attempts'] < 0:
            raise PlanError(f'Invalid attempt count: {tid}')
        for key in ('title', 'acceptance', 'first_step', 'fallback'):
            nonempty(t.get(key), f'{tid}.{key}')
        if not isinstance(t.get('depends_on'), list) or not isinstance(t.get('inputs'), list):
            raise PlanError(f'List dependencies/inputs required: {tid}')
        if not set(t['inputs']) <= source_ids:
            raise PlanError(f'Unknown input: {tid}')
        if t['status'] == 'done':
            for key in ('proof', 'review'):
                nonempty(t.get(key), f'{tid}.{key}')
                local_file(root, t[key])
        if t['kind'] != 'primary' and t['status'] not in {'suspended', 'conditional'}:
            nonempty(t.get('activation_reason'), f'{tid}.activation_reason')
            nonempty(t.get('dependency_for'), f'{tid}.dependency_for')
    def visit(tid: str, active: set[str], finished: set[str]) -> None:
        if tid in active:
            raise PlanError(f'Dependency cycle at {tid}')
        if tid in finished:
            return
        active.add(tid)
        for dep in by_id[tid]['depends_on']:
            if dep not in by_id:
                raise PlanError(f'Unknown dependency {dep} for {tid}')
            visit(dep, active, finished)
        active.remove(tid); finished.add(tid)
    finished: set[str] = set()
    for tid in by_id:
        visit(tid, set(), finished)
    for t in tasks:
        if t['status'] == 'done' and any(by_id[d]['status'] != 'done' for d in t['depends_on']):
            raise PlanError(f'Done task has unfinished prerequisite: {t["id"]}')
        if t.get('dependency_for') is not None:
            parent = by_id.get(t['dependency_for'])
            if parent is None or parent['kind'] != 'primary':
                raise PlanError(f'Activation must support a primary task: {t["id"]}')
    m = plan.get('mathematical_status', {})
    lo, r4 = m.get('global_lower_bound'), m.get('rank4_lower_bound')
    if type(lo) is not int or type(r4) is not int or not 0 <= r4 <= lo:
        raise PlanError('Invalid reported component lower bounds')
    upper = m.get('global_upper_bound')
    if upper is not None and (type(upper) is not int or upper < lo):
        raise PlanError('Invalid upper bound')
    if m.get('global_exhaustion'):
        if by_id.get('F', {}).get('status') != 'done' or upper is None:
            raise PlanError('Global exhaustion needs completed F and an explicit upper bound')
    if m.get('generic_five_exhaustion') and by_id.get('G', {}).get('status') != 'done':
        raise PlanError('Generic five-point claim needs completed G')
    overrides = []
    if check_inputs:
        for folder in ('', 'research', 'research/classification', 'scripts', 'papers', 'tests'):
            p = root / folder / 'AGENTS.override.md'
            if p.exists():
                overrides.append(str(p.relative_to(root)))
        if overrides:
            raise PlanError('Instruction override requires inspection: ' + ', '.join(overrides))
        # Root guidance is intentionally small; do not increase the Codex global cap.
        if local_file(root, 'AGENTS.md').stat().st_size > 12 * 1024:
            raise PlanError('Root AGENTS.md exceeds the local 12 KiB design budget')
    return by_id

def available(by_id: dict[str, dict]) -> list[dict]:
    return sorted((t for t in by_id.values()
                   if t['status'] in {'ready', 'blocked'}
                   and all(by_id[d]['status'] == 'done' for d in t['depends_on'])),
                  key=lambda t: (t['priority'], t['id']))

def validate_job(job: dict, by_id: dict[str, dict], root: Path) -> None:
    tid = job.get('task_id')
    eligible = {t['id'] for t in available(by_id)}
    if tid not in eligible:
        raise PlanError(f'Task is not active/ready: {tid}')
    for key in ('target_statement', 'implication_to_classification', 'coefficient_ring',
                'hypotheses', 'success_certificate', 'independent_check', 'negative_result_next_step'):
        nonempty(job.get(key), key)
    for key in ('inverted_elements', 'complements'):
        if not isinstance(job.get(key), list):
            raise PlanError(f'Explicit list required: {key}')
        for v in job[key]:
            nonempty(v, key)
    if job['inverted_elements'] and not job['complements']:
        raise PlanError('Localizations require a complement disposition')
    inputs = job.get('exact_inputs_and_hashes')
    if not isinstance(inputs, list) or not inputs:
        raise PlanError('Bound exact inputs required')
    for binding in inputs:
        p = local_file(root, binding.get('path'))
        if binding.get('sha256') != digest(p):
            raise PlanError(f'Input SHA256 mismatch: {binding.get("path")}')
    limit = job.get('resource_limit', {})
    for key in ('timeout_seconds', 'memory_mb'):
        if type(limit.get(key)) is not int or limit[key] <= 0:
            raise PlanError(f'Positive integer resource limit required: {key}')
    if by_id[tid]['nonclosing_attempts'] >= 2:
        nonempty(job.get('route_change_reason'), 'route_change_reason after two nonclosing attempts')

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    sub = ap.add_subparsers(dest='command', required=True)
    c = sub.add_parser('check', help='Check plan, graph, source bindings and evidence paths')
    c.add_argument('--no-inputs', action='store_true', help='Package-only check; NOT a repository preflight')
    n = sub.add_parser('next', help='Print the next eligible task; execute nothing')
    n.add_argument('--count', type=int, default=1)
    g = sub.add_parser('gate', help='Check a bounded research-job contract; execute nothing')
    g.add_argument('job', type=Path)
    args = ap.parse_args(argv)
    try:
        root = args.root.resolve()
        plan = read_json(root / PLAN_PATH)
        sources = read_json(local_file(root, plan['source_registry']))
        checks = not (args.command == 'check' and args.no_inputs)
        tasks = validate(plan, sources, root, check_inputs=checks)
        if args.command == 'check':
            print('PASS: plan schema, dependencies and declared scope are consistent.')
            print('Repository inputs checked.' if checks else 'Repository inputs NOT CHECKED (package-only mode).')
            print('This is NOT a verification of any mathematical theorem.')
        elif args.command == 'next':
            if args.count < 1:
                raise PlanError('--count must be positive')
            for t in available(tasks)[:args.count]:
                print(f'{t["id"]}: {t["title"]}\n  First: {t["first_step"]}\n  Accept: {t["acceptance"]}\n  Fallback: {t["fallback"]}')
                if t['nonclosing_attempts'] >= 2:
                    print('  ROUTE CHANGE REQUIRED before another job.')
            if not available(tasks):
                print('No ready task. Inspect blockers; do not auto-resume suspended secondary work.')
        else:
            validate_job(read_json(args.job), tasks, root)
            print('PASS: job contract and input bindings. Mathematics NOT verified. No job was started.')
        return 0
    except (PlanError, KeyError, TypeError, OSError) as exc:
        print(f'FAIL: {exc}', file=sys.stderr)
        return 2

if __name__ == '__main__':
    raise SystemExit(main())

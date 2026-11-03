#!/usr/bin/env python3
"""Freeze bounded new generic projection attempts, without solving claims."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / 'work/affine-infinity-lift-v62'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    inputs = [Path(__file__)]
    dependencies = ['S61_AFFINE_FINITE_WRONSKI_CRITERION_20261103_v62',
                    'S61_AFFINE_MONIC_DIMENSION_CRITERION_20261103_v62']
    for dependency in dependencies:
        directory = ROOT / 'reports' / dependency
        report_path = directory / 'report.json'
        saved = json.loads(report_path.read_text())
        assert saved['all_passed']
        inputs.append(report_path)
        for name, value in saved['input_sha256'].items():
            assert sha(ROOT / name) == value
        for name, value in saved['evidence_sha256'].items():
            assert sha(directory / name) == value
    specs = [
        ('algebraic87-polynomial-unit', 'm2', 'M2 --script',
         'QQ(ss1,...,ss5,cc1,cc2)[zz,q45]/(zz^87-2)', 46, 137, 124),
        ('monic-six-exact', 'sing', 'run-singular.sh --cpus=0 -q',
         'QQ(ss1,...,ss5,cc)[q44], q32 eliminated exactly', 44, 134, 14),
        ('monic-polynomial-cc', 'sing', 'run-singular.sh --cpus=0 -q',
         'QQ(ss1,...,ss5)[q45,cc], block elimination q45 before cc', 46, 135, 124),
    ]
    attempts = []
    for stem, suffix, executable, coefficient_ring, variables, equations, expected in specs:
        source = WORK / (stem + '.' + suffix)
        log = WORK / (stem + '.log')
        receipt = WORK / (stem + '-exit.json')
        value = json.loads(receipt.read_text())
        assert value['exit_code'] == expected
        contents = log.read_text()
        assert '_START' in contents
        assert '_GB_DONE' not in contents and '_COFACTOR_PASS' not in contents
        if expected == 14:
            assert 'no more memory' in contents
        for path in [source, log, receipt]:
            target = args.output / (path.name + '.gz')
            target.write_bytes(gzip.compress(path.read_bytes(), mtime=0))
        attempts.append(dict(id=stem, source_sha256=sha(source), output_sha256=sha(log),
                             receipt_sha256=sha(receipt), executable=executable,
                             coefficient_ring=coefficient_ring, polynomial_unknowns=variables,
                             equations=equations, **value,
                             result='memory_limit' if expected == 14 else 'timeout',
                             exact_unit_or_eliminant_obtained=False))
    # The five-parameter polynomial-domain version preserves every original row.
    full = (ROOT / 'reports/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62/wronski-exact.sing').read_text()
    original = full.split('ideal I=', 1)[1].split(';\n', 1)[0].split(',')[:134]
    polynomial = (WORK / 'monic-polynomial-cc.sing').read_text()
    rows = polynomial.split('ideal I=', 1)[1].split(';\n', 1)[0].split(',')
    assert len(rows) == 135 and rows[:134] == original
    assert 'diff(G[j],var(i))!=0' in polynomial
    # Ordinary87 and monic6 runs use immutable complete producer inputs.
    assert (WORK / 'algebraic87-polynomial-unit.m2').read_bytes() == (
        ROOT / 'reports/S61_AFFINE_FINITE_WRONSKI_CRITERION_20261103_v62/algebraic87-polynomial-exact.m2').read_bytes()
    assert (WORK / 'monic-six-exact.sing').read_bytes() == (
        ROOT / 'reports/S61_AFFINE_MONIC_DIMENSION_CRITERION_20261103_v62/monic-six-exact.sing').read_bytes()
    inputs.append(ROOT / 'research/notes/S61_AFFINE_NEW_PROJECTION_ATTEMPTS_20261103_v62.md')
    report = dict(id='S61_AFFINE_NEW_PROJECTION_ATTEMPTS_20261103_v62',
                  all_receipts_verified=True, checked_at_utc=datetime.now(timezone.utc).isoformat(),
                  attempts=attempts, setup_repeated=False, network_required=False,
                  full_affine_dimension_or_length_result_obtained=False,
                  blocker='Exact parameter-dependent elimination exceeds the bounded memory or runtime budgets; no identical computation is automatically repeated.',
                  restart_requirements='A new structural reduction or a justified alternative exact algorithm; otherwise additional explicit compute budget is needed for these same models.',
                  limitation='This report verifies unsuccessful bounded computations, not the mathematical unit conclusion.',
                  input_sha256={str(path.relative_to(ROOT)): sha(path) for path in inputs},
                  evidence_sha256={path.name: sha(path) for path in args.output.iterdir() if path.is_file()},
                  seconds=time.monotonic() - started)
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS receipts only:', [(a['id'], a['result']) for a in attempts])


if __name__ == '__main__':
    main()

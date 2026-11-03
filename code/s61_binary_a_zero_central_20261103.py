#!/usr/bin/env python3
"""Rebuild the marked a=0 central algebras and certify their reduced supports.

Uses the existing full135 extractor, retaining both positive beta coordinates.
Macaulay2 is selected through PATH; source cloud-setup/activate.sh on Cloud.
No historical input/output is overwritten. A nonzero exit is a failed check.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

import sympy as sp
from s61_binary_a_zero_atlas_20261102 import transversal
from s61_binary_endpoint_full_transverse_20261102 import CUTS, load_report, ROOT


SLOPES = ['2', '5/3', '5/4', '-10', '5/2', '25/14', '5', 'd0']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def m2_scalar(value, prime):
    value = sp.Rational(value)
    if prime:
        return str(int(value.p) * pow(int(value.q), -1, prime) % prime)
    return f'({value.p}/{value.q})'


def case_source(slope, points, prime):
    if slope == 'd0':
        b, d, e = 6, 0, 16
    else:
        ratio = sp.Rational(slope)
        b = 7 / (ratio - sp.Rational(8, 3))
        d, e = 1, b * ratio
    rows, weights, grades = transversal(b, 1, d, e)
    variables = [i for i in range(47) if i not in CUTS]
    base = f'(ZZ/{prime})' if prime else 'QQ'
    source = (f'R={base}[' + ','.join(f'q{i}' for i in variables)
              + ',Degrees=>{' + ','.join(str(weights[i]) for i in variables)
              + '}];\n')
    def polynomial(row):
        return '+'.join(m2_scalar(c, prime) + ''.join(f'*q{i}' for i in mon)
                        for mon, c in row.items()) or '0'
    source += 'I=ideal(' + ','.join(map(polynomial, rows)) + ');\n'
    chart = load_report('S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101')['linear']
    ideals = []
    for k, point in enumerate(pt for pt in points if pt['slope'] == slope):
        q = [sp.S.Zero] * 45
        beta = [sp.S.Zero] * 16
        for i, value in point['q'].items():
            q[int(i)] = sp.Rational(value)
        for i, value in point['beta'].items():
            beta[int(i)] = sp.Rational(value)
        assert q[chart['fixed']] == 1
        central = [q[i] for i in chart['cols']] + [0, beta[2], beta[3]]
        assert not any(central[i] for i in CUTS)
        anchor = min((i for i in variables if central[i]), key=lambda i: weights[i])
        equations = []
        for i in variables:
            if i == anchor:
                continue
            if not central[i]:
                equations.append(f'q{i}')
            else:
                power, remainder = divmod(weights[i], weights[anchor])
                assert remainder == 0
                equations.append(f'q{i}-' + m2_scalar(central[i] / central[anchor]**power, prime)
                                 + f'*q{anchor}^{power}')
        source += f'K{k}=ideal(' + ','.join(equations) + ');\n'
        # Each Kk is a polynomial ray; verify the original equations on it.
        source += f'assert isSubset(I,K{k});\n'
        ideals.append(f'K{k}')
    support = ('intersect(' + ','.join(ideals) + ')' if len(ideals) > 1
               else ideals[0] if ideals else 'ideal vars R')
    source += f'K={support}; G=gb I;\n'
    # K is radical (intersection of explicit polynomial rays, or the origin).
    # Thus these two containments certify sqrt(I)=K without guessing from HF.
    source += ('assert isSubset(I,K);\n'
               'assert all(flatten entries gens K,f->any(1..8,n->f^n%G==0));\n'
               'print("DIM="|toString dim(R/I));\n'
               'print("HS="|toString hilbertSeries(R/I,Reduce=>true));\n')
    if not ideals:
        length = 70 if slope == '5' else 65
        # degree is not vector-space length in this weighted grading.
        source += f'assert(dim(R/I)==0); assert(numColumns basis(R/I)=={length});\n'
    else:
        source += 'assert(dim(R/I)==1);\n'
    source += 'print "CENTRAL SUPPORT PASS"; exit 0;\n'
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--prime', type=int, default=32003, help='0 means QQ')
    parser.add_argument('--timeout', type=int, default=120)
    parser.add_argument('--slopes', nargs='+', choices=SLOPES, default=SLOPES)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    executable = shutil.which('M2')
    if not executable:
        raise RuntimeError('M2 missing from PATH; activate the existing Cloud tool prefix')
    points_path = ROOT / 'code/certificates/s61_binary_a_zero_points_20261102.json'
    points = json.loads(points_path.read_text())
    records = []
    for slope in args.slopes:
        name = slope.replace('/', '_')
        script, log = args.output / f'{name}.m2', args.output / f'{name}.log'
        script.write_text(case_source(slope, points, args.prime))
        started = time.monotonic()
        with log.open('w') as stream:
            result = subprocess.run([executable, '--script', str(script.resolve())],
                                    stdout=stream, stderr=subprocess.STDOUT,
                                    timeout=args.timeout, check=True)
        content = log.read_text()
        assert 'CENTRAL SUPPORT PASS' in content
        record = dict(slope=slope, prime=args.prime, seconds=time.monotonic()-started,
                      script_sha256=digest(script), log_sha256=digest(log),
                      output=content, returncode=result.returncode)
        records.append(record)
        print(slope, content.strip(), flush=True)
    inputs = [Path(__file__), points_path,
              ROOT/'code/s61_binary_a_zero_atlas_20261102.py',
              ROOT/'code/s61_binary_endpoint_full_transverse_20261102.py']
    report = dict(all_passed=True, coefficient_field='QQ' if not args.prime else f'GF({args.prime})',
                  m2_executable=executable, cases=records,
                  input_sha256={str(p.relative_to(ROOT)): digest(p) for p in inputs})
    (args.output / 'report.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()

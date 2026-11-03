#!/usr/bin/env python3
"""Reproduce Cloud a=0 certificates, preserving full inputs and checked outputs."""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import sympy as sp
from s61_binary_a_zero_atlas_20261102 import transversal, reduce_line
from s61_binary_a_zero_ray_transverse_20261102 import all_rays
from s61_binary_endpoint_full_transverse_20261102 import ROOT, PRIME, B, CUTS
from s61_binary_zero_pair_stratum_removal_20261102 import write_quartics, compile_span


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command, log, timeout):
    with log.open('w') as stream:
        subprocess.run(list(map(str, command)), stdout=stream, stderr=subprocess.STDOUT,
                       timeout=timeout, check=True)
    return log.read_text()


def infinity_source():
    # The compactification b=u^2, x=u^2,z=u^-1 has endpoint
    # (a,b,c,d,e,f)=(0,1,0,1,8/3,0), not the separate d=0 point.
    rows, weights, _ = transversal(1, 0, 1, sp.Rational(8, 3))
    variables = [i for i in range(47) if i not in CUTS]
    def polynomial(row):
        return '+'.join(f'({c})'+''.join(f'*q{i}' for i in mon)
                        for mon, c in row.items()) or '0'
    return ('R=QQ['+','.join(f'q{i}' for i in variables)+',Degrees=>{'
            +','.join(str(weights[i]) for i in variables)+'}];\n'
            +'I=ideal('+','.join(map(polynomial, rows))+');\n'
            +'assert(dim(R/I)==0); print("LENGTH="|toString numColumns basis(R/I));\n'
            +'print("HS="|toString hilbertSeries(R/I,Reduce=>true));\n'
            +'print "INFINITY FINITE PASS"; exit 0;\n')


def verify(output):
    start = time.monotonic()
    output.mkdir(parents=True, exist_ok=False)
    central = ROOT/'code/s61_binary_a_zero_central_20261103.py'
    for prime, name in [(0, 'central_QQ'), (PRIME, 'central_modp')]:
        run([sys.executable, central, '--prime', prime, '--output', output/name],
            output/f'{name}.log', 300)
    raw, weights, grades = transversal(B, 1, 1, 8*B/3+7)
    residuals, free_weights, determinants = reduce_line(raw, weights, grades)
    row_input = output/'weight4.in'
    nrows = write_quartics(row_input, residuals, free_weights, grades)
    executable = output/'row_determinant'
    source = ROOT/'code/s61_binary_a_zero_row_determinant_20261102.cpp'
    compile_span(source, executable)
    values = list(map(int, run([executable, row_input], output/'determinant.log', 90).split()))
    nr, nc, operations, degree = values[:4]
    assert nr == nrows == 531 and nc == 386 and degree == 17
    determinant = sp.Poly.from_list(list(reversed(values[4:])), B, modulus=PRIME)
    assert determinant.degree() == degree
    roots = {int(-factor.nth(0)/factor.nth(1)) % PRIME: multiplicity
             for factor, multiplicity in determinant.factor_list()[1]}
    assert roots == {857:1,15991:1,31996:1,28233:2,2526:2,3:3,31961:7}
    heavy = [i for i, weight in enumerate(free_weights) if weight == 3]
    def coefficient(polynomial, index):
        monomial = tuple(int(i == index) for i in range(len(free_weights)))
        return sp.Poly.from_dict({(mon[0],):int(c) for mon,c in polynomial.items()
                                  if mon[1:] == monomial}, B, modulus=PRIME).as_expr()
    matrix = sp.Matrix([[coefficient(row,i) for i in heavy]
                        for row,grade in zip(residuals,grades) if grade == 3 and row])
    gcd = sp.Poly(0, B, modulus=PRIME)
    for rows in itertools.combinations(range(matrix.rows), 3):
        gcd = sp.gcd(gcd, sp.Poly(matrix[list(rows),:].det(), B, modulus=PRIME))
        if gcd.degree() == 0:
            break
    expected_gcd = sp.Poly((B-857)*(B-3)**2, B, modulus=PRIME)
    assert gcd.monic() == expected_gcd.monic()
    (output/'weight3_matrix.json').write_text(json.dumps(
        dict(matrix=[[str(x) for x in row] for row in matrix.tolist()],
             gcd=str(gcd.monic().as_expr())), indent=2)+'\n')
    points_path = ROOT/'code/certificates/s61_binary_a_zero_points_20261102.json'
    rays = all_rays(json.loads(points_path.read_text()))
    (output/'rays.json').write_text(json.dumps(rays, indent=2)+'\n')
    script = output/'infinity_QQ.m2'
    script.write_text(infinity_source())
    m2 = shutil.which('M2')
    assert m2
    infinity = run([m2, '--script', script.resolve()], output/'infinity_QQ.log', 90)
    assert 'INFINITY FINITE PASS' in infinity
    inputs = [Path(__file__), central, source, points_path,
              ROOT/'code/s61_binary_a_zero_atlas_20261102.py',
              ROOT/'code/s61_binary_a_zero_ray_transverse_20261102.py',
              ROOT/'code/s61_binary_endpoint_full_transverse_20261102.py',
              ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
              ROOT/'code/s61_model.py', ROOT/'code/toric_closed_jet_coordinates_20261102.py',
              ROOT/'code/weight1_structure.py',
              ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json']
    for stem in ['S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102',
                 'S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101',
                 'S61_UNIVERSAL_WEIGHT6_20261031']:
        inputs.append(ROOT/'reports'/f'{stem}.json')
    # The executable is reproducible from the preserved source; do not publish it.
    executable.unlink()
    report = dict(id='S61_BINARY_A_ZERO_CLOUD_20261103', all_passed=True,
                  start_commit='4c3f623fb698f9adc6077015b25e2325b11072af',
                  m2_version=subprocess.run([m2,'--version'], capture_output=True,text=True,
                                            check=True,timeout=10).stdout.strip(),
                  shape=[nr,nc], row_operations=operations,
                  determinant_degree_mod_32003=degree,
                  determinant_coefficients_low_to_high_mod_32003=values[4:],
                  determinant_roots_mod_32003=roots,
                  heavy_cubic_gcd_mod_32003=str(gcd.monic().as_expr()),
                  constant_linear_pivot_determinants=determinants,
                  ray_transverse_lengths_upper_QQ=[r['transverse_length'] for r in rays],
                  infinity_finite_QQ=True,
                  limitations=['Modular determinant factors do not identify the QQ determinant.',
                               'Ambient nilpotence and whole-stratum removal require the proof note.'],
                  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},
                  evidence_sha256={str(p.relative_to(output)):sha(p)
                                   for p in sorted(output.rglob('*')) if p.is_file()},
                  seconds=time.monotonic()-start)
    (output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.output.resolve())
    print('PASS', result['shape'], result['ray_transverse_lengths_upper_QQ'], result['seconds'])

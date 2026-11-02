#!/usr/bin/env python3
"""Exact no-divisor rank-drop certificate for the 42 cubic products on K=0."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import sympy as s
from sympy.polys.rings import ring

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT/'research/notes/S61_BINARY_HYPERPLANE_PRODUCTS_CODIMENSION_20261102.md'
SOURCE = ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json'
PAIR = ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_PAIR_20261102.json'
CANDIDATE = ROOT/'code/certificates/s61_binary_hyperplane_products_20261102.json'
NATIVE = ROOT/'code/s61_binary_primitive_minor_support_20261102.cpp'


def verify():
    start = time.time()
    data = json.loads(SOURCE.read_text())
    pair = json.loads(PAIR.read_text())
    for report in [data, pair]:
        assert report['all_passed']
        for path, expected in report['input_sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected
    assert pair['generic_product_rank'] == 41
    candidate = json.loads(CANDIDATE.read_text())
    R, a, e, d = ring('a,e,d', s.QQ)
    K = R.to_field()
    known = [a, 500*a+90*e**2+3*e-24,
             50*a*d+18*d*e**2-9*d*e-2*d+60*e-40,
             100*a+18*e**2-9*e-2]
    def strip(poly):
        assert poly
        for factor in known:
            while poly.rem(factor) == 0:
                poly = poly.exquo(factor)
        return poly
    def unit(poly):
        return bool(poly) and max(map(sum, strip(poly))) == 0
    qs = [{tuple(mon): R.from_dict({tuple(mm): s.QQ(c) for mm, c in terms})
           for mon, terms in row} for row in data['quadratic_rows_sparse']]
    products = []
    for q in qs:
        for j in [0, 1, 2, 3, 4, 6, 7]:
            row = {}
            for mon, poly in q.items():
                shifted = list(mon)
                shifted[j] += 1
                row[tuple(shifted)] = poly
            products.append(row)
    mons = sorted(set().union(*(row.keys() for row in products)))
    assert len(products) == 42 and len(mons) == 84
    M = [[K(row.get(mon, R.zero)) for mon in mons] for row in products]
    used = []
    for rr, (i, j, expected) in enumerate(candidate['unit_pivots']):
        assert i >= rr and j not in used
        M[rr], M[i] = M[i], M[rr]
        pivot = M[rr][j]
        assert pivot == K.from_expr(s.sympify(expected))
        assert unit(pivot.numer) and unit(pivot.denom)
        M[rr] = [x/pivot for x in M[rr]]
        for i in range(rr+1, len(M)):
            coefficient = M[i][j]
            if coefficient:
                M[i] = [x-coefficient*y for x, y in zip(M[i], M[rr])]
        used.append(j)
    assert len(used) == 21
    cols = [j for j in range(84) if j not in used]
    residual = []
    for row in M[21:]:
        den = R.one
        for j in cols:
            den = den.lcm(row[j].denom)
        assert unit(den)
        cleared = [row[j].numer*den.exquo(row[j].denom) for j in cols]
        content = R.zero
        for poly in cleared:
            content = content.gcd(poly)
        if content:
            known_part = content.exquo(strip(content))
            assert unit(known_part)
            cleared = [poly.exquo(known_part) for poly in cleared]
        residual.append(cleared)
    assert len(residual) == 21 and len(cols) == 63
    def write_poly(stream, poly):
        stream.write(str(len(poly)))
        for mon, coefficient in poly.items():
            stream.write(' '+str(coefficient)+' '+' '.join(map(str, mon)))
        stream.write('\n')
    with tempfile.TemporaryDirectory(prefix='poisson-hyperplane-products-') as temp:
        work = Path(temp)
        inp, out, binary = work/'input', work/'output', work/'native'
        with inp.open('w') as stream:
            stream.write('21 63 20\n')
            for row in residual:
                for poly in row:
                    write_poly(stream, poly)
            stream.write('2\n')
            for rows, columns in candidate['selected_minors']:
                assert len(rows) == len(columns) == 20
                stream.write(' '.join(map(str, rows+columns))+'\n')
        command = ['clang++', '-O3', '-std=c++17',
                   '-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),
                   '-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),
                   '-I/opt/homebrew/include', str(NATIVE),
                   str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),
                   '-L/opt/homebrew/lib', '-lgmp', '-o', str(binary)]
        subprocess.run(command, check=True, capture_output=True, text=True)
        env = os.environ.copy()
        env['DYLD_LIBRARY_PATH'] = str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'
        run = subprocess.run([str(binary), str(inp), str(out), '2'],
                             env=env, check=True, capture_output=True, text=True, timeout=120)
        output_sha = hashlib.sha256(out.read_bytes()).hexdigest()
        def read_poly(line):
            tokens = line.split()
            n = int(tokens[0])
            assert len(tokens) == 1+4*n
            return R.from_dict({tuple(map(int, tokens[2+4*k:5+4*k])):s.QQ(tokens[1+4*k])
                                for k in range(n)})
        factors, cross = [], []
        with out.open() as stream:
            for line in stream:
                words = line.split()
                if words[0] == 'MINOR':
                    factors.append([])
                elif words[0] == 'EXP':
                    exponent = int(words[1])
                    assert exponent > 0
                    poly = read_poly(stream.readline())
                    factors[-1].append({'exponent': exponent,
                                        'degree': max(map(sum, poly)), 'terms': len(poly)})
                else:
                    assert words == ['CROSS_GCD']
                    poly = read_poly(stream.readline())
                    assert unit(poly)
                    cross.append(str(poly.as_expr()))
        assert len(factors) == 2 and len(cross) > 0
    paths = [Path(__file__), NOTE, SOURCE, PAIR, CANDIDATE, NATIVE]
    return {'id':'S61_BINARY_HYPERPLANE_PRODUCTS_CODIMENSION_20261102',
            'all_passed':True, 'product_matrix_shape':[42,84],
            'unit_pivot_count':21, 'residual_matrix_shape':[21,63],
            'minor_size':20, 'selected_minors':candidate['selected_minors'],
            'minor_factor_sizes':factors, 'cross_gcd_polynomials':cross,
            'all_cross_gcd_are_a_F_B_C_units':True,
            'native_checks':{'all_content_divisions_exact':True,
                             'all_factorizations_reexpanded':True,
                             'final_factor_exponents_nonnegative':True,
                             'all_cross_factor_gcd_checked':True},
            'native_output_sha256':output_sha, 'native_log':run.stderr,
            'generic_rank':41, 'bad_support_dimension_upper_on_three_core_open':1,
            'limitations':['Only nuisance-zero on a*F*B*C nonzero.',
                           'No assertion of rank41 at all points.',
                           'No complete cubic/quartic rank or ambient removal.'],
            'proof_note':str(NOTE.relative_to(ROOT)),
            'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                            for path in paths}, 'seconds':time.time()-start}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    result = verify()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({key:result[key] for key in ['all_passed', 'generic_rank',
                      'bad_support_dimension_upper_on_three_core_open', 'seconds']}))

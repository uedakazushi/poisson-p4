#!/usr/bin/env python3
"""Lift the 274 AB quadrics to the actual T01 ideal using polynomial associativity."""
import argparse
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import combinations_with_replacement, product
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time

import sympy as sp

from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from toric_algebra_fiber import TRIPLES, associativity_quadrics, multiplication_from_jet
from weight1_structure import COORDS

ROOT = Path(__file__).resolve().parents[1]
STEM = 'TORIC_ZERO_JET_ASSOCIATIVE_LIFTS_20261102'
PLANE = 'TORIC_UNIVERSAL_THREE_PLANE_20261101'
QUADRICS = 'AB_EXPLICIT_QUADRICS_CUBIC_GENERATION_20261102'
NATIVE = ROOT / 'code' / (STEM.lower() + '.cpp')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_report(name):
    report = json.loads((ROOT/'reports'/(name+'.json')).read_text())
    assert report['all_passed']
    return report


def accumulate(poly, monomial, value):
    result = poly.get(monomial, 0) + value
    if result:
        poly[monomial] = result
    else:
        poly.pop(monomial, None)


def multiply_monomials(left, right):
    return tuple(sorted(left + right))


def coordinate_data(plane):
    model = S61()
    y, beta, form = primitive_beta(model)
    cubic = list(sp.symbols('c0:20'))
    variables = list(beta) + cubic
    assert list(map(str, beta)) == plane['beta_variables']
    assert list(map(str, cubic)) == plane['c_variables']
    for weight in (1, 2):
        actual = [i for i,w in enumerate(model.reduced['parameter_weights']) if w == weight]
        assert actual == plane[f'weight{weight}_u_indices_zero_based']
    tensor = sp.Matrix([form[i,j].diff(y[k]) for (i,j),k in COORDS])
    factors = tensor.jacobian(beta)
    standard = [(1,0), (-1,0), (0,1), (0,-1)]
    def covariant_weight(indices):
        return tuple(-sum(standard[j][k] for j in indices) for k in (0,1))
    weights = []
    for h in range(len(beta)):
        w = {covariant_weight((i,j,k)) for row,((i,j),k) in enumerate(COORDS)
             if factors[row,h] != 0}
        assert len(w) == 1
        weights.append(w.pop())
    weights.extend(covariant_weight(triple) for triple in TRIPLES)

    @lru_cache(None)
    def polynomial(expression):
        result = {}
        for powers, coefficient in sp.Poly(expression,*variables,domain=sp.QQ).terms():
            if coefficient:
                term = tuple(i for i,power in enumerate(powers) for _ in range(power))
                result[term] = Fraction(int(coefficient.p), int(coefficient.q))
        return result

    def bilinear_rows(matrix):
        result = []
        for values in matrix:
            row = {}
            assert len(values) == 20
            for j, entry in enumerate(values):
                if entry == '0':
                    continue
                for term, coeff in polynomial(sp.sympify(entry)).items():
                    assert len(term) == 1 and term[0] < 16
                    row[(term[0],16+j)] = coeff
            result.append(row)
        return result

    multiplications = multiplication_from_jet(tensor, dict(zip(TRIPLES,cubic,strict=True)))
    _, equations = associativity_quadrics(multiplications)
    E = [polynomial(f) for f in equations]
    N = bilinear_rows(plane['q_linear_matrix_N'])
    L = bilinear_rows(plane['all51_homogeneous_linear_consequences'])
    centers = [polynomial(sp.sympify(f)) for f in plane['q_center']]
    assert len(E) == 96 and len(L) == 51 and len(N) == len(centers) == 45
    assert all(len(m) == 2 and max(m) < 16 for q in centers for m in q)
    for expected, combination in zip(L,plane['all51_constant_combinations_of96'],strict=True):
        row = {}
        for j,c in combination:
            for monomial, value in E[j].items():
                accumulate(row,monomial,Fraction(c)*value)
        assert row == expected
    return weights, E, L, N


def weight_of(monomial, weights):
    return tuple(sum(weights[j][k] for j in monomial) for k in (0,1))


def weight_blocks(E, L, N, ab, weights):
    def unique_weight(poly):
        values = {weight_of(m,weights) for m in poly}
        assert len(values) == 1
        return values.pop()

    assert [unique_weight(q) for q in N] == list(map(tuple,ab['coordinate_weights']))
    targets = defaultdict(list)
    for index, relation in enumerate(ab['quadratic_relations']):
        row = {}
        for j,k,coeff in relation:
            scalar = Fraction(coeff)
            for left,x in N[j].items():
                for right,y in N[k].items():
                    accumulate(row,multiply_monomials(left,right),scalar*x*y)
        targets[unique_weight(row)].append((str(index),row))

    generators = defaultdict(list)
    pairs = list(combinations_with_replacement(range(16),2))
    mixed = list(product(range(16), range(16,36)))
    for tag, source, multipliers in [('E',E,pairs), ('L',L,mixed), ('Q',N,pairs)]:
        for i, polynomial in enumerate(source):
            if not polynomial:
                continue
            w = unique_weight(polynomial)
            for j,k in multipliers:
                v = weight_of((j,k),weights)
                total_weight = (w[0]+v[0],w[1]+v[1])
                if total_weight not in targets:
                    continue
                projected = {multiply_monomials(m,(j,k)):value for m,value in polynomial.items()
                             if max(k,*m) >= 16}
                if projected:
                    generators[total_weight].append((f'{tag}{i}.{j}.{k}',projected))
    return generators, targets


def write_matrix(path, left, right):
    columns = sorted({m for label,row in left+right for m in row})
    assert all(len(m) == 4 and sum(j>=16 for j in m) in (1,2) for m in columns)
    column_index = {m:j for j,m in enumerate(columns)}
    with path.open('w') as out:
        print(len(left),len(right),len(columns),file=out)
        for label,row in left+right:
            entries = ' '.join(f'{column_index[m]} {v.numerator} {v.denominator}'
                               for m,v in row.items())
            print(label,len(row),entries,file=out)
    return len(columns)


def compile_solver(directory):
    # Use the repository's pinned FLINT distribution, as in the AB degree-four certificate.
    library = ROOT/'work/m2-runtime/lib'
    executable = directory/'exact'
    subprocess.run(['c++','-O3','-std=c++17',
                    '-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),
                    '-I/opt/homebrew/include',str(NATIVE),str(library/'libflint.24.0.dylib'),
                    '-L/opt/homebrew/lib','-lgmp','-o',str(executable)],check=True,timeout=60)
    env = dict(os.environ, DYLD_LIBRARY_PATH=str(library)+':/opt/homebrew/lib')
    return executable, env


def rational_terms(poly):
    return [[*monomial,str(coefficient)] for monomial,coefficient in sorted(poly.items())]


def read_verified_lifts(path, left, right, E, lift_by_index):
    matrix = dict(left)
    count_by_type = dict(E=0,L=0,Q=0)
    with path.open() as stream:
        header = stream.readline().split()
        assert header[0] == 'EXACT' and header[4] == 'PIVOTS'
        assert list(map(int,header[1:3])) == [len(left),len(right)]
        for label, target in right:
            title = stream.readline().split()
            assert title[:2] == ['TARGET',label] and int(label) not in lift_by_index
            residual = dict(target)
            linear = {}
            constant = {}
            for _ in range(int(title[2])):
                name, coeff = stream.readline().split()
                scalar = Fraction(coeff)
                # Recheck the emitted coefficients, after the native code undoes row scalings.
                for monomial, value in matrix[name].items():
                    accumulate(residual,monomial,-scalar*value)
                tag = name[0]
                i,j,k = map(int,name[1:].split('.'))
                count_by_type[tag] += 1
                if tag == 'E':
                    for m,value in E[i].items():
                        if max(m) < 16:
                            accumulate(constant,multiply_monomials((j,k),m),scalar*value)
                if tag == 'Q':
                    accumulate(linear,(i,j,k),-scalar)
            assert residual == {}
            lift_by_index[int(label)] = dict(
                quadratic_index=int(label), beta2_z_terms=rational_terms(linear),
                beta4_terms=rational_terms(constant))
        assert not stream.read().strip()
    return int(header[5]),count_by_type


def run():
    start = time.monotonic()
    ab = load_report(QUADRICS)
    plane = load_report(PLANE)
    assert ab['coordinate_order'] == 'The ordered weight-two parameters of S61, q0..q44.'
    assert ab['quadratic_count'] == len(ab['quadratic_relations']) == 274
    weights, E, L, N = coordinate_data(plane)
    generators, targets = weight_blocks(E,L,N,ab,weights)
    lifts = {}
    blocks = []
    with tempfile.TemporaryDirectory(prefix='toric-associative-lifts-') as temp:
        wd = Path(temp)
        executable, env = compile_solver(wd)
        for weight in sorted(targets):
            left, right = generators[weight], targets[weight]
            source, witness = wd/'matrix',wd/'witness'
            columns = write_matrix(source,left,right)
            subprocess.run([str(executable),str(source),str(witness)],env=env,check=True,
                           capture_output=True,timeout=120)
            rank,counts = read_verified_lifts(witness,left,right,E,lifts)
            blocks.append(dict(weight=list(weight),generators=len(left),targets=len(right),
                               columns=columns,independent_rows_used=rank,witness_counts=counts,
                               matrix_sha256=digest(source),witness_sha256=digest(witness)))
    assert sorted(lifts) == list(range(274)) and len(blocks) == 29
    dependencies = [PLANE,'TORIC_UNIVERSAL_THREE_PLANE_REVIEW_20261101',
                    'TORIC_ALGEBRA_FIBER_20261031','TORIC_COFRAME_FIBER_20261031',
                    QUADRICS,'AB_QUADRICS_QUARTIC_MULTIPLICATION_20261102',
                    'AB_IDEAL_HIGH_GENERATOR_VANISHING_20261102','TORIC_ZERO_JET_EXHAUSTION_20261101']
    for name in dependencies:
        load_report(name)
    files = [Path(__file__),NATIVE,ROOT/'research/notes'/(STEM+'.md')]
    files.extend(ROOT/'reports'/(name+'.json') for name in dependencies)
    files.extend(ROOT/'research/notes'/(name+'.md') for name in dependencies)
    files.extend(ROOT/'code'/name for name in [
        's61_model.py','s61_infinity_flags.py','elliptic_weight1_dominance.py',
        'primitive_jet_elliptic_curve_20261101.py','toric_algebra_fiber.py',
        'toric_coframe_fiber.py','weight1_structure.py'])
    files.extend(ROOT/'legacy/v02/data'/name for name in [
        'symplectic_slice_certificate.json','symplectic_slice_equations.json',
        'symplectic_slice_reduced_equations.json'])
    return dict(id=STEM,all_passed=True,field='QQ',working_ring='QQ[b0..b15,c0..c19]',
                beta_variables=plane['beta_variables'],
                weight2_u_indices_zero_based=plane['weight2_u_indices_zero_based'],
                q_center=plane['q_center'],coordinate_weights=ab['coordinate_weights'],
                associative_equations=96,validated_linear_combinations=51,
                nonzero_linear_combinations=sum(bool(q) for q in L),
                exact_row_identities=274,monomial_blocks=blocks,
                lift_recipe='z=q-q_center(b). Use quadratic_relations[quadratic_index](z) + beta2_z_terms [i,j,k,v]: v*z_i*b_j*b_k + beta4_terms [i,j,k,l,v]: v*b_i*b_j*b_k*b_l. All indices zero-based.',
                associative_lifts=[lifts[i] for i in range(274)],
                actual_T01_zero_fiber_equals_reduced_AB_scheme=True,
                input_sha256={str(p.relative_to(ROOT)):digest(p) for p in files},
                seconds=time.monotonic()-start)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output already exists')
    report = run()
    with args.output.open('x') as stream:
        json.dump(report,stream,separators=(',',':'),ensure_ascii=False)
        stream.write('\n')
    print('PASS',report['exact_row_identities'],round(report['seconds'],3),flush=True)

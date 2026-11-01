#!/usr/bin/env python3
"""Generic seventh-jet upper bound from an exact good-prime rank certificate."""
import argparse
import hashlib
import json
import os
from itertools import combinations_with_replacement as cwr
from math import comb
from pathlib import Path
import subprocess
import tempfile
import time
import numpy as np
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from truncated_local_ideal_io_20261101 import write_input

ROOT = Path(__file__).resolve().parents[1]

def known_seventh(data):
    saved = {degree: json.loads((ROOT/('reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_'+degree+'_20261101.json')).read_text())
             for degree in ['QUADRATIC', 'CUBIC', 'QUARTIC', 'QUINTIC']}
    def polynomial(row, monomials):
        return '+'.join('('+c.replace('**', '^').replace('a1', 'a')+')*'
                        +'*'.join('v_'+str(i) for i in mon)
                        for c, mon in zip(row, monomials) if c != '0') or '0'
    polys = [polynomial(row, saved['QUADRATIC']['quadratic_monomials'])
             for row in saved['QUADRATIC']['constant_quadratic_rref']]
    v = s.symbols('v0:11')
    for expression in saved['CUBIC']['constant_cubic_polynomials']:
        poly = s.Poly(s.sympify(expression), *v)
        polys.append('+'.join('('+str(c)+')*'+'*'.join('v_'+str(i) for i, e in enumerate(mon) for _ in range(e))
                              for mon, c in poly.terms()))
    for label, rows, monomial_key in [('QUARTIC', 'new_quartic_rref', 'monomials'),
                                       ('QUINTIC', 'constant_fifth_rref', 'fifth_monomials')]:
        record = saved[label]
        monomials = [record[monomial_key][j] for j in record['quotient_free_columns']]
        polys.extend(polynomial(row, monomials) for row in record[rows])
    monomials = [data['sixth_monomials'][j] for j in data['quotient_free_columns']]
    sixth = [polynomial([str(coefficient) for column, coefficient in row],
                        [monomials[column] for column, coefficient in row])
             for row in data['constant_sixth_rref_sparse']]
    program = 'B=QQ[a];K=frac B;R=K[v_0..v_10];I=ideal('+','.join(polys)+');\n'
    program += 'assert(hilbertFunction(6,R/I)==1237);J=I+ideal('+','.join(sixth)+');\n'
    program += 'assert(apply(toList(0..6),d->hilbertFunction(d,R/J))=={1,11,54,168,391,740,1203});assert(hilbertFunction(7,R/J)==1797);print "KNOWN SEVENTH QUOTIENT PASS";\n'
    with tempfile.TemporaryDirectory(prefix='free-flag-sixth-review-') as directory:
        source = Path(directory)/'review.m2'
        source.write_text(program)
        env = os.environ.copy()
        env['DYLD_LIBRARY_PATH'] = str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
        process = subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary', '--script', str(source)],
                                 env=env, capture_output=True, text=True, check=True, timeout=60)
    assert 'KNOWN SEVENTH QUOTIENT PASS' in process.stdout
    return process.stdout+process.stderr


def run():
    started = time.time()
    prime = 32003
    assert s.isprime(prime)
    values = [1, 2, 3, 4, 5]
    parameters = s.symbols('a0:5')
    substitution = dict(zip(parameters, values))
    quadratic_path = ROOT/'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUADRATIC_20261101.json'
    normal_path = ROOT/'reports/S61_MINIMAL_FREE_FLAG_NORMAL_FORM_20261101.json'
    endpoint_path = ROOT/'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json'
    sextic_path = ROOT/'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_SEXTIC_20261101.json'
    quadratic = json.loads(quadratic_path.read_text())
    sextic = json.loads(sextic_path.read_text())
    assert quadratic['all_passed'] and sextic['all_passed']
    assert all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest for path,digest in sextic['input_sha256'].items())
    known_log = known_seventh(sextic)
    assert sextic['true_initial_hilbert_function'] == [1, 11, 54, 168, 391, 740, 1203]
    mat = lambda key: s.Matrix(quadratic[key]).applyfunc(s.sympify).subs(substitution)
    J, T, U = [mat(key) for key in ['jacobian', 'tangent_basis_after_shift', 'right_inverse']]
    T[:, 1] -= s.Rational(values[1], 6)*T[:, 0]
    rows = quadratic['pivot_rows']
    assert J*T == s.zeros(135, 11) and J[rows, :]*U == s.eye(34)
    def mod(value):
        value = s.Rational(value)
        assert int(value.q) % prime
        return int(value.p)*pow(int(value.q), -1, prime) % prime
    array = lambda M: np.array([[mod(value) for value in row] for row in M.tolist()], dtype=np.int64)
    j, u, tangent = [array(M) for M in [J, U, T]]
    JU = j@u % prime
    model = S61()
    one = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 1]
    beta = s.Matrix(json.loads(normal_path.read_text())['normalized_beta16']).applyfunc(s.sympify).subs(substitution)
    two, equations = specialize(model, dict(zip(one, beta)), True)
    point = json.loads(endpoint_path.read_text())['strata']['free_M0']['local_data']
    fixed, chart = two.index(point['chart_fixed_u']-1), point['chart_columns']
    H = np.zeros((135, 45, 45), dtype=np.int64)
    direct_J = s.zeros(135, 45)
    constants = [s.S.Zero]*135
    for row, equation in enumerate(equations):
        for monomial, coefficient in equation.items():
            multiplicity = monomial.count(fixed)
            if multiplicity == 2:
                constants[row] += coefficient
            elif multiplicity == 1:
                remaining = [i for i in monomial if i != fixed]
                column = chart.index(remaining[0]) if remaining else 44
                direct_J[row, column] -= coefficient
            else:
                i, k = [chart.index(v) for v in monomial]+[44]*(2-len(monomial))
                value = mod(coefficient)
                if i == k:
                    H[row, i, k] += value
                else:
                    half = value*pow(2, -1, prime) % prime
                    H[row, i, k] += half
                    H[row, k, i] += half
    assert direct_J == J and not any(constants)
    H %= prime
    monomials = {degree: list(cwr(range(11), degree)) for degree in range(1, 8)}
    indices = {degree: {mon: i for i, mon in enumerate(mons)} for degree, mons in monomials.items()}
    W, jets = {1: tangent}, {}
    for degree in range(2, 8):
        Q = np.zeros((135, len(monomials[degree])), dtype=np.int64)
        for da in range(1, degree//2+1):
            db = degree-da
            factor = 1 if da == db else 2
            target = np.array([indices[degree][tuple(sorted(x+y))]
                               for x in monomials[da] for y in monomials[db]])
            # All integer dot products are far below int64 overflow. The
            # floating bincount sums integers smaller than 2^53 exactly.
            assert 45*(prime-1)**2 < 2**63
            assert len(target)*(prime-1) < 2**53
            for row in range(135):
                cross = W[da].T@(H[row]@W[db] % prime) % prime
                accumulated = np.bincount(target, weights=cross.ravel(), minlength=len(monomials[degree])).astype(np.int64)
                Q[row] = (Q[row]+factor*accumulated) % prime
        W[degree] = -u@Q[rows, :] % prime
        jets[degree] = (Q-JU@Q[rows, :]) % prime
        assert not np.any(jets[degree][rows, :])
    sources = ['truncated_local_ideal_filtered_trace_20261101.cpp', 'truncated_local_standard_basis_20261101.cpp']
    with tempfile.TemporaryDirectory(prefix='minimal-seventh-upper-') as tmp:
        directory = Path(tmp)
        data = directory/'original-jets.txt'
        write_input(data, prime, 11, jets, maximum_degree=7)
        input_hash = hashlib.sha256(data.read_bytes()).hexdigest()
        binaries = []
        for source in sources:
            binary = directory/Path(source).stem
            subprocess.run(['c++', '-O3', '-std=c++17', str(ROOT/'code'/source), '-o', str(binary)], check=True, capture_output=True, timeout=60)
            binaries.append(binary)
        filtered = subprocess.run([str(binaries[0]), str(data), str(directory/'trace'), str(directory/'order'), str(directory/'initial')],
                                  text=True, capture_output=True, check=True, timeout=120)
        assert filtered.stdout.strip() == '27489'
        standard = subprocess.run([str(binaries[1]), str(data), str(directory/'standard')],
                                  text=True, capture_output=True, check=True, timeout=60)
        assert 'basis 411 length 4335' in standard.stdout
        lengths = [int(line.split()[-1]) for line in standard.stdout.splitlines() if line.startswith('modm ')]
        assert lengths == [1, 12, 66, 234, 625, 1365, 2568, 4335]
        trace_hash = hashlib.sha256((directory/'trace').read_bytes()).hexdigest()
    assert comb(18, 7) == 31824 and 31824-27489 == 4335
    assert 4335-sum(sextic['true_initial_hilbert_function']) == 1767
    paths = [Path(__file__), quadratic_path, normal_path, endpoint_path, sextic_path]
    paths += [ROOT/'code'/name for name in sources]
    paths += [ROOT/f'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_{name}_20261101.json' for name in ['CUBIC','QUARTIC','QUINTIC']]
    paths += [ROOT/'code'/name for name in ['s61_model.py', 's61_flag_transverse_obstruction_20261101.py', 'truncated_local_ideal_io_20261101.py']]
    return dict(id='S61_MINIMAL_FREE_FLAG_SEVENTH_UPPER_BOUND_20261101', all_passed=True,
                good_prime=prime, specialization=values, exact_modular_ideal_rank=27489,
                truncated_ring_columns=31824, modular_hilbert_function=[1, 11, 54, 168, 391, 740, 1203, 1767],
                generic_hilbert_function_through_six=sextic['true_initial_hilbert_function'],
                generic_seventh_hilbert_upper_bound=1767, generic_quotient_by_m8_length_upper_bound=4335,
                generic_new_seventh_initial_rank_lower_bound=30,
                known_seventh_quotient_before_new_relations=1797, known_seventh_m2_log=known_log,
                modular_rank_algorithms=['Complete sparse filtered ideal closure.', 'Truncated local standard basis with product/chain criteria.'],
                modular_input_sha256=input_hash, filtered_trace_sha256=trace_hash,
                field_argument='A good-prime minor of the universal finite-jet Macaulay matrix is nonzero over QQ(a0,a1,a2,a3,a4).',
                limitations=['This is a generic upper bound, not universal seventh-row equality or a proof of exact generic HF7.',
                             'No complete associated-graded ideal, homogenizer nilpotence, or flag-stratum saturation removal.'],
                proof_note='research/notes/S61_MINIMAL_FREE_FLAG_SEVENTH_UPPER_BOUND_20261101.md',
                input_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
                seconds=time.time()-started)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print('MINIMAL FREE FLAG SEVENTH UPPER BOUND PASS', result['seconds'])

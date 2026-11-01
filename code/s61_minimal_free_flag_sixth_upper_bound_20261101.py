#!/usr/bin/env python3
"""Generic sixth-jet upper bound from an exact good-prime rank certificate."""
import argparse
import hashlib
import json
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
    quintic_path = ROOT/'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUINTIC_20261101.json'
    quadratic = json.loads(quadratic_path.read_text())
    quintic = json.loads(quintic_path.read_text())
    assert quadratic['all_passed'] and quintic['all_passed']
    assert quintic['true_initial_hilbert_function'] == [1, 11, 54, 168, 391, 740]
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
    monomials = {degree: list(cwr(range(11), degree)) for degree in range(1, 7)}
    indices = {degree: {mon: i for i, mon in enumerate(mons)} for degree, mons in monomials.items()}
    W, jets = {1: tangent}, {}
    for degree in range(2, 7):
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
    with tempfile.TemporaryDirectory(prefix='minimal-sixth-upper-') as tmp:
        directory = Path(tmp)
        data = directory/'original-jets.txt'
        write_input(data, prime, 11, jets, maximum_degree=6)
        input_hash = hashlib.sha256(data.read_bytes()).hexdigest()
        binaries = []
        for source in sources:
            binary = directory/Path(source).stem
            subprocess.run(['c++', '-O3', '-std=c++17', str(ROOT/'code'/source), '-o', str(binary)], check=True, capture_output=True, timeout=60)
            binaries.append(binary)
        filtered = subprocess.run([str(binaries[0]), str(data), str(directory/'trace'), str(directory/'order'), str(directory/'initial')],
                                  text=True, capture_output=True, check=True, timeout=60)
        assert filtered.stdout.strip() == '9808'
        standard = subprocess.run([str(binaries[1]), str(data), str(directory/'standard')],
                                  text=True, capture_output=True, check=True, timeout=60)
        assert 'basis 314 length 2568' in standard.stdout
        lengths = [int(line.split()[-1]) for line in standard.stdout.splitlines() if line.startswith('modm ')]
        assert lengths == [1, 12, 66, 234, 625, 1365, 2568]
        trace_hash = hashlib.sha256((directory/'trace').read_bytes()).hexdigest()
    assert comb(17, 6) == 12376 and 12376-9808 == 2568
    assert 2568-sum(quintic['true_initial_hilbert_function']) == 1203
    paths = [Path(__file__), quadratic_path, normal_path, endpoint_path, quintic_path]
    paths += [ROOT/'code'/name for name in sources]
    paths += [ROOT/'code'/name for name in ['s61_model.py', 's61_flag_transverse_obstruction_20261101.py', 'truncated_local_ideal_io_20261101.py']]
    return dict(id='S61_MINIMAL_FREE_FLAG_SIXTH_UPPER_BOUND_20261101', all_passed=True,
                good_prime=prime, specialization=values, exact_modular_ideal_rank=9808,
                truncated_ring_columns=12376, modular_hilbert_function=[1, 11, 54, 168, 391, 740, 1203],
                generic_hilbert_function_through_five=quintic['true_initial_hilbert_function'],
                generic_sixth_hilbert_upper_bound=1203, generic_quotient_by_m7_length_upper_bound=2568,
                generic_new_sixth_initial_rank_lower_bound=34,
                modular_rank_algorithms=['Complete sparse filtered ideal closure.', 'Truncated local standard basis with product/chain criteria.'],
                modular_input_sha256=input_hash, filtered_trace_sha256=trace_hash,
                field_argument='A good-prime minor of the universal finite-jet Macaulay matrix is nonzero over QQ(a0,a1,a2,a3,a4).',
                limitations=['This is a generic upper bound, not universal sixth-row equality or a proof of exact generic HF6.',
                             'No complete associated-graded ideal, homogenizer nilpotence, or flag-stratum saturation removal.'],
                proof_note='research/notes/S61_MINIMAL_FREE_FLAG_SIXTH_UPPER_BOUND_20261101.md',
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
    print('MINIMAL FREE FLAG SIXTH UPPER BOUND PASS', result['seconds'])

#!/usr/bin/env python3
"""Exact true sixth initial relations over the five-parameter free-M0 flag chart."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import time
import sympy as s
from sympy.polys.rings import ring
from s61_minimal_free_flag_filtered_data_20261101 import ROOT
from s61_minimal_free_flag_sextic_backend_20261101 import Backend
from s61_minimal_free_flag_sextic_data_20261101 import seed_data
from s61_minimal_free_flag_sextic_filtered_20261101 import phase_a, phase_b, times_monomial


def core_inverse(matrix, parameters):
    """Invert using only ground pivots in QQ(a1)[a0,a4], then verify exactly."""
    K = s.QQ.frac_field(parameters[1])
    P, _, _ = ring('a0,a4', K)
    n = matrix.rows
    A = [[P.from_expr(value) for value in row]+[P.one if i == j else P.zero for j in range(n)]
         for i, row in enumerate(matrix.tolist())]
    permutation = list(range(n))
    for k in range(n):
        choices = [(i, j) for i in range(k, n) for j in range(k, n) if A[i][j] and A[i][j].is_ground]
        assert choices
        i, j = choices[0]
        A[k], A[i] = A[i], A[k]
        for row in A:
            row[k], row[j] = row[j], row[k]
        permutation[k], permutation[j] = permutation[j], permutation[k]
        pivot = A[k][k].LC
        A[k] = [value.quo_ground(pivot) for value in A[k]]
        for i in range(n):
            if i != k and A[i][k]:
                coefficient = A[i][k]
                A[i] = [x-coefficient*y for x, y in zip(A[i], A[k])]
    inverse = s.zeros(n)
    for i, j in enumerate(permutation):
        inverse[j, :] = s.Matrix([[A[i][n+k].as_expr() for k in range(n)]])
    assert (matrix.to_DM()*inverse.to_DM()).to_Matrix() == s.eye(n)
    return inverse


def determinant_mod(matrix, prime):
    A = [row.copy() for row in matrix]
    value = 1
    for k in range(len(A)):
        choices = [i for i in range(k, len(A)) if A[i][k]]
        assert choices
        i = choices[0]
        if i != k:
            A[k], A[i] = A[i], A[k]
            value = -value
        pivot = A[k][k]
        value = value*pivot % prime
        inverse = pow(pivot, -1, prime)
        for i in range(k+1, len(A)):
            coefficient = A[i][k]*inverse % prime
            for j in range(k+1, len(A)):
                A[i][j] = (A[i][j]-coefficient*A[k][j]) % prime
    return value % prime


def verify_sixth(data, basis, backend, progress):
    R, a, d, mm, ix, nf, G23, C34, C45, p, q, E5, E6 = [data[key] for key in
        ['R', 'a', 'd', 'mm', 'ix', 'nf', 'G23', 'C34', 'C45', 'p', 'q', 'E5', 'E6']]
    rp = lambda matrix: [[R.from_expr(value) for value in row] for row in matrix.tolist()]
    def dense(rows, columns, negative=False):
        result = []
        for row in rows:
            out = [R.zero]*columns
            for j, coefficient in row.items():
                out[j] = -coefficient if negative else coefficient
            result.append(out)
        return result
    product_tails = []
    for i, monomial in basis['metadata']:
        if i < 12:
            row = times_monomial(R, G23[i], 3, monomial, 6, mm, ix, nf, True)
            row = [value*d**(2+q) for value in row]
        elif i < 23:
            row = times_monomial(R, C34[i-12], 4, monomial, 6, mm, ix, nf, True)
            row = [value*d**q for value in row]
        else:
            row = times_monomial(R, C45[i-23], 5, monomial, 6, mm, ix, nf, True)
        product_tails.append(row)
    normalized_tails = backend.linear(R, [[R.zero]*1237 for _ in basis['combinations']],
                                      dense(basis['combinations'], 4323, True), product_tails)
    syzygy_tails = backend.linear(R, [[R.zero]*1237 for _ in basis['syzygies']],
                                  dense(basis['syzygies'], 4323, True), product_tails)
    progress('all degree-five product tails')
    reduced = backend.linear(R, E6, [[row[j] for j in basis['pivots']] for row in E5],
                             normalized_tails, d**(3+q))
    fifth_path = ROOT/'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUINTIC_20261101.json'
    fifth = json.loads(fifth_path.read_text())
    selected = fifth['selected_rows']
    M = s.Matrix(fifth['all_row_coefficients']).applyfunc(s.sympify)
    inverse = core_inverse(M[selected, :], a)
    for j, i in enumerate(selected):
        if i >= 267:
            inverse[:, j] = inverse[:, j]/a[1]**(4+p)
    inverse = inverse.applyfunc(s.cancel)
    denominator = s.lcm([s.denom(value) for value in inverse])
    denominator = s.Poly(denominator, a[1], domain=s.QQ).monic().as_expr()
    assert denominator.free_symbols <= {a[1]}
    numerator = inverse.applyfunc(lambda value: s.cancel(denominator*value))
    for i in range(267, 1180):
        M[i, :] = M[i, :]*a[1]**(4+p)
    # The final fifth-order normalization must use the same leading rows as
    # the frozen quintic certificate, including the common-denominator units.
    fifth_reduced = backend.linear(R, E5,
        [[row[j] for j in basis['pivots']] for row in E5], dense(basis['basis'], 3003))
    assert not any(row[j] for row in fifth_reduced for j in basis['pivots'])
    fifth_rref = [[(j, R.from_expr(s.Rational(c))) for j, c in enumerate(row) if c != '0']
                  for row in fifth['constant_fifth_rref']]
    fifth_coefficients = rp(M)
    for i, row in enumerate(fifth_reduced):
        remainder = [row[j] for j in basis['free']]
        for k, coefficient in enumerate(fifth_coefficients[i]):
            if coefficient:
                for j, value in fifth_rref[k]:
                    remainder[j] -= coefficient*value
        assert not any(remainder), ('fifth-leading compatibility', i)
    del fifth_reduced
    progress('full fifth leading compatibility with frozen Q5')
    coefficients = (M.to_DM()*numerator.to_DM()).to_Matrix()
    new = backend.linear(R, reduced, rp(coefficients), [reduced[i] for i in selected], R.from_expr(denominator))
    all_rows = new+syzygy_tails
    assert len(all_rows) == 3263
    progress('all 3263 genuine sixth relations')
    candidate_path = ROOT/'code/data/s61_minimal_free_flag_sextic_candidate_20261101.json'
    candidate = json.loads(candidate_path.read_text())
    pivots = candidate['pivots']
    canonical = [[R.zero]*1237 for _ in range(34)]
    for row, sparse in zip(canonical, candidate['rref_sparse']):
        for j, coefficient in sparse:
            row[j] = R.from_expr(s.Rational(coefficient))
    assert [[row[j] for j in pivots] for row in canonical] == [[R.one if i == j else R.zero for j in range(34)] for i in range(34)]
    sparse = [[(j, value) for j, value in enumerate(row) if value] for row in canonical]
    digest = hashlib.sha256()
    for i, source in enumerate(all_rows):
        remainder = source.copy()
        for k, pivot in enumerate(pivots):
            coefficient = source[pivot]
            if coefficient:
                for j, value in sparse[k]:
                    remainder[j] -= coefficient*value
                for exponent, value in sorted(coefficient.items()):
                    digest.update(f'{i} {k} {exponent} {value}\n'.encode())
        assert not any(remainder), ('candidate incomplete', i)
    prime, values = 32003, [1, 2, 3, 4, 5]
    def evaluate(polynomial):
        value = 0
        for exponent, coefficient in polynomial.items():
            term = int(coefficient.numerator)*pow(int(coefficient.denominator), -1, prime) % prime
            for x, power in zip(values, exponent):
                term = term*pow(x, power, prime) % prime
            value += term
        return value % prime
    minor = [[evaluate(all_rows[i][j]) for j in pivots] for i in candidate['selected_rows']]
    determinant = determinant_mod(minor, prime)
    assert determinant
    progress('universal membership, completeness, and independent 34-row minor')
    return dict(new_sixth_rank=34, known_sixth_rank=6771, quotient_dimension_before_new_rows=1237,
                relation_rows=3263, old_residual_rows=1180, independent_product_syzygies=2083,
                constant_sixth_rref_sparse=candidate['rref_sparse'], pivots=pivots,
                full_fifth_leading_compatibility_checked=True,
                selected_rows=candidate['selected_rows'], coefficient_matrix_sha256=digest.hexdigest(),
                sixth_monomials=[list(mon) for mon in mm[6]],
                selected_minor_prime=prime, selected_minor_point=values, selected_minor_determinant=determinant,
                normalization_core_denominator=str(s.factor(denominator)),
                quartic_row_denominator_power=p, quartic_generator_tail_denominator_power=q,
                nonzero_canonical_coefficients=sum(len(row) for row in candidate['rref_sparse']))


def run():
    started = time.time()
    stages = []
    def progress(message):
        elapsed = time.time()-started
        stages.append(dict(stage=message, seconds=elapsed))
        print(message, round(elapsed, 3), flush=True)
    with tempfile.TemporaryDirectory(prefix='minimal-universal-sixth-', dir='/tmp') as temp:
        backend = Backend(temp)
        seed = seed_data(backend, progress)
        basis, free6 = seed['product_basis'], seed['free6']
        previous = phase_a(seed, backend, progress)
        del seed
        data = phase_b(previous, backend, progress)
        del previous
        result = verify_sixth(data, basis, backend, progress)
        native_timings = backend.timings
    paths = [Path(__file__)]
    paths += [ROOT/'code'/name for name in [
        's61_minimal_free_flag_sextic_backend_20261101.py', 's61_minimal_free_flag_sextic_data_20261101.py',
        's61_minimal_free_flag_sextic_filtered_20261101.py', 'multivariate_jet_contraction_20261101.cpp',
        'multivariate_polynomial_matrix_20261101.cpp', 'data/s61_minimal_free_flag_sextic_candidate_20261101.json',
        's61_minimal_free_flag_filtered_data_20261101.py', 's61_minimal_free_flag_universal_quintic_20261101.py',
        's61_minimal_free_flag_universal_quartic_20261101.py', 's61_minimal_free_flag_universal_cubic_20261101.py',
        's61_model.py', 's61_flag_transverse_obstruction_20261101.py']]
    paths += [ROOT/f'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_{name}_20261101.json' for name in
              ['QUADRATIC', 'CUBIC', 'QUARTIC', 'QUINTIC']]
    paths += [ROOT/'reports'/name for name in ['S61_MINIMAL_FREE_FLAG_NORMAL_FORM_20261101.json',
              'S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json', 'S61_UNIVERSAL_WEIGHT6_20261031.json']]
    return dict(id='S61_MINIMAL_FREE_FLAG_UNIVERSAL_SEXTIC_20261101', all_passed=True,
                field='QQ(a0,a1,a2,a3,a4)', **result, quotient_free_columns=free6,
                true_initial_hilbert_function=[1, 11, 54, 168, 391, 740, 1203], quotient_by_m7_length=2568,
                stages_seconds=stages, native_timings=native_timings,
                scope=['Exact true initial ideal through degree six on the five-parameter function-field chart.',
                       'All 34 new sixth canonical rows have constant rational coefficients.',
                       'The prior fourth-order row still depends on a1.',
                       'No complete associated-graded equality, homogenizer nilpotence, or flag saturation removal.'],
                proof_note='research/notes/S61_MINIMAL_FREE_FLAG_UNIVERSAL_SEXTIC_20261101.md',
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
    print('UNIVERSAL SEXTIC PASS', result['seconds'])

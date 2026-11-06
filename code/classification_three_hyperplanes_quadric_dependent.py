#!/usr/bin/env python3
"""Full ambient ranks at one actual eight-parameter residual-source point."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import (
    COLUMNS, jacobi, sparse, vector, minor_certificate)

ROOT = Path(__file__).resolve().parents[1]


def projective_source():
    h, x, y, z, w = xx = s.symbols('h x y z w')
    a0, a1, *rs = params = s.symbols('a0 a1 r0:6')
    f = x*(x-h)
    q = (a0*h+a1*x)*y+f
    magnetic = sum(c*m for c, m in zip(rs, (h*h, h*x, h*y, x*x, x*y, y*y)))
    P = s.zeros(5)
    for ij, value in {(1, 3): f, (2, 3): y*(a1*y+2*x-h),
                      (2, 4): q, (3, 4): a1*y*w-magnetic}.items():
        P[ij] = s.expand(value)
        P[ij[::-1]] = -P[ij]
    return P, xx, tuple(params), f, q


def ambient_matrices(P, xx, directions):
    ps = sparse(P, xx)
    dcols = []
    for ij, exponents in COLUMNS:
        variation = {key: dict(value) for key, value in ps.items()}
        variation.setdefault(ij, {})[exponents] = variation.get(ij, {}).get(exponents, 0)+1
        dcols.append(jacobi(variation))
    D = s.Matrix.hstack(*dcols)
    X = s.Matrix(xx)
    vertical, orbit = [], []
    for a, b in itertools.product(range(5), repeat=2):
        field = s.zeros(5, 1)
        field[a] = xx[b]
        vertical.append(vector(sparse(X*field.T-field*X.T, xx)))
        E = s.zeros(5)
        E[a, b] = 1
        orbit.append(vector(sparse(E*P+P*E.T-xx[b]*P.diff(xx[a]), xx)))
    V = s.Matrix.hstack(*vertical)
    M = V.row_join(s.Matrix.hstack(*orbit)).row_join(vector(ps))
    F = M.row_join(s.Matrix.hstack(*(vector(sparse(R, xx)) for R in directions)))
    assert D*F == s.zeros(350, F.cols)
    return D, V, M, F


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU, (55, 56))
    resource.setrlimit(resource.RLIMIT_AS, (1536*1024**2, 1536*1024**2))
    begun = time.monotonic()
    source, xx, params, f, q = projective_source()
    universal_jacobi = [s.expand(sum(
        source[i, ell]*s.diff(source[j, k], xx[ell])
        + source[j, ell]*s.diff(source[k, i], xx[ell])
        + source[k, ell]*s.diff(source[i, j], xx[ell])
        for ell in range(5))) for i, j, k in itertools.combinations(range(5), 3)]
    assert universal_jacobi == [0]*10
    pf = s.expand(source[1, 2]*source[3, 4]-source[1, 3]*source[2, 4]
                  + source[1, 4]*source[2, 3])
    assert s.expand(pf+f*q) == 0
    values = dict(zip(params, (1, 1, 1, 2, 3, 5, 7, 11)))
    P = source.subs(values)
    ps = sparse(P, xx)
    assert all(not c for c in jacobi(ps, raw=True).values())
    assert jacobi(ps) == s.zeros(350, 1)
    assert (s.hessian(q.subs(values), xx)/2).rank() == 3
    directions = [source.diff(parameter).subs(values) for parameter in params]
    D, V, M, F = ambient_matrices(P, xx, directions)
    mats = {'full350_D': D, 'vertical': V, 'orbit_vertical_scale': M,
            'actual_full_source': F}
    ranks = {name: len(DomainMatrix.from_Matrix(A).convert_to(s.QQ).rref()[1])
             for name, A in mats.items()}
    certificates = {name: minor_certificate(A, ranks[name]) for name, A in mats.items()}
    lower = ranks['actual_full_source']-25
    upper = 125-ranks['full350_D']
    paths = ['code/classification_three_hyperplanes_quadric_dependent.py',
             'code/classification_rank3_mixed_component_review.py',
             'research/classification/THREE_HYPERPLANES_QUADRIC_DEPENDENT.md']
    report = dict(
        task_id='E_reduced_nonnormal', all_passed=True,
        source_parameters=list(map(str, params)),
        source_open='a0*(a0+a1) != 0; omitted factor choices retained by permutation of the three hyperplanes',
        universal_all_ten_raw_Jacobi_zero=True,
        universal_affine_pfaffian=str(-f*q),
        universal_actual_quintic=str(-xx[0]*f*q),
        all_eight_actual_source_directions_retained=True,
        point_substitution={str(k): str(v) for k, v in values.items()},
        point_brackets={str(ij): str(P[ij]) for ij in itertools.combinations(range(5), 2) if P[ij]},
        all8_actual_parameter_directions=[{str(ij): str(R[ij])
            for ij in itertools.combinations(range(5), 2) if R[ij]} for R in directions],
        full350_projective_Jacobi_zero=True,
        full350_times_all_actual_directions_zero=True,
        ranks_QQ=ranks, certificates=certificates,
        projective_family_dimension_lower=lower,
        projective_tangent_dimension_upper=upper,
        true_transverse_tangent_dimension=upper-lower,
        smooth_full_component_proved=(lower == upper),
        known_component_membership_not_checked=True,
        global_exhaustion_not_claimed=True,
        input_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        matrix_sha256={name: hashlib.sha256(json.dumps(A.tolist(), default=str).encode()).hexdigest()
                       for name, A in mats.items()},
        seconds=time.monotonic()-begun)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x') as handle:
        json.dump(report, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(ranks_QQ=ranks, lower=lower, upper=upper,
                         transverse=upper-lower, seconds=report['seconds'])))


def normal_quadrics(output):
    resource.setrlimit(resource.RLIMIT_CPU, (55, 56))
    resource.setrlimit(resource.RLIMIT_AS, (1536*1024**2, 1536*1024**2))
    begun = time.monotonic()
    rank_path = 'reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_20261107.json'
    snapshot_path = 'reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_SOURCE_20261107.json'
    saved = json.loads((ROOT/rank_path).read_text())
    snapshots = json.loads((ROOT/snapshot_path).read_text())
    for record in snapshots['sources']:
        assert hashlib.sha256(record['source'].encode()).hexdigest() == record['sha256']
        assert record['sha256'] == saved['input_sha256'][record['path']]
    xx = s.symbols('h x y z w')

    def decode(encoded):
        A = s.zeros(5)
        for ij in itertools.combinations(range(5), 2):
            A[ij] = s.sympify(encoded.get(str(ij), '0'))
            A[ij[::-1]] = -A[ij]
        return A

    P = decode(saved['point_brackets'])
    directions = [decode(data) for data in saved['all8_actual_parameter_directions']]
    assert jacobi(sparse(P, xx)) == s.zeros(350, 1)
    D, V, M, F = ambient_matrices(P, xx, directions)
    for name, matrix in [('full350_D', D), ('vertical', V),
                         ('orbit_vertical_scale', M), ('actual_full_source', F)]:
        assert hashlib.sha256(json.dumps(matrix.tolist(), default=str).encode()).hexdigest() == saved['matrix_sha256'][name]
    dm = DomainMatrix.from_Matrix(D).convert_to(s.QQ)
    reduced, pivots = dm.rref()
    reduced = reduced.to_Matrix()
    assert len(pivots) == 99
    kernel = []
    for j in range(150):
        if j in pivots:
            continue
        col = s.zeros(150, 1)
        col[j] = 1
        for i, pivot in enumerate(pivots):
            col[pivot] = -reduced[i, j]
        denominator = s.ilcm(1, *[s.denom(v) for v in col])
        kernel.append(col*denominator)
    K = s.Matrix.hstack(*kernel)
    assert D*K == s.zeros(350, 51)
    combined_pivots = list(DomainMatrix.from_Matrix(F.row_join(K)).convert_to(s.QQ).rref()[1])
    source_columns = [j for j in combined_pivots if j < F.cols]
    normal_indices = [j-F.cols for j in combined_pivots if j >= F.cols]
    assert len(source_columns) == 49 and len(normal_indices) == 2
    normals = []
    for index in normal_indices:
        A = s.zeros(5)
        for (ij, exponent), value in zip(COLUMNS, K[:, index]):
            A[ij] += value*s.prod(v**power for v, power in zip(xx, exponent))
        for ij in itertools.combinations(range(5), 2):
            A[ij[::-1]] = -A[ij]
        normals.append(A)
    squares = [jacobi(sparse(A, xx)) for A in normals]
    Q = s.Matrix.hstack(squares[0], jacobi(sparse(normals[0]+normals[1], xx))-sum(squares, s.zeros(350, 1)), squares[1])
    independent = D[:, list(pivots)]
    rows = list(DomainMatrix.from_Matrix(independent.T).convert_to(s.QQ).rref()[1])
    block = DomainMatrix.from_Matrix(D.extract(rows, list(pivots))).convert_to(s.QQ)
    inverse = block.inv()
    rhs = DomainMatrix.from_Matrix(Q.extract(rows, range(3))).convert_to(s.QQ)
    obstruction = Q-independent*inverse.matmul(rhs).to_Matrix()
    assert obstruction.extract(rows, range(3)) == s.zeros(99, 3)
    qr, qpivots = DomainMatrix.from_Matrix(obstruction).convert_to(s.QQ).rref()
    qr = qr.to_Matrix()
    n0, n1 = s.symbols('n0 n1')
    mons = [n0*n0, n0*n1, n1*n1]
    polynomials = [s.expand(sum(qr[i, j]*mons[j] for j in range(3))) for i in range(len(qpivots))]
    gb = s.groebner(polynomials, n0, n1, order='grevlex') if polynomials else None
    zero_dimensional = bool(gb and gb.is_zero_dimensional)
    selected_rows = list(DomainMatrix.from_Matrix(obstruction.T).convert_to(s.QQ).rref()[1])
    functionals = []
    for row in selected_rows:
        weights = DomainMatrix.from_Matrix(independent[row, :]).convert_to(s.QQ).matmul(inverse).to_Matrix()
        functional = s.zeros(1, 350)
        functional[row] = 1
        for j, selected in enumerate(rows):
            functional[selected] -= weights[j]
        assert functional*D == s.zeros(1, 150)
        coefficients = list(functional*Q)
        assert coefficients == list(obstruction[row, :])
        functionals.append(dict(row=row, nonzero_coefficients={str(j): str(value) for j, value in enumerate(functional) if value},
                                quadratic_coefficients=list(map(str, coefficients))))
    inputs = [rank_path, snapshot_path, 'code/classification_three_hyperplanes_quadric_dependent.py',
              'code/classification_rank3_mixed_component_review.py']
    report = dict(task_id='E_reduced_nonnormal', all_passed=True,
        frozen_all_four_matrix_hashes_reproduced=True, full_D_rank=99,
        actual_source_tangent_rank=49, full_kernel_dimension=51, normal_count=2,
        normal_kernel_indices=normal_indices, source_independent_column_indices=source_columns,
        normals=[{str(ij): str(A[ij]) for ij in itertools.combinations(range(5), 2) if A[ij]} for A in normals],
        normal_monomial_order=list(map(str, mons)), D_times_all_normals_zero=True,
        eliminated_D_rows=rows, eliminated_D_columns=list(pivots),
        invertible_D_minor_determinant=str(block.det()),
        quadratic_obstruction_rank=len(qpivots), quadratic_initial_generators=list(map(str, polynomials)),
        left_cokernel_functionals=functionals,
        grevlex_basis=[str(poly.as_expr()) for poly in gb.polys] if gb else [],
        quadratic_support_zero_dimensional=zero_dimensional,
        local_projective_dimension_upper_if_zero_dimensional=24 if zero_dimensional else None,
        source_componenthood_if_formal_slice_argument_accepted=zero_dimensional,
        local_uniqueness_and_full_initial_ideal_not_claimed=True,
        known_component_membership_not_checked=True,
        input_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
        producer_source=Path(__file__).read_text(),
        obstruction_matrix_sha256=hashlib.sha256(json.dumps(obstruction.tolist(), default=str).encode()).hexdigest(),
        seconds=time.monotonic()-begun)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x') as handle:
        json.dump(report, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(normal_count=2, quadratic_rank=len(qpivots), quadrics=list(map(str, polynomials)),
                         zero_dimensional=zero_dimensional, seconds=report['seconds'])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('ranks', 'normal'), default='ranks')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    (main if args.mode == 'ranks' else normal_quadrics)(args.output)

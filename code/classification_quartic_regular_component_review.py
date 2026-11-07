#!/usr/bin/env python3
"""Independent Fraction check of a nil4 quartic one-normal ambient slice.

Reuses only the older independent sparse polynomial and rational arithmetic.
No producer, legacy producer helper, SymPy or NumPy is imported.
"""
import argparse
import ast
import hashlib
import json
import math
import resource
import time
from itertools import product
from pathlib import Path

import classification_nil22_weighted_slice_review as q

ROOT = Path(__file__).resolve().parents[1]


def read_poly(text, parameters=None):
    parameters = parameters or {}
    coordinates = dict(zip(('h', 'x', 'y', 'z', 'w'), range(5)))
    def walk(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return q.atom(node.value)
        if isinstance(node, ast.Name):
            if node.id in coordinates:
                return q.atom(1, coordinates[node.id])
            return q.atom(parameters[node.id])
        if isinstance(node, ast.UnaryOp):
            return q.scale(walk(node.operand), -1 if isinstance(node.op, ast.USub) else 1)
        if isinstance(node, ast.BinOp):
            a = walk(node.left)
            if isinstance(node.op, ast.Pow):
                out = q.atom(1)
                for _ in range(node.right.value):
                    out = q.mul(out, a)
                return out
            b = walk(node.right)
            if isinstance(node.op, ast.Add):
                return q.add(a, b)
            if isinstance(node.op, ast.Sub):
                return q.add(a, q.scale(b, -1))
            if isinstance(node.op, ast.Mult):
                return q.mul(a, b)
            if isinstance(node.op, ast.Div):
                assert set(b) == {(0,)*q.N}
                return q.scale(a, 1/b[(0,)*q.N])
        raise ValueError(ast.dump(node))
    return walk(ast.parse(text, mode='eval').body)


def plus(A, B):
    return {ij: q.add(A.get(ij, {}), B.get(ij, {})) for ij in set(A)|set(B)}


def jvec(P):
    return q.jvec(q.jacobi(P))


def from_vector(vector):
    out = {}
    for value, (ij, exponent) in zip(vector, q.COLS):
        if value:
            out[ij] = q.add(out.get(ij, {}), {exponent+(0,): q.F(value)})
    return q.antisymmetric(out)


def decode_vector(data, length):
    return [q.F(data.get(str(i), '0')) for i in range(length)]


def encode_vector(vector):
    return {str(i): str(v) for i, v in enumerate(vector) if v}


def mixed(A, B, JA=None, JB=None):
    return [c-a-b for c, a, b in zip(jvec(plus(A, B)),
                                    jvec(A) if JA is None else JA,
                                    jvec(B) if JB is None else JB)]


def solve(A, b):
    """Independent square rational elimination, retaining an exact solution."""
    rows = [[q.F(v) for v in row]+[q.F(value)] for row, value in zip(A, b)]
    n = len(rows)
    assert all(len(row) == n+1 for row in rows)
    for j in range(n):
        k = next(i for i in range(j, n) if rows[i][j])
        rows[j], rows[k] = rows[k], rows[j]
        pivot = rows[j][j]
        rows[j] = [v/pivot for v in rows[j]]
        for i in range(j+1, n):
            value = rows[i][j]
            if value:
                rows[i] = [u-value*v for u, v in zip(rows[i], rows[j])]
    answer = [q.F(0)]*n
    for i in reversed(range(n)):
        answer[i] = rows[i][-1]-sum(rows[i][j]*answer[j] for j in range(i+1, n))
    assert [sum(a*v for a, v in zip(row, answer)) for row in A] == list(b)
    return answer


def certificate(A, data, compare_determinant=True):
    rows, columns = data['rows'], data['columns']
    value = q.determinant([[A[i][j] for j in columns] for i in rows])
    assert value
    if compare_determinant:
        assert value == q.F(data['determinant'])
    return dict(shape=[len(A), len(A[0])], rows=rows, columns=columns,
                determinant=str(value), size=len(columns))


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU, (165, 170))
    resource.setrlimit(resource.RLIMIT_AS, (1536*1024**2,)*2)
    started = time.monotonic()
    paths = [ROOT/'reports'/name for name in (
        'CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json',
        'CLASSIFICATION_HYPERPLANE_QUARTIC_NIL_REGULAR_COMPONENT_20261107.json')]
    source = json.loads(paths[0].read_text())['full_nilpotent_sources'][0]
    saved = json.loads(paths[1].read_text())['nil_regular_component']
    parameters = {f'n{i}': q.F(0) for i in range(36)}
    parameters.update(n1=q.F(1), n8=q.F(-1, 4), n12=q.F(-1, 3))
    P = {}
    for ij, text in zip(((1,2),(1,3),(1,4),(2,3),(2,4),(3,4)), source['brackets']):
        affine = read_poly(text, parameters)
        homogeneous = {}
        for e, c in affine.items():
            assert e[0] == e[5] == 0 and sum(e) <= 2
            homogeneous[(2-sum(e),)+e[1:]] = c
        P[ij] = homogeneous
    denominator = math.lcm(*(c.denominator for p in P.values() for c in p.values()))
    assert denominator == int(saved['integer_scale']) == 36
    P = q.antisymmetric({ij: q.scale(p, denominator) for ij, p in P.items()})
    for i, j in q.PAIRS:
        assert P.get((i,j), {}) == read_poly(saved['witness_brackets'].get(f'{i},{j}', '0'))
    assert all(not p for p in q.jacobi(P).values())
    raw_jacobi = {}
    for i, j, k in q.combinations(range(5), 3):
        raw_jacobi[i,j,k] = q.add(*(q.mul(P.get((a,l), {}), q.der(P.get((b,c), {}), l))
            for a,b,c in ((i,j,k),(j,k,i),(k,i,j)) for l in range(5)))
    assert all(not p for p in raw_jacobi.values())

    Dcols = []
    for ij, exponent in q.COLS:
        variation = {key: dict(value) for key, value in P.items()}
        variation[ij] = q.add(variation.get(ij, {}), {exponent+(0,): q.F(1)})
        variation[ij[::-1]] = q.scale(variation[ij], -1)
        Dcols.append(jvec(variation))
    D = q.transpose(Dcols)
    vertical, orbit = [], []
    for a, b in product(range(5), repeat=2):
        V, O = {}, {}
        for i, j in q.PAIRS:
            V[i,j] = q.add(q.atom(1,i,b) if j == a else {},
                           q.atom(-1,j,b) if i == a else {})
            O[i,j] = q.add(P.get((b,j), {}) if i == a else {},
                           P.get((i,b), {}) if j == a else {},
                           q.scale(q.mul(q.atom(1,b), q.der(P.get((i,j), {}), a)), -1))
        vertical.append(q.vec(V))
        orbit.append(q.vec(O))
    Mcols = vertical+orbit+[q.vec(P)]
    M = q.transpose(Mcols)
    Rvec = decode_vector(saved['independent_normal_sparse_vector'], 150)
    R = from_vector(Rvec)
    Kcols = Mcols+[Rvec]
    K = q.transpose(Kcols)
    assert all(sum(a*b for a,b in zip(row,col)) == 0 for row in D for col in Kcols)
    assert q.rank(q.transpose(vertical))[0] == 24
    assert q.rank(D)[0] == 102
    assert q.rank(M)[0] == 47 and q.rank(K)[0] == 48
    certificates = dict(D=certificate(D, saved['derivative']),
                        actual_orbit_vertical_scale=certificate(M, saved['orbit_vertical_scale'], False),
                        actual_orbit_plus_normal=certificate(K, saved['whole_kernel'], False))
    Q = jvec(R)
    assert Q == decode_vector(saved['quadratic_coefficient_sparse'], 350)
    assert q.rank([row+[v] for row,v in zip(D,Q)])[0] == 102
    rows, columns = saved['derivative']['rows'], saved['derivative']['columns']
    A = [[D[i][j] for j in columns] for i in rows]
    independent_solution = solve(A, [-Q[i] for i in rows])
    Tvec = [q.F(0)]*150
    for j, v in zip(columns, independent_solution):
        Tvec[j] = v
    assert Tvec == decode_vector(saved['quadratic_correction_sparse'], 150)
    assert all(sum(a*b for a,b in zip(row,Tvec))+v == 0 for row,v in zip(D,Q))
    T = from_vector(Tvec)
    cubic = mixed(R, T, Q, jvec(T))
    assert cubic == decode_vector(saved['cubic_coefficient_sparse'], 350)
    # M+R is proved to span the ENTIRE 48-dimensional kernel. Therefore this
    # verifies correction ambiguity for every kernel vector, not one slice.
    mixed_columns = [mixed(R, from_vector(v), Q) for v in Kcols]
    augmented = [row+list(values) for row,values in zip(D,q.transpose(mixed_columns))]
    assert q.rank(augmented)[0] == 102
    assert q.rank([row+[v] for row,v in zip(augmented,cubic)])[0] == 103
    # A literal left-cokernel witness, supported on 103 original equation rows.
    # Express the additional row through the independent 102x102 row minor.
    extra = next(i for i in saved['cubic_obstruction']['rows'] if i not in rows)
    coefficients = solve(q.transpose(A), [-D[extra][j] for j in columns])
    left = [q.F(0)]*350
    for i,v in zip(rows,coefficients):
        left[i] = v
    left[extra] = 1
    assert all(sum(a*b for a,b in zip(left,col)) == 0 for col in Dcols)
    assert all(sum(a*b for a,b in zip(left,col)) == 0 for col in mixed_columns)
    pairing = sum(a*b for a,b in zip(left,cubic))
    assert pairing
    cubic_minor_rows = rows+[extra]
    cubic_minor_cols = columns+[150]
    cubic_minor = [[D[i][j] if j < 150 else cubic[i] for j in cubic_minor_cols]
                   for i in cubic_minor_rows]
    value = q.determinant(cubic_minor)
    assert value and value == q.F(certificates['D']['determinant'])*pairing
    certificates['cubic'] = dict(rows=cubic_minor_rows, columns=cubic_minor_cols,
                                  determinant=str(value), size=103)
    # C150=M47 direct-sum (U102 direct-sum R). D injects on U, while
    # M+R spans kerD. Thus this is a genuine 47-cut affine slice, not a quotient.
    report = dict(task_id='E_reduced_nonnormal', all_passed=True,
        no_producer_import=True, no_CAS_import=True,
        method='Independent standard-library Fraction full150 sparse coefficient reconstruction; only the old independent nil22 arithmetic helper is reused.',
        exact_field='QQ', integer_scale=denominator,
        all_six_witness_brackets_reconstructed_from_full_source=True,
        all_ten_raw_Jacobi_zero=True, all350_projective_equations_zero=True,
        derivative_shape=[350,150], derivative_rank=102, vertical_rank=24,
        actual_orbit_vertical_scale_rank=47, entire_kernel_dimension=48,
        actual_plus_R_spans_entire_kernel=True,
        all_D_times_actual_and_normal_columns_zero=True,
        certificates=certificates,
        quadratic_vector_independently_reconstructed=encode_vector(Q),
        exact_pivot_correction_independently_solved=encode_vector(Tvec),
        full_D_T_plus_quadratic_zero=True,
        entire_kernel_mixed_correction_rank=102,
        cubic_augmented_rank=103,
        left_functional=encode_vector(left),
        left_annihilates_full_D_and_every_kernel_mixed_column=True,
        nonzero_cubic_pairing=str(pairing),
        true_slice_dimension=103, true_slice_linear_cut_count=47,
        literal_slice='U is the102 standard coordinate columns indexed by the verified D minor; S=U+R complements the47-dimensional actual orbit/vertical/scale space.',
        formal_support='The102 selected rows eliminate U by IFT. The actual quadratic correction solves every350 row. The displayed left functional gives a nonzero cubic coefficient in the remaining one variable. Every kernel ambiguity changes the cubic by an element of imD, so the true slice is supported only at the origin.',
        projective_orbit_dimension=22, projective_ambient_local_dimension=22,
        orbit_closure_is_full_component=True,
        reduced_ambient_component_at_witness_unique=True,
        uniqueness_reason='Every global component through the witness is stable under the connected PGL5 action, so it contains the same22-dimensional orbit closure; the verified local dimension bound22 and maximality force equality.',
        new_global_component_count_claimed=False,
        distinctness_from_previous89_checked=False,
        exact_nonreduced_structure_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in paths+[Path(q.__file__)]},
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        reviewer_source_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-started,
        maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with output.open('x') as handle:
        json.dump(report, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k:report[k] for k in ('all_passed','seconds','maximum_RSS_KiB','nonzero_cubic_pairing')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)

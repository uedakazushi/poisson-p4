#!/usr/bin/env python3
"""Exact true cubic initial relations on the five-parameter free-M0 flag chart."""
import argparse
import hashlib
import json
import time
from itertools import combinations_with_replacement
from pathlib import Path

import sympy as s
from sympy.polys.rings import ring
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize

ROOT = Path(__file__).resolve().parents[1]
QUADRATIC = 'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUADRATIC_20261101.json'
NORMAL = 'reports/S61_MINIMAL_FREE_FLAG_NORMAL_FORM_20261101.json'
POINT = 'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json'
enc = lambda matrix: [[str(v) for v in row] for row in matrix.tolist()]


def build():
    """Return polynomial numerator data; actual cubic residuals are P/a1**2."""
    report = json.loads((ROOT / QUADRATIC).read_text())
    a = s.symbols('a0:5')
    R, *_ = ring(','.join(map(str, a)), s.QQ)
    parse = lambda name: s.Matrix(report[name]).applyfunc(s.sympify)
    J, T, U = map(parse, ['jacobian', 'tangent_basis_after_shift', 'right_inverse'])
    rows, selected = report['pivot_rows'], report['selected_quadratic_rows']
    UP = U.applyfunc(lambda value: s.cancel(a[1] * value))
    JUP = (J.to_DM() * UP.to_DM()).to_Matrix()
    M = parse('quadratic_row_coefficients')
    selected_inverse = M[selected, :].inv(method='DM')
    N = (M.to_DM() * selected_inverse.to_DM()).to_Matrix().applyfunc(s.cancel)
    denominator = s.factor(s.lcm([s.denom(value) for value in N]))
    assert denominator.is_Rational
    NP = N.applyfunc(lambda value: s.cancel(value * denominator))
    model = S61()
    one = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 1]
    beta = s.Matrix(json.loads((ROOT / NORMAL).read_text())['normalized_beta16']).applyfunc(s.sympify)
    two, F = specialize(model, dict(zip(one, beta)), True)
    point = json.loads((ROOT / POINT).read_text())['strata']['free_M0']['local_data']
    fixed, chart = two.index(point['chart_fixed_u'] - 1), point['chart_columns']
    mon2 = list(combinations_with_replacement(range(11), 2))
    mon3 = list(combinations_with_replacement(range(11), 3))
    index3 = {mon: i for i, mon in enumerate(mon3)}
    rp = lambda matrix: [[R.from_expr(v) for v in row] for row in matrix.tolist()]

    def multiply(left, right):
        out = [[R.zero for _ in right[0]] for _ in left]
        for i, row in enumerate(left):
            for k, coefficient in enumerate(row):
                if coefficient:
                    for j, value in enumerate(right[k]):
                        if value:
                            out[i][j] += coefficient * value
        return out

    TT = rp(T)
    H2 = [[R.zero for _ in mon2] for _ in F]
    terms = []
    for row, polynomial in enumerate(F):
        terms.append([])
        for mon, coefficient in polynomial.items():
            if fixed in mon:
                continue
            i, j = [chart.index(k) for k in mon] + [44] * (2 - len(mon))
            coefficient = R.from_expr(coefficient)
            terms[-1].append((i, j, coefficient))
            for k, (v, w) in enumerate(mon2):
                H2[row][k] += coefficient * (TT[i][v] * TT[j][w] +
                    (TT[i][w] * TT[j][v] if v != w else R.zero))
    W2num = multiply(rp(-UP), [H2[i] for i in rows])
    H3num = [[R.zero for _ in mon3] for _ in F]
    for row, row_terms in enumerate(terms):
        for i, j, coefficient in row_terms:
            for left, right in [(i, j), (j, i)]:
                for k, value in enumerate(TT[left]):
                    if value:
                        for ell, correction in enumerate(W2num[right]):
                            if correction:
                                H3num[row][index3[tuple(sorted((k,) + mon2[ell]))]] += coefficient * value * correction
    C2 = parse('constant_quadratic_rref')
    products = s.zeros(132, len(mon3))
    for row in range(12):
        for k in range(11):
            for j, mon in enumerate(mon2):
                products[11 * row + k, index3[tuple(sorted(mon + (k,)))]] += C2[row, j]
    quotient, pivots = products.to_DM().rref()
    quotient = quotient.to_Matrix()[:len(pivots), :]
    assert len(pivots) == 107
    free = [i for i in range(len(mon3)) if i not in pivots]
    qrows = rp(quotient)
    projected = []
    for row in H3num:
        row = row.copy()
        for pivot, qrow in zip(pivots, qrows):
            coefficient = row[pivot]
            if coefficient:
                for j, value in enumerate(qrow):
                    if value:
                        row[j] -= coefficient * value
        projected.append([row[j] for j in free])
    correction = multiply(rp(JUP), [projected[i] for i in rows])
    P = [[R.from_expr(a[1]) * v - w for v, w in zip(row, corr)]
         for row, corr in zip(projected, correction)]
    correction = multiply(rp(NP), [P[i] for i in selected])
    C = [[R.from_expr(denominator) * v - w for v, w in zip(row, corr)]
         for row, corr in zip(P, correction)]
    CC = s.Matrix([[v.as_expr() for v in row] for row in C])
    return dict(a=a, J=J, T=T, U=U, UP=UP, M=M, rows=rows, selected=selected,
                selected_inverse=selected_inverse, C2=C2, mon2=mon2, mon3=mon3,
                quotient=quotient, quotient_pivots=list(pivots), free=free,
                W2num=W2num, projected=projected, P=P, CC=CC,
                row_cancellation_denominator=denominator, multiply=multiply, rp=rp,
                terms=terms, R=R, H2=H2, H3num=H3num)


def run():
    started = time.time()
    data = build()
    a, CC, free, mon3 = [data[k] for k in ['a', 'CC', 'free', 'mon3']]
    sample = dict(zip(a, range(1, 6)))
    selected = list(CC.subs(sample).T.to_DM().rref()[1])
    assert len(selected) == 11
    rr, pivots = CC[selected, :].to_DM().rref()
    rr = rr.to_Matrix().applyfunc(s.factor)
    assert len(pivots) == 11
    assert (CC.to_DM() - CC[:, list(pivots)].to_DM() * rr.to_DM()).is_zero_matrix
    v = s.symbols('v0:11')
    cubic_polynomials = [sum(row[j] * s.prod(v[k] for k in mon3[free[j]])
                             for j in range(len(free))) for row in rr.tolist()]
    constant_cubics = [s.expand(poly.subs(v[0], v[0] - a[1] * v[1] / 6)) for poly in cubic_polynomials]
    assert all(not poly.free_symbols.intersection(a) for poly in constant_cubics)
    C2 = data['C2']
    assert all(C2[:, j] == s.zeros(12, 1) for j, mon in enumerate(data['mon2']) if 0 in mon or 1 in mon)
    cubic_matrix = s.zeros(11, len(mon3))
    for i, poly in enumerate(constant_cubics):
        coefficients = s.Poly(poly, *v)
        for j, mon in enumerate(mon3):
            cubic_matrix[i, j] = coefficients.coeff_monomial(s.prod(v[k] for k in mon))
    assert data['quotient'].col_join(cubic_matrix).to_DM().rank() == 118
    # Correct t's cubic term by the FULL residuals whose quadratic terms remove t2.
    p2 = list(C2.to_DM().rref()[1])
    t2 = s.Matrix([[x.as_expr() / a[1] for x in data['W2num'][44]]])
    assert (t2.to_DM() - t2[:, p2].to_DM() * C2.to_DM()).is_zero_matrix
    lam = (t2[:, p2].to_DM() * data['selected_inverse'].to_DM()).to_Matrix()
    w3t = data['multiply']([data['rp'](-data['UP'])[44]],
                          [data['projected'][i] for i in data['rows']])[0]
    t3 = s.Matrix([[x.as_expr() / a[1]**2 for x in w3t]]) - lam * s.Matrix(
        [[x.as_expr() / a[1]**2 for x in data['P'][i]] for i in data['selected']])
    remainder = (t3.to_DM() - t3[:, list(pivots)].to_DM() * rr.to_DM()).to_Matrix().applyfunc(s.factor)
    assert any(remainder) and not any(x.free_symbols for x in remainder)
    first = next(i for i, value in enumerate(remainder) if value)
    assert remainder[first] == -s.Rational(25, 81)
    sources = [Path(__file__), ROOT / QUADRATIC, ROOT / NORMAL, ROOT / POINT,
               ROOT / 'code/s61_model.py', ROOT / 'code/s61_flag_transverse_obstruction_20261101.py',
               ROOT / 'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_MINIMAL_FREE_FLAG_UNIVERSAL_CUBIC_20261101', all_passed=True,
                field='QQ(a0,a1,a2,a3,a4)', quadratic_product_rank=107,
                new_cubic_rank=11, true_initial_hilbert_function=[1, 11, 54, 168],
                cubic_monomials=[list(mon) for mon in mon3],
                quadratic_product_pivots=data['quotient_pivots'], quotient_free_columns=free,
                cubic_selected_rows=selected, cubic_pivots=list(pivots),
                cubic_row_coefficients=enc(CC[:, list(pivots)]), cubic_rref=enc(rr),
                row_cancellation_denominator=str(data['row_cancellation_denominator']),
                second_tangent_change='v0_old = v0_new - a1*v1_new/6',
                constant_cubic_polynomials=list(map(str, constant_cubics)),
                t_order=3, t2_zero_mod_quadrics=True,
                t3_normal_form=list(map(str, remainder)),
                t3_first_nonzero_monomial=list(mon3[free[first]]),
                t3_first_nonzero_coefficient=str(remainder[first]),
                failed_distinct_probe=dict(field='GF(32003)', variables=45,
                    order='Weights=>toList(45:-1), Global=>false', command='M2-binary --script INPUT',
                    input_sha256='883a7f86f620286f0956956d3b6b7c30d81e9955b6e8a0fb4302205b78cd09d2',
                    timeout_seconds=180, result='INPUT_READY only; no completed gb or dimension result'),
                scope=['Exact function-field TRUE initial ideal through degree3; both tangent changes make it constant.',
                       'The generic completed ring has t in m^3 but not m^4.',
                       'No higher initial-ideal equality, ambient nilpotence, flatness or saturation-removal conclusion.'],
                proof_note='research/notes/S61_MINIMAL_FREE_FLAG_UNIVERSAL_CUBIC_20261101.md',
                input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                seconds=time.time()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(all_passed=True, seconds=result['seconds'])))

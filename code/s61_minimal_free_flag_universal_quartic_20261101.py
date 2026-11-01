#!/usr/bin/env python3
"""True quartic initial ideal over the five-parameter free-M0 flag chart."""
import argparse
import hashlib
import json
import time
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_minimal_free_flag_universal_cubic_20261101 import build as cubic_data

ROOT = Path(__file__).resolve().parents[1]
CUBIC = 'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_CUBIC_20261101.json'
enc = lambda matrix: [[str(v) for v in row] for row in matrix.tolist()]


def quartic_quotient(C2, mon2, mon4):
    report = json.loads((ROOT / CUBIC).read_text())
    v = s.symbols('v0:11')
    index = {mon: j for j, mon in enumerate(mon4)}
    products = s.zeros(12 * 66 + 11 * 11, len(mon4))
    for i in range(12):
        for j, left in enumerate(mon2):
            for k, right in enumerate(mon2):
                products[i * 66 + j, index[tuple(sorted(left + right))]] += C2[i, k]
    for i, expression in enumerate(report['constant_cubic_polynomials']):
        polynomial = s.Poly(s.sympify(expression), *v)
        for exponents, coefficient in polynomial.terms():
            mon = tuple(j for j, exponent in enumerate(exponents) for _ in range(exponent))
            for j in range(11):
                products[792 + 11 * i + j, index[tuple(sorted(mon + (j,)))]] += coefficient
    rr, pivots = products.to_DM().rref()
    assert len(pivots) == 595
    return rr.to_Matrix()[:len(pivots), :], list(pivots)


def run():
    started = time.time()
    stages = {}
    data = cubic_data()
    R, a, rp, multiply = [data[k] for k in ['R', 'a', 'rp', 'multiply']]
    d = R.from_expr(a[1])
    mon2, mon3 = data['mon2'], data['mon3']
    mon4 = list(combinations_with_replacement(range(11), 4))
    index4 = {mon: j for j, mon in enumerate(mon4)}

    def shift(matrix, monomials):
        """Exact v0_old=v0_new-a1*v1_new/6, also applied to implicit jets."""
        indices = {mon: j for j, mon in enumerate(monomials)}
        out = [[R.zero for _ in monomials] for _ in matrix]
        for j, mon in enumerate(monomials):
            power = mon.count(0)
            for k in range(power + 1):
                target = list(mon)
                for _ in range(k):
                    target.remove(0)
                    target.append(1)
                column = indices[tuple(sorted(target))]
                coefficient = R.from_expr(s.binomial(power, k) * (-a[1] / 6)**k)
                for row, source in enumerate(matrix):
                    if source[j]:
                        out[row][column] += coefficient * source[j]
        return out

    T = data['T'].copy()
    T[:, 1] -= a[1] * T[:, 0] / 6
    TT = rp(T)
    W2 = shift(data['W2num'], mon2)
    H3 = shift(data['H3num'], mon3)
    W3 = multiply(rp(-data['UP']), [H3[i] for i in data['rows']])
    stages['cubic_jets'] = time.time() - started

    def convolution(left, right, left_mons, right_mons):
        out = [R.zero] * len(mon4)
        for i, x in enumerate(left):
            if x:
                for j, y in enumerate(right):
                    if y:
                        out[index4[tuple(sorted(left_mons[i] + right_mons[j]))]] += x * y
        return out

    pairs = sorted({(min(i, j), max(i, j)) for row in data['terms'] for i, j, _ in row})
    cache = {}
    mon1 = [(i,) for i in range(11)]
    for i, j in pairs:
        first = convolution(TT[i], W3[j], mon1, mon3)
        second = convolution(TT[j], W3[i], mon1, mon3)
        third = convolution(W2[i], W2[j], mon2, mon2)
        cache[i, j] = [x + y + z for x, y, z in zip(first, second, third)]
    H4 = [[R.zero] * len(mon4) for _ in range(135)]
    for row, terms in enumerate(data['terms']):
        for i, j, coefficient in terms:
            for k, value in enumerate(cache[min(i, j), max(i, j)]):
                if value:
                    H4[row][k] += coefficient * value
    stages['raw_quartic'] = time.time() - started
    raw_term_count = sum(len(value) for row in H4 for value in row)
    quotient, piv4 = quartic_quotient(data['C2'], mon2, mon4)
    free4 = [j for j in range(len(mon4)) if j not in piv4]
    sparse_quotient = [[(j, value) for j, value in enumerate(row) if value] for row in rp(quotient)]

    def project(matrix):
        out = []
        for row in matrix:
            row = row.copy()
            for pivot, relation in zip(piv4, sparse_quotient):
                coefficient = row[pivot]
                if coefficient:
                    for j, value in relation:
                        row[j] -= coefficient * value
            out.append([row[j] for j in free4])
        return out

    H4p = project(H4)
    JUP = rp((data['J'].to_DM() * data['UP'].to_DM()).to_Matrix())
    correction = multiply(JUP, [H4p[i] for i in data['rows']])
    A4 = [[d*x-y for x, y in zip(row, corr)] for row, corr in zip(H4p, correction)]
    correction = multiply(JUP, [H3[i] for i in data['rows']])
    A3 = [[d*x-y for x, y in zip(row, corr)] for row, corr in zip(H3, correction)]
    selected = data['selected']
    G23 = multiply(rp(data['selected_inverse']), [A3[i] for i in selected])
    N = rp((data['M'].to_DM() * data['selected_inverse'].to_DM()).to_Matrix())
    correction = multiply(N, [A3[i] for i in selected])
    G33 = [[x-y for x, y in zip(row, corr)] for row, corr in zip(A3, correction)]
    correction = multiply(N, [A4[i] for i in selected])
    G34 = [[x-y for x, y in zip(row, corr)] for row, corr in zip(A4, correction)]
    stages['projected_jets'] = time.time() - started
    # Keep the actual fourth-order tails of the 132 linear multiples of G2.
    index3 = {mon: j for j, mon in enumerate(mon3)}
    P3 = s.zeros(132, len(mon3))
    P4 = [[R.zero] * len(mon4) for _ in range(132)]
    for i in range(12):
        for j in range(11):
            for k, mon in enumerate(mon2):
                P3[11*i+j, index3[tuple(sorted(mon+(j,)))]] += data['C2'][i, k]
            for k, mon in enumerate(mon3):
                P4[11*i+j][index4[tuple(sorted(mon+(j,)))]] += G23[i][k]
    P4 = project(P4)
    product_rows = list(P3.T.to_DM().rref()[1])
    product_pivots = data['quotient_pivots']
    assert len(product_rows) == 107
    product_inverse = P3[product_rows, product_pivots].inv(method='DM')
    B4 = multiply(rp(product_inverse), [P4[i] for i in product_rows])
    correction = multiply(rp(P3[:, product_pivots]), B4)
    product_syzygies = [[x-y for x, y in zip(row, corr)] for row, corr in zip(P4, correction)]
    # Eliminate cubic quadratic-products, then the eleven true cubic generators.
    G3pivots = [[row[j] for j in product_pivots] for row in G33]
    correction = multiply(G3pivots, rp(data['quotient']))
    G33reduced = [[x-y for x, y in zip(row, corr)] for row, corr in zip(G33, correction)]
    C3 = s.Matrix([[value.exquo(d*d).as_expr() for value in row] for row in G33reduced])
    cubic_rows = list(C3.T.to_DM().rref()[1])
    rr, cubic_pivots = C3[cubic_rows, :].to_DM().rref()
    assert len(cubic_pivots) == 11 and not any(v.free_symbols for v in rr.to_Matrix())
    cubic_inverse = C3[cubic_rows, list(cubic_pivots)].inv(method='DM')
    assert not any(s.denom(s.cancel(v)).free_symbols for v in cubic_inverse)
    correction = multiply(G3pivots, B4)
    G34reduced = [[d*x-y for x, y in zip(row, corr)] for row, corr in zip(G34, correction)]
    C4 = multiply(rp(cubic_inverse), [G34reduced[i] for i in cubic_rows])
    correction = multiply(rp(C3[:, list(cubic_pivots)]), C4)
    new4 = [[x-y for x, y in zip(row, corr)] for row, corr in zip(G34reduced, correction)] + product_syzygies
    matrix = s.Matrix([[value.as_expr() for value in row] for row in new4])
    sample = dict(zip(a, range(1, 6)))
    selected4 = list(matrix.subs(sample).T.to_DM().rref()[1])
    assert len(selected4) == 15
    rr, pivots = matrix[selected4, :].to_DM().rref()
    rr = rr.to_Matrix().applyfunc(s.factor)
    assert len(pivots) == 15
    assert (matrix.to_DM() - matrix[:, list(pivots)].to_DM()*rr.to_DM()).is_zero_matrix
    parameter_entries = [(i, j) for i in range(rr.rows) for j in range(rr.cols) if rr[i, j].free_symbols]
    assert len(parameter_entries) == 1
    i, j = parameter_entries[0]
    assert rr[i, j].free_symbols == {a[1]} and s.degree(rr[i, j], a[1]) == 1
    assert mon4[free4[j]] == (3, 3, 3, 10)
    stages['all_quartic_relations'] = time.time() - started
    paths = [Path(__file__)] + [ROOT / name for name in [
        'code/s61_minimal_free_flag_universal_cubic_20261101.py', CUBIC,
        'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUADRATIC_20261101.json',
        'reports/S61_MINIMAL_FREE_FLAG_NORMAL_FORM_20261101.json',
        'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json',
        'reports/S61_UNIVERSAL_WEIGHT6_20261031.json', 'code/s61_model.py',
        'code/s61_flag_transverse_obstruction_20261101.py']]
    return dict(id='S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUARTIC_20261101', all_passed=True,
                field='QQ(a0,a1,a2,a3,a4)', lower_quartic_rank=595, new_quartic_rank=15,
                true_initial_hilbert_function=[1, 11, 54, 168, 391],
                monomials=[list(mon) for mon in mon4], quotient_pivots=piv4, quotient_free_columns=free4,
                quadratic_product_rows=product_rows, quadratic_product_pivots=product_pivots,
                cubic_rows=cubic_rows, cubic_pivots=list(cubic_pivots),
                new_quartic_selected_rows=selected4, new_quartic_pivots=list(pivots),
                new_quartic_coefficients=enc(matrix[:, list(pivots)]), new_quartic_rref=enc(rr),
                parameter_entry=dict(row=i, column=j, monomial=list(mon4[free4[j]]), coefficient=str(rr[i, j])),
                nuisance_parameters_absent=['a0', 'a2', 'a3', 'a4'],
                raw_quartic_parameter_term_count=raw_term_count, stages_seconds=stages,
                scope=['TRUE initial ideal through degree4 over the generic five-variable chart.',
                       'Four parameters disappear completely from the quartic rowspace; a1 occurs linearly in one coefficient.',
                       'No full associated-graded equality, ambient t-nilpotence, flatness or saturation-removal assertion.'],
                proof_note='research/notes/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUARTIC_20261101.md',
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
    print(json.dumps(dict(all_passed=True, seconds=result['seconds'], stages=result['stages_seconds'])))

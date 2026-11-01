#!/usr/bin/env python3
"""True fifth-order initial relations over the five-parameter free-M0 chart."""
import argparse
import hashlib
import json
import time
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_minimal_free_flag_filtered_data_20261101 import build_four, ROOT

enc = lambda matrix: [[str(v) for v in row] for row in matrix.tolist()]


def fifth_quotient(state, monomials):
    """Constant rank2234, followed by six rows over QQ[a1]."""
    R, rp = state['R'], state['rp']
    C2 = state['data']['C2']
    v = s.symbols('v0:11')
    index = {mon: j for j, mon in enumerate(monomials[5])}
    generators = [(2, {mon: C2[i, j] for j, mon in enumerate(monomials[2]) if C2[i, j]}) for i in range(12)]
    cubic = json.loads((ROOT / 'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_CUBIC_20261101.json').read_text())
    for expression in cubic['constant_cubic_polynomials']:
        polynomial = s.Poly(s.sympify(expression), *v)
        generators.append((3, {tuple(j for j, exponent in enumerate(ex) for _ in range(exponent)): co for ex, co in polynomial.terms()}))
    quartics = [{monomials[4][column]: state['rr'][i, j] for j, column in enumerate(state['free4']) if state['rr'][i, j]} for i in range(15)]
    generators += [(4, row) for row in quartics[1:]]
    entries, count = {}, 0
    for degree, polynomial in generators:
        for mon in monomials[5-degree]:
            for exponent, coefficient in polynomial.items():
                key = count, index[tuple(sorted(mon+exponent))]
                entries[key] = entries.get(key, 0) + coefficient
            count += 1
    matrix = s.SparseMatrix(count, 3003, entries)
    rr, pivots = matrix.to_DM().rref()
    assert len(pivots) == 2234
    free = [j for j in range(3003) if j not in pivots]
    domain = rr.domain
    rows = [[(j, R.from_expr(domain.to_sympy(value))) for j, value in row.items()] for _, row in sorted(rr.to_sdm().items())]

    def constant_project(row):
        row = row.copy()
        for pivot, relation in zip(pivots, rows):
            coefficient = row[pivot]
            if coefficient:
                for j, value in relation:
                    row[j] -= coefficient * value
        return [row[j] for j in free]

    variable_rows = []
    for j in range(11):
        row = [R.zero]*3003
        for mon, coefficient in quartics[0].items():
            row[index[tuple(sorted(mon+(j,)))]] += R.from_expr(coefficient)
        variable_rows.append(constant_project(row))
    variable = s.Matrix([[value.as_expr() for value in row] for row in variable_rows])
    vr, vp = variable.to_DM().rref()
    assert len(vp) == 6
    vr = vr.to_Matrix()[:len(vp), :].applyfunc(s.factor)
    assert not any(s.denom(value).free_symbols for value in vr)
    variable_basis = [[(j, value) for j, value in enumerate(row) if value] for row in rp(vr)]
    variable_free = [j for j in range(len(free)) if j not in vp]

    def project(row):
        row = constant_project(row)
        for pivot, relation in zip(vp, variable_basis):
            coefficient = row[pivot]
            if coefficient:
                for j, value in relation:
                    row[j] -= coefficient * value
        return [row[j] for j in variable_free]

    return project, [free[j] for j in variable_free], dict(constant_rank=2234, variable_rank=6, variable_pivots=list(vp), variable_rref=enc(vr))


def reduce_rows(matrix, basis, pivots, free, R):
    sparse = [[(j, value) for j, value in enumerate(row) if value] for row in basis]
    out = []
    for row in matrix:
        row = row.copy()
        for pivot, relation in zip(pivots, sparse):
            coefficient = row[pivot]
            if coefficient:
                for j, value in relation:
                    row[j] -= coefficient * value
        out.append([row[j] for j in free])
    return out


def raw_fifth(state, monomials, project):
    """Reduce Wk only for the highest-degree bilinear contraction."""
    data, R, rp, multiply = [state[k] for k in ['data', 'R', 'rp', 'multiply']]
    C2 = data['C2']
    p2 = list(C2.to_DM().rref()[1]); f2 = [j for j in range(66) if j not in p2]
    W2 = reduce_rows(state['W2'], rp(C2), p2, f2, R)
    p3, f3 = data['quotient_pivots'], data['free']
    W3 = reduce_rows(state['W3'], rp(data['quotient']), p3, f3, R)
    C3 = (state['cubic_inverse'].to_DM()*state['C3'][state['cubic_rows'], :].to_DM()).to_Matrix()[:, f3]
    pc3 = list(C3.to_DM().rref()[1]); fc3 = [j for j in range(len(f3)) if j not in pc3]
    W3 = reduce_rows(W3, rp(C3), pc3, fc3, R); f3 = [f3[j] for j in fc3]
    p4 = list(state['pivots']); f4 = [j for j in range(406) if j not in p4]
    H4 = reduce_rows(state['H4p'], rp(state['rr']), p4, f4, R)
    W4 = multiply(rp(-data['UP']), [H4[i] for i in data['rows']])
    f4 = [state['free4'][j] for j in f4]
    bases = {2: [monomials[2][j] for j in f2], 3: [monomials[3][j] for j in f3], 4: [monomials[4][j] for j in f4]}
    index = {mon: j for j, mon in enumerate(monomials[5])}

    def convolution_add(out, left, right, left_mons, right_mons):
        for i, x in enumerate(left):
            if x:
                for j, y in enumerate(right):
                    if y:
                        out[index[tuple(sorted(left_mons[i]+right_mons[j]))]] += x*y

    pairs = sorted({(min(i, j), max(i, j)) for row in data['terms'] for i, j, _ in row})
    cache = {}
    T = state['TT']
    for i, j in pairs:
        out = [R.zero]*3003
        for left, right, lm, rm in [(T[i], W4[j], monomials[1], bases[4]), (T[j], W4[i], monomials[1], bases[4]),
                                    (W2[i], W3[j], bases[2], bases[3]), (W2[j], W3[i], bases[2], bases[3])]:
            convolution_add(out, left, right, lm, rm)
        cache[i, j] = project(out)
    H5 = [[R.zero]*763 for _ in range(135)]
    for row, terms in enumerate(data['terms']):
        for i, j, coefficient in terms:
            for k, value in enumerate(cache[min(i, j), max(i, j)]):
                if value:
                    H5[row][k] += coefficient*value
    correction = multiply(state['JUP'], [H5[i] for i in data['rows']])
    return [[state['d']*x-y for x, y in zip(row, corr)] for row, corr in zip(H5, correction)]


def run():
    started = time.time(); stages = {}
    state = build_four()
    data, R, a, rp, multiply, d = [state[k] for k in ['data', 'R', 'a', 'rp', 'multiply', 'd']]
    monomials = {k: list(combinations_with_replacement(range(11), k)) for k in range(1, 6)}
    indices = {k: {mon: j for j, mon in enumerate(monomials[k])} for k in monomials}
    project, free5, quotient_report = fifth_quotient(state, monomials)
    stages['fourth_data_and_quotient'] = time.time()-started
    A5 = raw_fifth(state, monomials, project)
    stages['raw_fifth_projected'] = time.time()-started
    # IMPORTANT: retain the original FULL third and fourth tails below.
    correction = multiply(state['JUP'], [state['H4'][i] for i in data['rows']])
    A4 = [[d*x-y for x, y in zip(row, corr)] for row, corr in zip(state['H4'], correction)]
    A3, C2, M, inverse, selected = state['A3'], data['C2'], data['M'], data['selected_inverse'], data['selected']
    N = rp((M.to_DM()*inverse.to_DM()).to_Matrix())
    subtract = lambda left, right: [[x-y for x, y in zip(row, corr)] for row, corr in zip(left, right)]
    scale = lambda matrix, coefficient: [[coefficient*x for x in row] for row in matrix]

    def times_monomial(row, degree, mon, total, highest=False):
        out = [R.zero]*len(monomials[total])
        for k, value in enumerate(row):
            if value:
                out[indices[total][tuple(sorted(monomials[degree][k]+mon))]] += value
        return project(out) if highest else out

    G23 = multiply(rp(inverse), [A3[i] for i in selected]); G24 = multiply(rp(inverse), [A4[i] for i in selected])
    G33 = subtract(A3, multiply(N, [A3[i] for i in selected])); G34 = subtract(A4, multiply(N, [A4[i] for i in selected]))
    G35 = subtract(A5, multiply(N, [A5[i] for i in selected]))
    P3 = s.zeros(132, 286); P34 = []; P35 = []
    for i in range(12):
        for j in range(11):
            for k, mon in enumerate(monomials[2]):
                P3[11*i+j, indices[3][tuple(sorted(mon+(j,)))]] += C2[i, k]
            P34.append(times_monomial(G23[i], 3, (j,), 4))
            P35.append(times_monomial(G24[i], 4, (j,), 5, True))
    pr = list(P3.T.to_DM().rref()[1]); p3 = list(P3.to_DM().rref()[1])
    inverse3 = P3[pr, p3].inv(method='DM'); B3 = (inverse3.to_DM()*P3[pr, :].to_DM()).to_Matrix()
    B34 = multiply(rp(inverse3), [P34[i] for i in pr]); B35 = multiply(rp(inverse3), [P35[i] for i in pr])
    GP = [[row[j] for j in p3] for row in G33]
    G3red = subtract(G33, multiply(GP, rp(B3)))
    C3 = s.Matrix([[value.exquo(d*d).as_expr() for value in row] for row in G3red])
    cr = list(C3.T.to_DM().rref()[1]); rf, pc = C3[cr, :].to_DM().rref(); rf = rf.to_Matrix()
    assert len(pc) == 11 and not any(value.free_symbols for value in rf)
    inverse_cubic = C3[cr, list(pc)].inv(method='DM')
    C34raw = subtract(scale(G34, d), multiply(GP, B34)); C35raw = subtract(scale(G35, d), multiply(GP, B35))
    C34 = multiply(rp(inverse_cubic), [C34raw[i] for i in cr]); C35 = multiply(rp(inverse_cubic), [C35raw[i] for i in cr])
    G44 = subtract(C34raw, multiply(rp(C3[:, list(pc)]), C34)); G45 = subtract(C35raw, multiply(rp(C3[:, list(pc)]), C35))
    G44 += scale(subtract(P34, multiply(rp(P3[:, p3]), B34)), d*d)
    G45 += scale(subtract(P35, multiply(rp(P3[:, p3]), B35)), d*d)
    stages['full_lower_tail_elimination'] = time.time()-started
    # The complete 913 quartic multiples have 595 independent leading rows.
    P4 = s.zeros(913, 1001); P45 = []
    for i in range(12):
        for j, mon in enumerate(monomials[2]):
            for k, right in enumerate(monomials[2]):
                P4[66*i+j, indices[4][tuple(sorted(mon+right))]] += C2[i, k]
            P45.append([d*d*x for x in times_monomial(G23[i], 3, mon, 5, True)])
    for i in range(11):
        for j in range(11):
            for k, mon in enumerate(monomials[3]):
                P4[792+11*i+j, indices[4][tuple(sorted(mon+(j,)))]] += rf[i, k]
            P45.append(times_monomial(C34[i], 4, (j,), 5, True))
    r4 = list(P4.T.to_DM().rref()[1]); p4 = list(P4.to_DM().rref()[1]); assert len(r4) == len(p4) == 595
    inverse4 = P4[r4, p4].inv(method='DM'); B4 = (inverse4.to_DM()*P4[r4, :].to_DM()).to_Matrix()
    B45 = multiply(rp(inverse4), [P45[i] for i in r4]); product_syzygies = subtract(P45, multiply(rp(P4[:, p4]), B45))
    g4p = [[row[j] for j in p4] for row in G44]; G4red = subtract(G44, multiply(g4p, rp(B4)))
    free4 = [j for j in range(1001) if j not in p4]
    Q4 = s.Matrix([[row[j].as_expr() for j in free4] for row in G4red])
    previous = json.loads((ROOT/'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUARTIC_20261101.json').read_text())
    c4rows, c4piv = previous['new_quartic_selected_rows'], previous['new_quartic_pivots']
    assert free4 == previous['quotient_free_columns']
    rref4 = s.Matrix(previous['new_quartic_rref']).applyfunc(s.sympify)
    assert (Q4.to_DM()-Q4[:, c4piv].to_DM()*rref4.to_DM()).is_zero_matrix
    N5 = subtract(scale(G45, d**3), multiply(g4p, B45))
    N4 = (Q4[:, c4piv].to_DM()*Q4[c4rows, c4piv].inv(method='DM').to_DM()).to_Matrix().applyfunc(s.cancel)
    denominator = s.lcm([s.denom(value) for value in N4]); power = int(s.degree(denominator, a[1]))
    assert s.Poly(denominator, *a).monoms() == [(0, power, 0, 0, 0)]
    NP = N4.applyfunc(lambda value: s.cancel(value*a[1]**power))
    new5 = subtract(scale(N5, d**power), multiply(rp(NP), [N5[i] for i in c4rows])) + product_syzygies
    matrix = s.Matrix([[value.as_expr() for value in row] for row in new5])
    sample = dict(zip(a, range(1, 6))); selected5 = list(matrix.subs(sample).T.to_DM().rref()[1]); assert len(selected5) == 23
    rr, pivots = matrix[selected5, :].to_DM().rref(); rr = rr.to_Matrix().applyfunc(s.factor)
    assert len(pivots) == 23 and not any(value.free_symbols for value in rr)
    assert (matrix.to_DM()-matrix[:, list(pivots)].to_DM()*rr.to_DM()).is_zero_matrix
    stages['all_fifth_relations'] = time.time()-started
    inputs = [Path(__file__)] + [ROOT/name for name in [
        'code/s61_minimal_free_flag_filtered_data_20261101.py',
        'code/s61_minimal_free_flag_universal_quartic_20261101.py',
        'code/s61_minimal_free_flag_universal_cubic_20261101.py',
        'code/s61_model.py', 'code/s61_flag_transverse_obstruction_20261101.py',
        'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUADRATIC_20261101.json',
        'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_CUBIC_20261101.json',
        'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUARTIC_20261101.json',
        'reports/S61_MINIMAL_FREE_FLAG_NORMAL_FORM_20261101.json',
        'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json',
        'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']]
    return dict(id='S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUINTIC_20261101', all_passed=True,
                field='QQ(a0,a1,a2,a3,a4)', known_fifth_rank=2240, new_fifth_rank=23,
                true_initial_hilbert_function=[1, 11, 54, 168, 391, 740],
                fifth_monomials=[list(mon) for mon in monomials[5]], quotient_free_columns=free5,
                known_fifth_quotient=quotient_report, selected_rows=selected5, pivots=list(pivots),
                all_row_coefficients=enc(matrix[:, list(pivots)]), constant_fifth_rref=enc(rr),
                relation_rows=matrix.rows, stages_seconds=stages,
                scope=['Exact TRUE initial ideal through degree5 over the five-parameter function field.',
                       'All23 new fifth rows have rational constant coefficients; the earlier fourth row still depends on a1.',
                       'No full associated-graded equality, ambient t-nilpotence or saturation-removal conclusion.'],
                proof_note='research/notes/S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUINTIC_20261101.md',
                input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
                seconds=time.time()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); result = run()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2); stream.write('\n')
    print(json.dumps(dict(all_passed=True, seconds=result['seconds'], stages=result['stages_seconds'])))

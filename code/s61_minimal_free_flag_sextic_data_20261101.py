"""Exact lower-order spaces and native implicit-jet inputs for the sextic proof."""
import json
from itertools import combinations_with_replacement as cwr
import sympy as s
from s61_minimal_free_flag_filtered_data_20261101 import build_four, ROOT
from s61_minimal_free_flag_universal_quintic_20261101 import reduce_rows


def initial_generators(ring, monomials):
    v = s.symbols('v0:11')
    reports = [json.loads((ROOT/f'reports/S61_MINIMAL_FREE_FLAG_UNIVERSAL_{name}_20261101.json').read_text())
               for name in ['QUADRATIC', 'CUBIC', 'QUARTIC', 'QUINTIC']]
    generators = []
    for row in reports[0]['constant_quadratic_rref']:
        generators.append((2, {mon: ring.from_expr(s.Rational(c)) for mon, c in zip(monomials[2], row) if c != '0'}))
    for expression in reports[1]['constant_cubic_polynomials']:
        polynomial = s.Poly(s.sympify(expression), *v)
        generators.append((3, {tuple(j for j, n in enumerate(exponent) for _ in range(n)): ring.from_expr(c)
                               for exponent, c in polynomial.terms()}))
    for degree, report, key in [(4, reports[2], 'new_quartic_rref'), (5, reports[3], 'constant_fifth_rref')]:
        for row in report[key]:
            generators.append((degree, {monomials[degree][j]: ring.from_expr(s.sympify(c))
                                       for j, c in zip(report['quotient_free_columns'], row) if c != '0'}))
    return generators, reports


def fifth_product_basis(ring, monomials, generators):
    """Tracked RREF, keeping a basis of all degree-five product syzygies."""
    index = {mon: j for j, mon in enumerate(monomials[5])}
    products, metadata = [], []
    for i, (degree, polynomial) in enumerate(generators[:38]):
        for mon in monomials[5-degree]:
            products.append({index[tuple(sorted(mon+exponent))]: c for exponent, c in polynomial.items()})
            metadata.append((i, mon))
    assert len(products) == 4323
    variable = [i for i, (j, mon) in enumerate(metadata) if j == 23]
    assert len(variable) == 11
    basis, syzygies = {}, []
    def subtract(row, other, coefficient):
        for j, value in other.items():
            row[j] = row.get(j, ring.zero)-coefficient*value
            if not row[j]:
                del row[j]
    for i in [j for j in range(4323) if j not in variable]+variable:
        row, combination = products[i].copy(), {i: ring.one}
        while row:
            pivot = min(row)
            coefficient = row[pivot]
            if pivot not in basis:
                assert coefficient.is_ground
                row = {j: value.quo_ground(coefficient.LC) for j, value in row.items()}
                combination = {j: value.quo_ground(coefficient.LC) for j, value in combination.items()}
                basis[pivot] = row, combination
                break
            subtract(row, basis[pivot][0], coefficient)
            subtract(combination, basis[pivot][1], coefficient)
        if not row:
            syzygies.append(combination)
    pivots = sorted(basis)
    assert len(pivots) == 2240 and len(syzygies) == 2083
    for pivot in reversed(pivots):
        row, combination = basis[pivot]
        for earlier in pivots:
            if earlier >= pivot:
                break
            left, tracked = basis[earlier]
            coefficient = left.get(pivot, ring.zero)
            if coefficient:
                subtract(left, row, coefficient)
                subtract(tracked, combination, coefficient)
    free = [j for j in range(3003) if j not in basis]
    free_index = {j: i for i, j in enumerate(free)}
    normal_forms = []
    for j in range(3003):
        normal_forms.append([(free_index[j], ring.one)] if j in free_index else
                            [(free_index[k], -c) for k, c in basis[j][0].items() if k != j])
    return dict(products=products, metadata=metadata, pivots=pivots, free=free,
                basis=[basis[j][0] for j in pivots], combinations=[basis[j][1] for j in pivots],
                syzygies=syzygies, normal_forms=normal_forms)


def sixth_normal_forms(ring, monomials, generators):
    index = {mon: j for j, mon in enumerate(monomials[6])}
    entries, count = {}, 0
    for i, (degree, polynomial) in enumerate(generators):
        if i == 23:
            continue
        for mon in monomials[6-degree]:
            for exponent, coefficient in polynomial.items():
                assert coefficient.is_ground
                key = count, index[tuple(sorted(mon+exponent))]
                entries[key] = entries.get(key, 0)+coefficient.as_expr()
            count += 1
    matrix = s.SparseMatrix(count, 8008, entries)
    reduced, pivots = matrix.to_DM().rref()
    assert len(pivots) == 6757
    rows = reduced.to_sdm()
    free = [j for j in range(8008) if j not in pivots]
    free_index = {j: i for i, j in enumerate(free)}
    rational = reduced.domain.to_sympy
    variable = []
    for mon in monomials[2]:
        row = {index[tuple(sorted(mon+exponent))]: c.as_expr() for exponent, c in generators[23][1].items()}
        for i, pivot in enumerate(pivots):
            coefficient = row.get(pivot, 0)
            if coefficient:
                for j, value in rows[i].items():
                    row[j] = row.get(j, 0)-coefficient*rational(value)
        variable.append([row.get(j, 0) for j in free])
    vr, vp = s.Matrix(variable).to_DM().rref()
    vr = vr.to_Matrix()[:len(vp), :].applyfunc(s.factor)
    assert len(vp) == 14 and not any(s.denom(c).free_symbols for c in vr)
    varfree = [j for j in range(len(free)) if j not in vp]
    vi = {j: i for i, j in enumerate(varfree)}
    vrows = [[(j, ring.from_expr(c)) for j, c in enumerate(row) if c] for row in vr.tolist()]
    pivot_index = {j: i for i, j in enumerate(pivots)}
    normal_forms = []
    for j in range(8008):
        row = {free_index[j]: ring.one} if j in free_index else {
            free_index[k]: -ring.from_expr(rational(c)) for k, c in rows[pivot_index[j]].items() if k != j}
        for pivot, relation in zip(vp, vrows):
            coefficient = row.get(pivot, ring.zero)
            if coefficient:
                for k, value in relation:
                    row[k] = row.get(k, ring.zero)-coefficient*value
        normal_forms.append([(vi[k], value) for k, value in row.items() if value])
    return normal_forms, [free[j] for j in varfree]


def seed_data(backend, progress):
    state = build_four()
    progress('degree-four exact data')
    data, R, a, rp, multiply, d = [state[k] for k in ['data', 'R', 'a', 'rp', 'multiply', 'd']]
    ar = R.gens
    mm = {k: list(cwr(range(11), k)) for k in range(1, 7)}
    ix = {k: {mon: j for j, mon in enumerate(mm[k])} for k in mm}
    generators, reports = initial_generators(R, mm)
    product_basis = fifth_product_basis(R, mm, generators)
    progress('tracked degree-five product space')
    assert product_basis['free'] == reports[3]['quotient_free_columns']
    nf, free6 = sixth_normal_forms(R, mm, generators)
    C2, M, inverse, selected = data['C2'], data['M'], data['selected_inverse'], data['selected']
    N = rp((M.to_DM()*inverse.to_DM()).to_Matrix())
    A3 = state['A3']
    correction = multiply(state['JUP'], [state['H4'][i] for i in data['rows']])
    A4 = [[d*x-y for x, y in zip(row, corr)] for row, corr in zip(state['H4'], correction)]
    pairs = sorted({(min(i, j), max(i, j)) for row in data['terms'] for i, j, c in row})
    pair_index = {pair: j for j, pair in enumerate(pairs)}
    F = [[R.zero]*len(pairs) for _ in range(135)]
    for row, terms in enumerate(data['terms']):
        for i, j, coefficient in terms:
            F[row][pair_index[min(i, j), max(i, j)]] += coefficient
    correction = multiply(state['JUP'], [F[i] for i in data['rows']])
    effective = [[d*x-y for x, y in zip(row, corr)] for row, corr in zip(F, correction)]
    correction = multiply(N, [effective[i] for i in selected])
    canceled = [[x-y for x, y in zip(row, corr)] for row, corr in zip(effective, correction)]
    W4full = multiply(rp(-data['UP']), [state['H4'][i] for i in data['rows']])
    identity5 = [[(j, R.one)] for j in range(3003)]
    G35 = backend.contract(R, [state['TT'], state['W2'], state['W3'], W4full],
                          [mm[k] for k in range(1, 5)], 5, identity5, 3003, pairs, canceled)
    progress('full fifth tail')
    # For the highest fifth/sixth contraction only, quotient each lower Wk.
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
    bases = [mm[1], [mm[2][j] for j in f2], [mm[3][j] for j in f3], [mm[4][j] for j in f4]]
    stacked = backend.contract(R, [state['TT'], W2, W3, W4], bases, 5,
                               product_basis['normal_forms'], 763, pairs, F+effective)
    H5, A5p = stacked[:135], stacked[135:]
    R5 = s.Matrix(reports[3]['constant_fifth_rref']).applyfunc(s.Rational)
    p5 = reports[3]['pivots']; f5 = [j for j in range(763) if j not in p5]
    H5 = reduce_rows(H5, rp(R5), p5, f5, R)
    W5 = multiply(rp(-data['UP']), [H5[i] for i in data['rows']])
    free5 = product_basis['free']
    bases.append([mm[5][free5[j]] for j in f5])
    A6 = backend.contract(R, [state['TT'], W2, W3, W4, W5], bases, 6, nf, 1237, pairs, effective)
    progress('projected fifth and sixth tails')
    return {key: value for key, value in locals().items() if key in
            ['R', 'a', 'ar', 'd', 'mm', 'ix', 'nf', 'free6', 'C2', 'M', 'inverse', 'selected', 'N',
             'A3', 'A4', 'A5p', 'free5', 'G35', 'A6', 'product_basis']}

#!/usr/bin/env python3
"""True fourth initial relations on the seven-parameter flag chart, b8 != 0."""
import argparse
import hashlib
import json
import time
from itertools import combinations_with_replacement
from pathlib import Path

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_initial_quartic_from_jets_20261101 import initial_quartic

ROOT = Path(__file__).resolve().parents[1]


def run():
    start = time.time()
    low_path = ROOT / 'reports/S61_FLAG_ZERO12_LOW_JETS_20261101.json'
    low = json.loads(low_path.read_text())
    assert low['all_passed']
    for path, expected in low['input_sha256'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected, path
    old = json.loads((ROOT / 'reports/S61_FLAG_UNIFORM_OPEN_20261101.json').read_text())
    model = S61()
    b = s.symbols('b0:15')
    parameters = list(map(s.Symbol, low['parameters']))
    PR, *_ = ring(parameters, QQ)
    K = QQ.frac_field(*parameters)
    cv = lambda x: PR.from_expr(s.cancel(x))
    W = {int(d): [{tuple(mon): cv(s.sympify(c)) for mon, c in row}
                  for row in mat] for d, mat in low['lift'].items()}
    stored = {int(d): [{tuple(mon): cv(s.sympify(c)) for mon, c in row}
                       for row in mat] for d, mat in low['residual'].items()}
    H = s.Matrix(old['flag']['beta_kernel']).applyfunc(s.Rational)
    sub = {s.Symbol(k): s.sympify(v) for k, v in low['substitutions'].items()}
    one = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 1]
    two = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 2]
    q0 = s.Matrix(old['flag']['point61']).applyfunc(s.Rational)[two, :]
    fixed = two.index(37)
    chart = [i for i in range(45) if i != fixed]
    _, F = specialize(model, dict(zip(one, (H * s.Matrix(b)).subs(sub))), True)
    J = s.zeros(135, 45)
    terms = []
    for row, f in enumerate(F):
        tr = []
        for mon, c in f.items():
            if len(mon) == 2:
                i, j = mon
                if i != fixed:
                    J[row, chart.index(i)] += c * q0[j]
                if j != fixed:
                    J[row, chart.index(j)] += c * q0[i]
            elif len(mon) == 1:
                J[row, 44] += c * q0[mon[0]]
            if fixed not in mon:
                i, j = [chart.index(k) for k in mon] + [44] * (2 - len(mon))
                cc = cv(c)
                if cc:
                    tr.append((i, j, cc))
        terms.append(tr)
    rows, cols = low['pivot_rows'], low['pivot_columns']
    block = J[rows, cols].to_DM().to_field()
    assert s.cancel(block.domain.to_sympy(block.det())) == -s.Rational(1400832, 25)
    inverse = block.inv().to_Matrix()
    UP = [[PR.zero] * 36 for _ in range(45)]
    for i, c in enumerate(cols):
        UP[c] = [cv(v) for v in inverse[i, :]]
    JP = [[cv(v) for v in row] for row in J.tolist()]

    def clean(f):
        return {mon: c for mon, c in f.items() if c}

    def product(out, f, g, c):
        for a, ca in f.items():
            for d, cd in g.items():
                mon = tuple(sorted(a + d))
                out[mon] = out.get(mon, PR.zero) + c * ca * cd

    def combine(matrix, polys):
        result = []
        for row in matrix:
            out = {}
            for c, f in zip(row, polys):
                if c:
                    for mon, value in f.items():
                        out[mon] = out.get(mon, PR.zero) + c * value
            result.append(clean(out))
        return result

    assert all(not f for f in combine(JP, W[1]))
    residual = {}
    for degree in [2, 3, 4]:
        Q = []
        for tr in terms:
            out = {}
            for da in range(1, degree // 2 + 1):
                db = degree - da
                for i, j, c in tr:
                    product(out, W[da][i], W[db][j], c)
                    if da != db:
                        product(out, W[db][i], W[da][j], c)
            Q.append(clean(out))
        lift = combine([[-v for v in row] for row in UP], [Q[i] for i in rows])
        if degree < 4:
            assert lift == W[degree]
        else:
            W[degree] = lift
        correction = combine(JP, lift)
        residual[degree] = []
        for q, correction_row in zip(Q, correction):
            out = dict(q)
            for mon, value in correction_row.items():
                out[mon] = out.get(mon, PR.zero) + value
            residual[degree].append(clean(out))
        assert all(not residual[degree][i] for i in rows)
        if degree < 4:
            assert residual[degree] == stored[degree]
    assert not W[1][-1] and not W[2][-1]
    print(json.dumps(dict(stage='consistent_raw_four_jet', seconds=time.time() - start)), flush=True)

    mons = {d: list(combinations_with_replacement(range(9), d)) for d in [2, 3, 4]}

    def matrix(polys, degree):
        ix = {mon: i for i, mon in enumerate(mons[degree])}
        return s.SparseMatrix(len(polys), len(ix),
                              {(i, ix[mon]): c.as_expr() for i, f in enumerate(polys)
                               for mon, c in f.items()})

    result = initial_quartic(*[matrix(residual[d], d) for d in [2, 3, 4]], 9, K)
    assert result['quadratic_rank'] == 8 and result['cubic_initial_rank'] == 67
    L, E = result['quartic_lower'], result['quartic_extra']
    rr, piv = L.rref()
    _, lower_rows = L.transpose().rref()
    assert len(piv) == len(lower_rows) == 284
    lower_minor = K.to_sympy(L.extract(lower_rows, piv).det())
    lower_constant = s.cancel(lower_minor / b[8] ** 16)
    assert lower_constant.is_Rational and lower_constant != 0
    nonp = [j for j in range(495) if j not in piv]

    def project(X):
        return (X.extract(range(X.shape[0]), nonp)
                - X.extract(range(X.shape[0]), piv) * rr.extract(range(284), nonp))

    N = project(E)
    nr, pc = N.rref()
    _, extra_rows = N.transpose().rref()
    assert len(pc) == len(extra_rows) == 9
    extra_minor = K.to_sympy(N.extract(extra_rows, pc).det())
    extra_constant = s.cancel(extra_minor / b[8] ** 17)
    assert extra_constant.is_Rational and extra_constant != 0
    # Unlike the deeper-binary example, quadratic linear syzygies add no quartics here.
    assert N.extract(range(8 * 9), range(N.shape[1])).is_zero_matrix

    B3 = result['degree3_leading'].extract(result['cubic_basis_source_rows'], range(165))
    tail4 = result['degree3_tail'].extract(result['cubic_basis_source_rows'], range(495))
    t3 = matrix([W[3][-1]], 3).to_DM().convert_to(K)
    t4 = matrix([W[4][-1]], 4).to_DM().convert_to(K)
    p3 = result['cubic_pivots']
    ct = t3.extract([0], p3) * B3.extract(range(B3.shape[0]), p3).inv()
    assert ct * B3 == t3
    nt = project(t4 - ct * tail4)
    last_columns = [j for j in range(nt.shape[1]) if j not in pc]
    tclass = (nt.extract([0], last_columns)
              - nt.extract([0], pc) * nr.extract(range(9), last_columns))
    witness = next((j, c) for (i, j), c in tclass.to_dok().items()
                   if K.to_sympy(c).is_Rational and c)
    assert not tclass.is_zero_matrix

    # Specialization is valid on ALL b8 != 0: no hidden nuisance denominator.
    b8index = parameters.index(b[8])

    def denominator_audit(D):
        powers = set()
        for value in D.to_dok().values():
            denominator = value.denom
            assert len(denominator) == 1
            mon = next(iter(denominator))
            assert all(exponent == 0 for i, exponent in enumerate(mon) if i != b8index)
            powers.add(mon[b8index])
        return sorted(powers)

    audits = {name: denominator_audit(D) for name, D in [
        ('quadratic_row_reduction', result['quadratic_row_reduction']),
        ('cubic_row_reduction', result['cubic_row_reduction']),
        ('lower_rref', rr), ('extra_rref', nr), ('t3_relation', ct), ('t4_class', tclass)]}
    sparse = lambda D: dict(shape=list(D.shape), entries=[
        [i, j, str(K.to_sympy(c))] for (i, j), c in sorted(D.to_dok().items())])
    paths = [Path(__file__), low_path, ROOT / 'code/s61_flag_zero12_low_jets_20261101.py',
             ROOT / 'code/s61_initial_quartic_from_jets_20261101.py', ROOT / 'code/s61_model.py',
             ROOT / 'code/s61_flag_transverse_obstruction_20261101.py',
             ROOT / 'reports/S61_FLAG_UNIFORM_OPEN_20261101.json',
             ROOT / 'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_FLAG_ZERO12_QUARTIC_20261101', all_passed=True,
                parameters=list(map(str, parameters)), substitutions=low['substitutions'],
                open='b8 != 0', linear_rank=36, tangent_dimension=9,
                low_jets_reconstructed_from_raw_135=True,
                raw_fourth_residual_terms=sum(map(len, residual[4])),
                raw_fourth_lift_terms=sum(map(len, W[4])),
                quadratic_rank=8, cubic_initial_rank=67, quartic_lower_rank=284,
                quartic_additional_rank=9, quartic_initial_rank=293,
                hilbert_function_through_degree4=[1, 9, 37, 98, 202],
                length_mod_m5=347,
                lower_minor=dict(rows=list(lower_rows), columns=list(piv),
                                 constant=str(lower_constant), b8_power=16),
                extra_minor=dict(rows=list(extra_rows), columns=list(pc),
                                 constant=str(extra_constant), b8_power=17),
                quartic_quotient_columns=nonp, quartic_extra_basis=sparse(N.extract(extra_rows, range(N.shape[1]))),
                quadratic_linear_syzygies_add_no_quartic_classes=True,
                denominator_b8_powers=audits,
                t_order_exactly=4, t3_relation=sparse(ct), t4_class=sparse(tclass),
                t4_remaining_monomials=[list(mons[4][nonp[j]]) for j in last_columns],
                nonzero_constant_t4_witness=dict(monomial=list(mons[4][nonp[last_columns[witness[0]]]]),
                                                coefficient=str(K.to_sympy(witness[1]))),
                limitations=['This theorem covers b8 != 0; no fourth-order rank assertion on b8=0.',
                             'Only initial degrees through four are identified. No full initial ideal or t-nilpotence claim.',
                             'No entire boundary saturation removal or affine exhaustion follows.'],
                proof_note='research/notes/S61_FLAG_ZERO12_QUARTIC_20261101.md',
                seconds=time.time() - start,
                input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(all_passed=True, seconds=result['seconds'])))

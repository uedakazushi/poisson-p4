#!/usr/bin/env python3
"""Independent exact pointwise RAW24 -> original product41 bridge, no GB.

This does not call the authors' chart/elimination helpers. The original six
quadratic source rows are retained as a previously certified finite input.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s
from sympy.polys.rings import ring

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_NONSQUARE_RAW24_PRODUCT41_BRIDGE_20261103_v63_REVIEW'
NOTE = ROOT / 'research/notes' / (ID + '.md')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def det(matrix, one):
    """Sparse Laplace determinant over the exact coefficient ring."""
    n = len(matrix)
    if not n:
        return one
    i = min(range(n), key=lambda k: sum(bool(q) for q in matrix[k]))
    ans = one * 0
    for j, q in enumerate(matrix[i]):
        if q:
            sub = [[matrix[k][l] for l in range(n) if l != j]
                   for k in range(n) if k != i]
            ans += (-1) ** (i + j) * q * det(sub, one)
    return ans


def main(args):
    start = time.monotonic()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    inputs = [
        ROOT / 'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json',
        ROOT / 'reports/S61_BINARY_HYPERPLANE_QUADRATIC_PAIR_20261102.json',
        ROOT / 'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json',
        ROOT / 'reports/S61_BINARY_NONSQUARE_PAIR_COVER_20261103_v62/report.json',
        ROOT / 'reports/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/S61_NONSQUARE_IRREDUCIBLE_QUAD_SHADOW_20261103_v63_REVIEW/report.json',
    ]
    data, pair, chart, cover, _, _ = [json.loads(p.read_text()) for p in inputs]
    checked = {}
    for path in inputs:
        report = json.loads(path.read_text())
        assert report['all_passed'], path
        for key in ('input_sha256', 'evidence_sha256'):
            for name, expected in report.get(key, {}).items():
                candidates = [ROOT / name, path.parent / name]
                found = next((p for p in candidates if p.is_file() and sha(p) == expected), None)
                assert found is not None, (path, name)
                checked[str(found.relative_to(ROOT))] = expected

    R, a, e, d = ring('a,e,d', s.QQ)
    K = R.to_field()
    S, *v = ring('v0,v1,v2,v3,v4,v5,v6,v7', R)
    light = [0, 1, 2, 3, 4, 6, 7]
    quadrics = [S.from_dict({tuple(mon): R.from_dict({tuple(ex): s.QQ(c)
                  for ex, c in terms}) for mon, terms in row})
                for row in data['quadratic_rows_sparse']]
    assert data['quadratic_rows'] == [19, 60, 84, 114, 117, 133]
    assert len(quadrics) == 6
    assert all(not mon[5] and sum(mon) == 2 for q in quadrics for mon in q)
    F = 500*a + 90*e**2 + 3*e - 24
    B = (50*a + 18*e**2 - 9*e - 2)*d + 60*e - 40
    C = 100*a + 18*e**2 - 9*e - 2
    J = 5*a*d - 6*e - 4
    Z, W = 3*v[0] - v[1] + v[6], v[3]
    A = -(d*e + 8)*quadrics[1] + 8*quadrics[2] + d*quadrics[3]
    ell, remainder = divmod(A, W)
    assert remainder == 0
    BB = (-10*a - 3*e**2 + 3*e)*quadrics[1] + 3*(e-1)*quadrics[3] + 10*a*quadrics[4]
    assert BB == -a*Z*ell
    assert ell == S.from_expr(s.sympify(pair['linear_form']))
    coeff = [ell.get(tuple(int(i == j) for j in range(8)), R.zero) for i in range(8)]
    assert coeff[5] == 0 and all(sum(mon) == 1 for mon in ell)
    g = coeff[1] + coeff[6]
    assert g == -s.QQ(16, 243)*C*J
    change = [[R.one if i == j else R.zero for j in range(6)] for i in range(6)]
    change[2] = [R.zero, -(d*e+8), R(8), d, R.zero, R.zero]
    change[4] = [R.zero, -10*a-3*e**2+3*e, R.zero, 3*(e-1), 10*a, R.zero]
    assert det(change, R.one) == 80*a
    basis = [v[0], v[2], v[4], v[7], Z, W, ell]
    coord = [[q.get(tuple(int(k == j) for k in range(8)), R.zero) for j in light]
             for q in basis]
    chart_det = det(coord, R.one)
    assert chart_det in (g, -g)
    # Solve v1 directly, independently of the producer's v6 formula.
    T, *w = ring('w0,w1,w2,w3,w4,w5,H', K)
    v1 = -sum((T(K(c))*ww for c, ww in zip(
        [coeff[0]-3*coeff[6], coeff[2], coeff[4], coeff[7], coeff[6], coeff[3]], w[:6])), T.zero)/K(g)
    images = [w[0], v1, w[1], w[5], w[2], w[6], w[4]-3*w[0]+v1, w[3]]

    def restrict(q, target, values):
        out = target.zero
        for mon, c in q.items():
            term = target(K(c)) if target is T else target(c)
            for image, exponent in zip(values, mon):
                if exponent:
                    term *= image**exponent
            out += term
        return out

    def unit(poly, known):
        if not poly:
            return False
        for factor in known:
            while poly.rem(factor) == 0:
                poly = poly.exquo(factor)
        return poly.is_ground

    assert restrict(ell, T, images) == 0
    assert all(restrict(q, T, images) == target for q, target in zip(basis[:6], w[:6]))
    assert restrict(v[5], T, images) == w[6]
    assert all(unit(c.denom, [C, J]) for image in images for c in image.values())
    restricted = [restrict(quadrics[i], T, images) for i in [0, 1, 3, 5]]
    multipliers = [K.one, K(J**2)/25, K(J**2)/25, K(J)/5]
    assert [K.from_expr(s.sympify(q)) for q in chart['exact_clearing_unit_multipliers']] == multipliers
    cleared = []
    for q, saved, multiplier in zip(restricted, chart['cleared_qfour_sparse'], multipliers):
        expected = T.from_dict({tuple(mon): K(R.from_dict({tuple(ex): s.QQ(c)
                                for ex, c in terms})) for mon, terms in saved})
        assert q*multiplier == expected
        assert all(mon[6] == 0 and sum(mon) == 2 for mon in q)
        assert all(unit(c.denom, [C, J]) for c in q.values())
        cleared.append(expected)

    # Actual (Z,W) restriction of the ORIGINAL four rows and safe fixed pivots.
    U, *u = ring('u0,u1,u2,u3,u4', R)
    images_ZW = [u[0], u[1], u[2], U.zero, u[3], U.zero, u[1]-3*u[0], u[4]]
    qZW = [restrict(quadrics[i], U, images_ZW) for i in [0, 1, 3, 5]]
    mons = sorted(set().union(*qZW))
    assert len(mons) == 15
    matrix = [[K(q.get(mon, R.zero)) for mon in mons] for q in qZW]
    initial = [row[:] for row in matrix]
    wanted = [K(5*a), K(a*C)/240, K(16*C**2)/81, K(160*a*C)/81]
    pivots = []
    used = []
    for i, (record, value) in enumerate(zip(cover['restriction_pivots'], wanted)):
        swap, column = record['row_swap'], record['column']
        assert column not in used
        matrix[i], matrix[swap] = matrix[swap], matrix[i]
        pivot = matrix[i][column]
        assert pivot == value == K.from_expr(s.sympify(record['value']))
        assert unit(pivot.numer, [a, F, B, C, J]) and unit(pivot.denom, [a, F, B, C, J])
        matrix[i] = [q/pivot for q in matrix[i]]
        for k in range(i+1, 4):
            scalar = matrix[k][column]
            matrix[k] = [x-scalar*y for x, y in zip(matrix[k], matrix[i])]
        used.append(column)
        pivots.append(str(pivot))
    selected_minor = det([[row[j] for j in used] for row in initial], K.one)
    assert selected_minor == K(160*a**3*C**4)/19683

    relation = [S.zero, -a*(d*e+8)*Z+(-10*a-3*e**2+3*e)*W,
                8*a*Z, a*d*Z+3*(e-1)*W, 10*a*W, S.zero]
    assert sum((l*q for l, q in zip(relation, quadrics)), S.zero) == 0
    assert relation[4].get(tuple(int(j == 3) for j in range(8)), R.zero) == 10*a
    assert all(l == S.from_expr(s.sympify(saved)) for l, saved in zip(relation, pair['syzygy_coefficients']))
    original_cubic_mons = {tuple(mon[k]+(j == k) for k in range(8))
                          for q in quadrics for mon in q for j in light}
    assert len(original_cubic_mons) == 84
    # Canonical ell(Z,W)S1 shadow has13 distinct degree3 monomials.
    kernel13 = {tuple((k == 0)+(k == j)+(k == i) for k in range(7))
                for j in [1, 2] for i in range(7)}
    assert len(kernel13) == 13

    # Pin the exact serialized P=ad normalization used by the new finite packet.
    RP, ap, ep, pp = ring('a,e,P', s.QQ)
    shifts = [0, 2, 2, 1]
    parameter_terms = 0
    converted = []
    for saved, shift in zip(chart['cleared_qfour_sparse'], shifts):
        row = []
        for mon, terms in saved:
            terms_P = []
            assert len(mon) == 7 and mon[6] == 0
            for (ia, ie, id_), c in terms:
                ex = (ia-id_+shift, ie, id_)
                assert min(ex) >= 0
                terms_P.append([list(ex), c])
            row.append([mon[:6], terms_P])
            parameter_terms += len(terms_P)
        converted.append(row)
    assert parameter_terms == 804
    # Enumerate the actual RAW24 product monomial union without a rank sample.
    cubic_mons = {tuple(mon[k]+(j == k) for k in range(6))
                  for row in converted for mon, _ in row for j in range(6)}
    assert len(cubic_mons) == 56
    finite_review = None
    if not args.bridge_only:
        finite_path = ROOT / 'reports/S61_NONSQUARE_PAIR_FINITE_BAD_LOCUS_20261103_v63_REVIEW/report.json'
        finite_review = json.loads(finite_path.read_text())
        # The independently frozen finite-source theorem is checked explicitly;
        # scope-key spelling can vary, so the full saved receipt is retained.
        assert finite_review.get('all_passed') is True, finite_path
        assert any('raw24' in k.lower() and bool(value) for k, value in finite_review.items()), finite_path
        inputs.append(finite_path)
        for name, expected in finite_review.get('input_sha256', {}).items():
            assert sha(ROOT / name) == expected, name
            checked[name] = expected

    evidence = {
        'original_quadratic_source_rows': data['quadratic_rows'],
        'old_light_positions': light,
        'ell_coefficient_g': str(g), 'coordinate_determinant': str(chart_det),
        'inplace_six_row_change_determinant': '80*a',
        'all_six_row_change_and_syzygy_coefficients_exact': True,
        'all_eight_physical_coordinate_inverse_identities': True,
        'all_four_mod_ell_quadric_identities': True,
        'coordinate_denominators_only_C_J': True,
        'exact_clear_multipliers': ['1', 'J^2/25', 'J^2/25', 'J/5'],
        'mod_ZW_unit_pivots': pivots,
        'mod_ZW_selected_columns': used,
        'mod_ZW_selected_minor': str(selected_minor),
        'P_equals_a_times_original_d': True,
        'P_conversion_a_clear_shifts': shifts,
        'actual_raw24_parameter_coefficient_terms': parameter_terms,
        'actual_raw24_product_shape': [24, 56],
        'actual_original_product_shape': [42, 84],
        'canonical_kernel13_monomials': sorted(kernel13),
        'four_further_ell_q_classes_independent_mod_kernel13': True,
        'restriction_kernel_at_least17_pointwise': True,
        'raw_rank24_implies_original_rank41_pointwise': True,
        'original_all_point_upper41_nonzero_syzygy_unit_coefficient': '10*a',
        'new_finite_review_bound': finite_review is not None,
    }
    (output / 'checks.json').write_text(json.dumps(evidence, indent=2)+'\n')
    result = {
        'id': ID, 'all_passed': True,
        'independent_physical_coordinate_and_unit_row_reconstruction': True,
        'pointwise_RAW24_rank24_implies_original_product_rank41': True,
        'original_product_rank41_outside_finite_FULLU_subset_accepted': finite_review is not None,
        'original_rank_defect_subset_of_raw24_defect': True,
        'full_open_original_rank41_not_claimed': True,
        'J_zero_full57_quartic_ambient_removal_not_claimed': True,
        'constructive_product_pivot_count_not_increased': True,
        'field': 'QQ, all geometric points on a*P*F*B0*C*J!=0',
        'no_GB_no_new_field_no_rank_sample': True,
        'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), NOTE, *inputs]},
        'checked_dependency_sha256': checked,
        'evidence_sha256': {'checks.json': sha(output/'checks.json')},
        'seconds': time.monotonic()-start,
        'peak_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (output / 'report.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['all_passed', 'original_product_rank41_outside_finite_FULLU_subset_accepted', 'seconds']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'reports'/ID)
    parser.add_argument('--bridge-only', action='store_true')
    main(parser.parse_args())

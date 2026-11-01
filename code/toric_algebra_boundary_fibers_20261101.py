#!/usr/bin/env python3
"""Exact non-injectivity and image dimensions for three algebra-boundary types.

Uses universal parameter formulas and independently computed first-jet ranks.
This is about the algebra-model image, not surjectivity onto the T01 boundary.
"""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as s

from weight1_structure import COORDS, J, constraints, infinitesimal_action, sp_basis

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'reports/TORIC_ALGEBRA_BOUNDARY_20261101.json'


def flatten_beta(beta, y):
    return s.Matrix([s.diff(beta[i, j], y[k]) for (i, j), k in COORDS])


def rank_certificate(matrix):
    columns = list(matrix.rref()[1])
    rows = list(matrix[:, columns].T.rref()[1])
    minor = matrix.extract(rows, columns)
    return {'rank': len(columns), 'rows': rows, 'columns': columns,
            'determinant': str(minor.det()),
            'matrix': [[str(x) for x in row] for row in minor.tolist()]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    y = s.Matrix(s.symbols('y0:4'))
    e = s.eye(4)[:, 3]
    eta = J * e
    assert eta == s.eye(4)[:, 2]
    entries = s.symbols('b00 b01 b02 b11 b12 b22')
    b00, b01, b02, b11, b12, b22 = entries
    b = s.Matrix([[b00, b01, b02, 0], [b01, b11, b12, 0],
                  [b02, b12, b22, 0], [0, 0, 0, 0]])
    m_y = e * (b * y).T
    f = s.eye(4) + m_y
    p = (f * J * f.T).applyfunc(s.expand)
    assert (p - J - m_y * J - J * m_y.T).applyfunc(s.expand) == s.zeros(4)
    beta = (-J * (p - J) * J).applyfunc(s.expand)
    assert beta == eta * (b * y).T - (b * y) * eta.T
    parameter = s.Symbol('a')
    shifted_b = b + parameter * eta * eta.T
    shifted_f = s.eye(4) + e * (shifted_b * y).T
    assert (shifted_f * J * shifted_f.T - p).applyfunc(s.expand) == s.zeros(4)
    delta = b00 * b11 - b01 ** 2
    remainder = b11 * b02 ** 2 - 2 * b01 * b02 * b12 + b00 * b12 ** 2
    assert s.expand(b[:3, :3].det() - delta * b22 + remainder) == 0
    rank_two_value = remainder / delta
    assert s.cancel(b[:3, :3].det().subs(b22, rank_two_value)) == 0
    vector = flatten_beta(beta, y)
    assert s.Matrix.vstack(*constraints()) * vector == s.zeros(8, 1)
    base2 = {b00: 1, b01: 0, b02: 0, b11: 1, b12: 0, b22: 0}
    base3 = {**base2, b22: 1}
    value2 = vector.subs(base2)
    value3 = vector.subs(base3)
    assert value2 == value3
    orbit = s.Matrix.hstack(*(infinitesimal_action(value2, a) for a in sp_basis()))
    free2 = [b00, b01, b02, b11, b12]
    tangent2 = orbit.row_join(vector.subs(b22, rank_two_value).jacobian(free2).subs(base2))
    tangent3 = orbit.row_join(vector.jacobian(entries).subs(base3))
    rank2 = rank_certificate(tangent2)
    rank3 = rank_certificate(tangent3)
    assert rank2['rank'] == rank3['rank'] == 8

    gamma = s.eye(4)[:, 0]
    rank_one_b = gamma * gamma.T
    rank_one_beta = eta * (rank_one_b * y).T - (rank_one_b * y) * eta.T
    rank_one_vector = flatten_beta(rank_one_beta, y)
    rank1 = rank_certificate(s.Matrix.hstack(*(infinitesimal_action(rank_one_vector, a)
                                              for a in sp_basis())))
    assert rank1['rank'] == 5
    shifted_e = e + parameter * J * gamma
    assert (rank_one_b * shifted_e) == s.zeros(4, 1)
    shifted_eta = J * shifted_e
    assert (shifted_eta * (rank_one_b * y).T - (rank_one_b * y) * shifted_eta.T
            - rank_one_beta).applyfunc(s.expand) == s.zeros(4)

    # A1 x A4_sq: fixed e=e4 is enough after Sp4 transport; ell(e)=1.
    l0, l1, l2 = s.symbols('l0 l1 l2')
    ell = s.Matrix([l0, l1, l2, 1])
    t = (ell.T * y)[0]
    w = J * ell
    m = t * s.eye(4) + y * ell.T - 5 * t * e * ell.T
    q = 4 * ell * ell.T
    product_f = s.eye(4) + m - y * (q * y).T
    assert (product_f - (1 + t) * s.eye(4) - ((1 - 4 * t) * y - 5 * t * e) * ell.T).applyfunc(s.expand) == s.zeros(4)
    product_p = (product_f * J * product_f.T).applyfunc(s.expand)
    product_shift_f = (1 + t) * s.eye(4) + ((1 - 4 * t) * y - 5 * t * (e + parameter * w)) * ell.T
    assert (product_shift_f * J * product_shift_f.T - product_p).applyfunc(s.expand) == s.zeros(4)
    pfaffian = (product_p[0, 1] * product_p[2, 3] - product_p[0, 2] * product_p[1, 3]
                + product_p[0, 3] * product_p[1, 2])
    assert s.expand(pfaffian - (1 - 4 * t) * (1 + t) ** 4) == 0
    linear_p = product_p.applyfunc(lambda value: sum(s.diff(value, variable).subs(dict.fromkeys(y, 0))
                                                    * variable for variable in y))
    product_beta = (-J * linear_p * J).applyfunc(s.expand)
    expected_beta = 2 * t * J + (J * y) * ell.T - ell * (J * y).T + 5 * t * (ell * eta.T - eta * ell.T)
    assert (product_beta - expected_beta).applyfunc(s.expand) == s.zeros(4)
    product_vector = flatten_beta(product_beta, y)
    assert s.Matrix.vstack(*constraints()) * product_vector == s.zeros(8, 1)
    sample = {l0: 1, l1: 2, l2: 3}
    product_tangent = s.Matrix.hstack(*(infinitesimal_action(product_vector.subs(sample), a)
                                       for a in sp_basis())).row_join(product_vector.jacobian((l0, l1, l2)).subs(sample))
    product_rank = rank_certificate(product_tangent)
    assert product_rank['rank'] == 6

    # Identify the previously saved multiplication tensors with these formulas.
    catalog = json.loads(SOURCE.read_text())
    catalog_checks = []
    for record in catalog['all_cases']:
        name = record['algebra_type']
        if name not in ('A1 x A4_sq', 'A5_31_r1', 'A5_31_r2', 'A5_31_r3'):
            continue
        regular = [s.Matrix(matrix) for matrix in record['trace_zero_regular_multiplication']]
        lower = [matrix[1:5, 1:5] for matrix in regular]
        trace_form = s.Matrix.vstack(*(matrix[0, 1:5] for matrix in regular))
        if name == 'A1 x A4_sq':
            i = next(i for i in range(4) if trace_form[i, i])
            li = -s.trace(regular[i] ** 3) / (15 * trace_form[i, i])
            covector = trace_form[i, :].T / (4 * li)
            distinguished = (2 * li * s.eye(4)[:, i] - lower[i][:, i]) / (5 * li ** 2)
            assert (covector.T * distinguished)[0] == 1
            assert trace_form == 4 * covector * covector.T
            for k in range(4):
                expected = covector[k] * s.eye(4) + s.eye(4)[:, k] * covector.T - 5 * covector[k] * distinguished * covector.T
                assert lower[k] == expected
            catalog_checks.append({'type': name, 'ell': list(map(str, covector)),
                                   'e': list(map(str, distinguished)), 'all_multiplication_entries_match': True})
        else:
            columns = s.Matrix.hstack(*lower)
            assert columns.rank() == 1 and trace_form == s.zeros(4)
            distinguished = columns[:, next(i for i in range(16) if any(columns[:, i]))]
            i = next(i for i in range(4) if distinguished[i])
            bilinear = s.Matrix.vstack(*(matrix[i, :] / distinguished[i] for matrix in lower))
            assert bilinear == bilinear.T and bilinear * distinguished == s.zeros(4, 1)
            assert all(lower[k] == distinguished * bilinear[k, :] for k in range(4))
            expected_rank = int(name[-1])
            assert bilinear.rank() == expected_rank
            catalog_checks.append({'type': name, 'e': list(map(str, distinguished)),
                                   'B': [[str(x) for x in row] for row in bilinear.tolist()],
                                   'B_rank': expected_rank, 'all_multiplication_entries_match': True})

    report = {
        'all_passed': True, 'field': 'QQ for identities and ranks; C for geometric images',
        'local_31_eta_convention': 'eta=J e=Omega(e,-), where Omega=-J',
        'local_31_beta_formula': 'beta(y)=eta wedge (B y), B=B^t, B e=0',
        'local_31_B_shift_preserves_entire_P': True,
        'rank2_representative_B22': str(rank_two_value), 'rank2_chart_Delta': str(delta),
        'rank2_image_dimension': 8, 'rank3_image_dimension': 8,
        'rank2_and_rank3_image_closures_equal': True,
        'rank2_differential': rank2, 'rank3_differential': rank3,
        'rank1_image_dimension': 5, 'rank1_image': 'the closed primitive rank-one flag cone',
        'rank1_generic_fiber': 'A1', 'rank1_differential': rank1,
        'product_A1_A4sq_m': 'm(u,v)=ell(u)v+ell(v)u-5ell(u)ell(v)e, ell(e)=1',
        'product_A1_A4sq_q': '4ell tensor ell',
        'product_A1_A4sq_norm': '(1-4ell(y))*(1+ell(y))^4',
        'product_A1_A4sq_fiber': 'e+a J ell, a in C, with ell fixed',
        'product_A1_A4sq_image_dimension': 6, 'product_A1_A4sq_differential': product_rank,
        'catalog_multiplication_identifications': catalog_checks,
        'scope': 'Dimensions and explicit non-injectivity of the indicated algebra-model images. '
                 'No assertion that the algebra model covers the entire T01 boundary, no ambient '
                 'S61 local-ring classification, and no classification of projective Poisson orbits.',
        'source_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in (Path(__file__), SOURCE, ROOT / 'code/weight1_structure.py')},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'image_dimensions': {'A1 x A4_sq': 6,
                      'A5_31_r1': 5, 'A5_31_r2': 8, 'A5_31_r3': 8},
                      'rank2_rank3_image_closures_equal': True}))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Independent differentiation and interpolation review of the C9 pencil proof."""
import argparse
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path
import time

import sympy as s

from weight1_structure import J, constraints, form_matrices

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'reports/TORIC_C9_POLYNOMIAL_20261101.json'
TRIPLES = list(combinations_with_replacement(range(4), 3))


def differentiated_B(vector):
    beta = form_matrices(vector)
    base = [s.Matrix(4, 4, lambda a, j: sum(J[a, k] * (beta[i][j, k] + beta[j][i, k]) / 3
                                           for k in range(4))) for i in range(4)]
    columns = []
    for triple in TRIPLES:
        direction = [s.Matrix(4, 4, lambda a, j: sum(J[a, k] * int(tuple(sorted((i, j, k))) == triple)
                                                    for k in range(4))) for i in range(4)]
        trace_derivative = s.Matrix(4, 4, lambda i, j: s.trace(direction[i] * base[j]
                                                                + base[i] * direction[j]) / 3)
        values = []
        for i, j in combinations(range(4), 2):
            matrix = (direction[i] * base[j] + base[i] * direction[j] - direction[j] * base[i]
                      - base[j] * direction[i] + s.eye(4)[:, i] * trace_derivative[j, :]
                      - s.eye(4)[:, j] * trace_derivative[i, :])
            values.extend(list(matrix))
        columns.append(s.Matrix(values))
    return s.Matrix.hstack(*columns)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    source = json.loads(SOURCE.read_text())
    basis = s.Matrix.hstack(*s.Matrix.vstack(*constraints()).nullspace())
    base = differentiated_B(basis * s.Matrix(source['base_weight_coordinates']))
    direction = differentiated_B(basis * s.Matrix(source['direction_weight_coordinates']))
    assert base == s.Matrix(source['B_base'])
    assert direction == s.Matrix(source['B_direction'])
    t = s.Symbol('t')
    polynomials = []
    nodes = list(range(-10, 11))
    for entry in source['maximal_minors']:
        polynomial = s.Poly(s.sympify(entry['polynomial']), t, domain=s.QQ)
        assert polynomial.degree() == 20
        rows = entry['rows']
        for point in nodes:
            determinant = (base + point * direction)[rows, :].det(method='domain-ge')
            assert determinant == polynomial.eval(point)
        # A 20-by-20 matrix pencil has determinant degree at most20; these
        # 21 exact values therefore verify the entire stored determinant.
        assert direction[rows, :].det(method='domain-ge') == polynomial.LC() != 0
        polynomials.append(polynomial)
    left, right = [s.Poly(s.sympify(value), t, domain=s.QQ) for value in source['bezout_coefficients']]
    assert left * polynomials[0] + right * polynomials[1] == s.Poly(1, t, domain=s.QQ)
    report = {
        'all_passed': True,
        'B_reconstruction': 'Direct derivative of commutator and trace-form equations; no polarization code imported',
        'both_96_by20_matrices_match': True,
        'determinants_reviewed': 2, 'exact_nodes_for_each_degree20_determinant': nodes,
        'degree_bound_justifies_interpolation': True,
        'rational_Bezout_identity_verified': True, 'rank20_at_infinity_verified': True,
        'geometric_review': 'The projective pencil avoids the homogeneous rank-drop locus, '
                            'so it has no codimension1 component. Unique rational V is regular '
                            'in every height1 local ring and hence polynomial on normal affine W16.',
        'finite_actual_projection_over_Delta_nonzero_review': 'Accepted: B_Delta is an R_Delta-submodule '
                            'of the finite module A_Delta, hence finite because R_Delta is Noetherian.',
        'scope': 'Global polynomiality of C9 and finiteness of actual T01 over the smooth-curve open. '
                 'Does not prove equality of its coordinate ring with the normal cubic ring.',
        'seconds': time.time() - started,
        'source_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in (Path(__file__), SOURCE, ROOT / 'code/weight1_structure.py')},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'both_B_matrices_match': True,
                      'exact_determinant_evaluations': 42, 'seconds': report['seconds']}))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Check the grading and specialization inputs for the boundary theorem.

Normality, non-finiteness and failure of boundary surjectivity are proved in
the accompanying note; these finite checks do not replace that argument.
"""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    names = ['TORIC_UNIVERSAL_THREE_PLANE_20261101', 'TORIC_JACOBIAN_TWO_TORSION_20261101',
             'TORIC_DEGREE3_SCALAR_20261101', 'PRIMITIVE_CURVE_DISCRIMINANT_REVIEW_20261101',
             'PRIMITIVE_CURVE_DOMINANT_SLICE_20261101']
    paths = [ROOT / f'reports/{name}.json' for name in names]
    reconstruction, cover, scalar, discriminant, sliced = [json.loads(path.read_text()) for path in paths]
    beta = s.symbols(' '.join(reconstruction['beta_variables']))
    center = [s.sympify(value) for value in reconstruction['q_center']]
    matrix = s.Matrix(reconstruction['q_linear_matrix_N'])
    assert len(center) == 45 and matrix.shape == (45, 20)
    assert all(s.Poly(value, *beta).total_degree() == 2 for value in center if value)
    assert all(s.Poly(value, *beta).total_degree() == 1 for value in matrix if value)
    origin = dict.fromkeys(beta, 0)
    assert all(value.subs(origin) == 0 for value in center)
    assert matrix.subs(origin) == s.zeros(45, 20)
    assert all(sum(monomial) == 3 for entry in scalar['D3_sparse'] for monomial, _ in entry)
    assert all(sum(monomial) == 8 for monomial, _ in sliced['c4_sparse'])
    assert all(sum(monomial) == 12 for monomial, _ in sliced['c6_sparse'])
    assert discriminant['cubic_discriminant_equals_136048896_Delta']
    assert cover['generic_cubic'] == 'X^3-27*c4(beta)*X-54*c6(beta)'
    t = s.Symbol('t')
    hilbert_y = (1 + t ** 4 + t ** 8) / (1 - t) ** 16
    assert s.expand(s.series(hilbert_y, t, 0, 3).removeO()) == 1 + 16 * t + 136 * t ** 2
    assert s.cancel(hilbert_y - (1 - t ** 12) / ((1 - t) ** 16 * (1 - t ** 4))) == 0
    report = {
        'all_passed': True,
        'weights': {'beta': 1, 'S61_q': 2, 'algebra_c': 1, 'D3': 3, 'X': 4, 'c4': 8, 'c6': 12},
        'normal_cubic_Hilbert_series': '(1+t^4+t^8)/(1-t)^16',
        'normal_cubic_degree2_equals_base_degree2': True,
        'base_degree2_dimension': 136,
        'universal_reconstruction_q_center_degree': 2,
        'universal_reconstruction_N_degree': 1,
        'every_finite_algebra_model_point_over_beta0_maps_to_q0': True,
        'proof_conclusions': [
            'The monogenic cubic Y is the normal finite flat degree3 cover with function field C(T01).',
            'Its reduced branch divisor is the irreducible degree24 canonical discriminant.',
            'Generic inertia along the branch divisor is a transposition.',
            'The actual reduced T01 projection to W16 is not finite.',
            'The actual beta=0 fiber of T01 has positive dimension.',
            'The finite-dimensional commutative-algebra model does not cover the whole T01 boundary.'],
        'proof_note': 'research/notes/TORIC_NORMALIZATION_BOUNDARY_20261101.md',
        'limitations': ['No description of every missing boundary point or its dimension beyond the positive lower bound.',
                        'No assertion that the actual T01 projection is finite or proper away from the discriminant.',
                        'No exhaustion of the full S61 generic fiber or full Poisson component space.'],
        'source_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in [Path(__file__), *paths]},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'grading_checks': True, 'beta0_algebra_image_q0': True}))


if __name__ == '__main__':
    main()

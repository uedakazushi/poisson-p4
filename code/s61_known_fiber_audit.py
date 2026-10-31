#!/usr/bin/env python3
"""Certify smooth, unramified S61 points on the three known dominant components.

This uses all 134 archived equations and all 45 weight-two variables. It is
independent of the first-jet parameter Jacobians; no exhaustion is asserted.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import sympy as s

from elliptic_weight1_dominance import sample_data
from s61_model import S61, pushforward, toric_family

ROOT = Path(__file__).resolve().parents[1]


def equation_jacobian(equations, point):
    """Differentiate the sparse encoded polynomials over QQ, including repeats."""
    rows = []
    for equation in equations:
        row = [s.S.Zero] * len(point)
        for monomial, numerator, denominator in equation:
            for column, multiplicity in Counter(monomial).items():
                remaining = list(monomial)
                remaining.remove(column)
                row[column] += (s.Rational(numerator * multiplicity, denominator)
                                * s.prod(point[i] for i in remaining))
        rows.append(row)
    return s.Matrix(rows)


def audit():
    model = S61()
    columns = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 2]
    assert len(columns) == 45
    results = []
    for label in ('T01', 'T39', 'T40'):
        if label == 'T01':
            bivector = toric_family(2, 3, 5, 7, 11, 13)[0]
            sample = {'log_poisson_parameters': [2, 3, 5, 7, 11, 13]}
        else:
            x, lam, original, z, affine, point, constant, derivatives, frame, linear = sample_data(
                1 if label == 'T39' else 2)
            bivector = pushforward(original.subs(lam, 2), s.diag(1, frame))
            sample = {'lambda': 2, 'projective_point': [1, 0, 0, 0, 0]}
        normalized = model.normalize(bivector)
        verification = model.verify_point(normalized)
        point = normalized['parameters']
        assert all(value.is_Rational for value in point)
        jacobian = equation_jacobian(model.reduced['equations'], point)
        fiber_jacobian = jacobian[:, columns]
        rows = list(fiber_jacobian.T.rref()[1])
        minor = fiber_jacobian[rows, :]
        # DomainMatrix avoids expensive symbolic zero testing on rational data.
        determinant = str(minor.to_DM().det())
        rank = len(jacobian.to_DM().rref()[1])
        assert len(rows) == rank == 45 and s.Rational(determinant) != 0
        results.append({
            'component_id': label, 'sample': sample,
            'S61_coordinates_u1_through_u61': [str(value) for value in point],
            'point_verification': verification,
            'full_equation_jacobian_rank': rank,
            'full_Zariski_tangent_dimension': 61 - rank,
            'fixed_first_jet_tangent_dimension': 45 - len(rows),
            'fiber_minor': {
                'rows_zero_based': rows, 'u_columns_zero_based': columns,
                'determinant': determinant,
                'matrix': [[str(value) for value in row] for row in minor.tolist()],
            },
        })
    return {
        'id': 'S61_KNOWN_FIBER_AUDIT_20261031', 'all_passed': True, 'field': 'QQ',
        'equation_count': 134, 'weight_two_variables': 45, 'points': results,
        'interpretation': 'Each specified point is smooth of dimension 16 on the full S61 '
                          'scheme and the first-jet projection is etale there, using '
                          'the established 16-dimensional component through that point.',
        'limitations': ['No assertion about other points of these fibers.',
                        'No generic degree, branch locus, or global component exhaustion.',
                        'Unknown components may still dominate W16.'],
        'proof_note': 'research/notes/S61_KNOWN_FIBER_AUDIT_20261031.md',
        'input_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in (Path(__file__), ROOT/'code/s61_model.py',
                                      ROOT/'code/elliptic_weight1_dominance.py',
                                      model.source/'symplectic_slice_reduced_equations.json',
                                      model.source/'symplectic_slice_equations.json',
                                      model.source/'symplectic_slice_certificate.json')},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f'Refusing to overwrite {args.output}')
    report = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False)
        stream.write('\n')
    print(json.dumps({'all_passed': True, 'full_ranks': [45, 45, 45],
                      'fiber_tangent_dimensions': [0, 0, 0], 'output': str(args.output)}))


if __name__ == '__main__':
    main()

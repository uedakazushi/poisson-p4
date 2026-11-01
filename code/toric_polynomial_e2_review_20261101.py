#!/usr/bin/env python3
"""Independent ordered-index contraction review of the polynomial E2 coordinate."""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from s61_model import S61, X
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import J

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'reports/TORIC_POLYNOMIAL_E2_COORDINATE_20261101.json'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(SOURCE.read_text())
    model = S61()
    variables, u, beta = primitive_beta(model)
    y = s.Matrix(variables)
    indices = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 2]
    q = s.symbols('q0:45')
    parameters = [q[indices.index(i)] if i in indices else 0 for i in range(61)]
    homogeneous = model.reconstruct(parameters)
    substitutions = dict(zip(X, [1, *y], strict=True))
    affine = s.Matrix(4, 4, lambda i, j: s.expand((homogeneous[i + 1, j + 1]
                      - X[i + 1] * homogeneous[0, j + 1]
                      + X[j + 1] * homogeneous[0, i + 1]).subs(substitutions)))
    delta = -J * affine * J
    vector_field = -J * beta * y
    jacobian = vector_field.jacobian(y)
    eta1 = jacobian.T * J * jacobian
    pfaffian = beta[0, 1] * beta[2, 3] - beta[0, 2] * beta[1, 3] + beta[0, 3] * beta[1, 2]
    alpha = -pfaffian * J * y
    eta2 = s.Matrix(4, 4, lambda i, j: s.diff(alpha[j], y[i]) - s.diff(alpha[i], y[j]))
    combined = -18 * eta1 - 9 * eta2
    all_indices = list(product(range(4), repeat=4))
    zero_y = dict.fromkeys(y, 0)
    # The coefficient tensor equals half the second derivative, for both
    # repeated and distinct polynomial indices. Keep all ordered form indices.
    tensor = {(i, j, a, b): s.expand(s.diff(combined[i, j], y[a], y[b]) / 2)
              for i, j, a, b in all_indices}
    target = {(i, j, a, b): s.expand(s.diff(delta[i, j], y[a], y[b]).subs(zero_y) / 2)
              for i, j, a, b in all_indices}
    total = s.S.Zero
    for i, j, a, b in all_indices:
        # J has one nonzero entry per row. This is the full ordered-index sum,
        # with no compressed antisymmetric/symmetric multiplicity factors.
        total += (J[i, i ^ 1] * J[j, j ^ 1] * J[a, a ^ 1] * J[b, b ^ 1]
                  * tensor[i, j, a, b] * target[i ^ 1, j ^ 1, a ^ 1, b ^ 1])
    row = s.Matrix([s.expand(total)]).jacobian(q)
    expected = s.Matrix([list(map(s.sympify, source['lambda_coefficients_q']))])
    assert (row - expected).applyfunc(s.expand) == s.zeros(1, 45)
    assert all(s.Poly(value, *u).total_degree() == 2 for value in row if value)

    # Independently re-evaluate the stored residual contractions on the full
    # dominant six-variable slice, using the just-verified scalar coefficients.
    c5_path = ROOT / 'reports/TORIC_DEGREE5_COVARIANT_20261101.json'
    slice_path = ROOT / 'reports/PRIMITIVE_CURVE_DOMINANT_SLICE_20261101.json'
    c9_path = ROOT / 'reports/TORIC_JACOBIAN_TWO_TORSION_20261101.json'
    c5 = json.loads(c5_path.read_text())
    sliced = json.loads(slice_path.read_text())
    c9 = json.loads(c9_path.read_text())
    r, *sv = ring('s0,s1,s2,s3,s4,s5', QQ)
    change = s.Matrix(c5['b_in_u']).inv() * s.Matrix(sliced['slice_matrix_b'])
    images = [sum((QQ.from_sympy(change[i, j]) * sv[j] for j in range(6)), r.zero) for i in range(16)]
    residual = []
    for expression in source['residual_row_lambda_N_minus_D']:
        value = r.zero
        for monomial, coefficient in s.Poly(s.sympify(expression), *u).terms():
            term = r(QQ.from_sympy(coefficient))
            for i, degree in enumerate(monomial):
                if degree:
                    term *= images[i] ** degree
            value += term
        residual.append(value)
    def decode(entries):
        return r.from_dict({tuple(mon): QQ.from_sympy(s.Rational(co)) for mon, co in entries})
    for polynomials in (sliced['C5_sparse'], c9['C9_sparse']):
        assert not sum((a * decode(b) for a, b in zip(residual, polynomials, strict=True)), r.zero)

    report = {
        'all_passed': True,
        'independent_pairing_method': 'Half second derivatives and all256 ordered index tuples; no compressed metric helper imported',
        'all45_lambda_polynomial_coefficients_match': True,
        'scalar_weight': 4, 'degree_in_q': 1,
        'both_six_variable_residual_contractions_reverified': True,
        'proof_review': 'Use covariance of the two intrinsic scalars on toric solutions, '
                        'not off-shell equivariance of the row-reduced N. Slice equality holds '
                        'in the full cubic algebra and extends by dominance and reduced closure.',
        'global_consequences_with_polynomial_C9': ['A embeds into B through the polynomial X coordinate.',
                        'B localized at Delta equals the normal cubic A localized at Delta.',
                        'Actual T01 over Delta nonzero is the finite etale degree3 Jacobian two-torsion cover.',
                        'Globally B is an affine modification of A along explicit discriminant-supported numerators.'],
        'proof_note': 'research/notes/TORIC_AFFINE_MODIFICATION_20261101.md',
        'scope': 'The reduced known component T01, not the full S61 scheme or unknown components.',
        'source_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in (Path(__file__), SOURCE, c5_path, slice_path, c9_path,
                                      ROOT / 'code/s61_model.py', ROOT / 'code/primitive_jet_elliptic_curve_20261101.py')},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'full_ordered_tensor_pairing_verified': True,
                      'all45_lambda_match': True, 'both_residual_contractions_zero': True}))


if __name__ == '__main__':
    main()

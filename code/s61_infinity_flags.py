#!/usr/bin/env python3
"""Rank-one flags at infinity of the UNSATURATED homogenized S61 fibers.

This verifies the representative identity, the incidence-section injectivity,
and an additional rational flag with all16 first-jet coordinates free. It does
not identify the true projective closure after homogenizing-variable saturation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as s

from elliptic_weight1_dominance import projective_affine
from s61_model import S61, X, J as J5
from weight1_structure import J, constraints, form_matrices, minor_witness

ROOT = Path(__file__).resolve().parents[1]


def string_matrix(matrix):
    return [[str(value) for value in row] for row in matrix.tolist()]


def incidence_injectivity():
    v = s.Matrix(s.symbols('v1:5'))
    w = s.Matrix(s.symbols('w1:5'))
    omega = (v.T*J*w)[0]
    kernel = s.Matrix.vstack(*constraints()).nullspace()
    polynomials = []
    for vector in kernel:
        beta_at_v = sum((v[k]*b for k, b in enumerate(form_matrices(vector))), s.zeros(4))
        polynomials.append(s.Poly(s.expand((v.T*beta_at_v*w)[0]), *v, *w))
    relations = [s.Poly(s.expand(omega*x), *v, *w) for x in v]
    monomials = sorted(set().union(*(set(p.monoms()) for p in polynomials+relations)))
    matrix = s.Matrix([[p.coeff_monomial(mon) for p in polynomials+relations]
                       for mon in monomials])
    witness = minor_witness(matrix)
    assert witness['rank'] == 20
    assert witness['determinant'] == '625'
    # Perspective invariance of the incidence condition modulo omega(v,w).
    alpha = s.Matrix(s.symbols('a1:5'))
    tangent = -2*alpha.dot(v)*J + v*(J*alpha).T - (J*alpha)*v.T
    beta_perspective = -J*tangent*J
    assert s.expand((v.T*beta_perspective*w)[0]+3*alpha.dot(v)*omega) == 0
    return {
        'primitive_closed_basis': string_matrix(s.Matrix.hstack(*kernel)),
        'coefficient_variables': [str(x) for x in (*v, *w)],
        'coefficient_monomials': [list(mon) for mon in monomials],
        'column_order': '16 primitive closed beta basis vectors; omega(v,w)*v1,...,v4',
        'coefficient_matrix': string_matrix(matrix),
        'rank20_witness': witness,
        'conclusion': 'No nonzero primitive closed beta gives a section identically zero on isotropic flags.',
        'perspective_identity': 'T_(beta+beta_Ta)(v,w)-T_beta(v,w)=-3*alpha(v)*omega(v,w)',
    }


def leading_equations_matrix(model, coordinates):
    weights = model.reduced['parameter_weights']
    one = [i for i, weight in enumerate(weights) if weight == 1]
    two = [i for i, weight in enumerate(weights) if weight == 2]
    linear_index = {old: i for i, old in enumerate(one)}
    rows = []
    for equation in model.reduced['equations']:
        row = [s.S.Zero]*17
        for monomial, numerator, denominator in equation:
            quadratic_indices = [i for i in monomial if i in two]
            if len(quadratic_indices) != 2:
                continue
            linear_indices = [i for i in monomial if i in linear_index]
            assert len(linear_indices) <= 1
            coefficient = s.Rational(numerator, denominator)*s.prod(
                coordinates[i] for i in quadratic_indices)
            column = 1+linear_index[linear_indices[0]] if linear_indices else 0
            row[column] += coefficient
        rows.append(row)
    return s.Matrix(rows)


def old_beta_basis(model):
    z = s.symbols('z1:5')
    origin = dict(zip(z, [0]*4))
    result = []
    for old, weight in enumerate(model.reduced['parameter_weights']):
        if weight != 1:
            continue
        parameters = [s.S.Zero]*61
        parameters[old] = s.S.One
        affine = projective_affine(model.reconstruct(parameters), X, z)
        linear = sum((affine.diff(zi).subs(origin)*zi for zi in z), s.zeros(4))
        result.append(-J*linear*J)
    return z, result


def flag_check(model, v, w, beta_basis):
    assert (v.T*J*w)[0] == 0 and v.row_join(w).rank() == 2
    ell = -J*v
    quadratic = ell.dot(s.Matrix(X[1:5]))**2*(v*w.T-w*v.T)
    bivector = X[0]**2*J5
    bivector[1:5, 1:5] += quadratic
    normalized = model.normalize(bivector)
    model.verify_point(normalized)
    coordinates = normalized['parameters']
    matrix = leading_equations_matrix(model, coordinates)
    assert matrix[:, 0] == s.zeros(134, 1)
    z, betas = beta_basis
    incidence = s.Matrix([[s.expand((v.T*beta.subs(dict(zip(z, v)))*w)[0])
                           for beta in betas]])
    assert matrix[:, 1:].rank() == matrix[:, 1:].col_join(incidence).rank() == 1
    return {
        'v': [str(x) for x in v], 'w': [str(x) for x in w],
        'quadratic_direction_S61': {f'u{i+1}': str(value)
                                    for i, value in enumerate(coordinates) if value},
        'leading_equation_matrix_constant_then_16_weight_one': string_matrix(matrix),
        'incidence_coefficients_in_old_weight_one_order': [str(x) for x in incidence],
        'rank': 1,
    }


def verify():
    model = S61()
    weights = model.reduced['parameter_weights']
    max_direction_degree = max(sum(weights[i] == 2 for i in monomial)
        for entry in model.data['eliminated_raw_variables']
        for monomial, numerator, denominator in entry['polynomial'])
    assert max_direction_degree == 1
    bianchi = json.loads((model.source/'symplectic_bianchi_identities.json').read_text())
    identities = bianchi['generator_redundancy_identities']
    assert len(identities) == 20
    assert all(weights[parameter] == 1 for identity in identities
        for old, parameter, numerator, denominator in identity['linear_coefficients_on_weight4_old_generators'])
    basis = old_beta_basis(model)
    representative = flag_check(model, s.Matrix([0, -1, 0, 0]),
                                s.Matrix([0, 0, 0, -1]), basis)
    assert representative['quadratic_direction_S61'] == {'u38': '1'}
    nonzero = [(i, j, value) for i, row in enumerate(
        representative['leading_equation_matrix_constant_then_16_weight_one'])
        for j, value in enumerate(row) if value != '0']
    assert nonzero == [(84, 1, '2')]
    additional = flag_check(model, s.Matrix([1, 2, 3, 5]),
                            s.Matrix([0, -16, 1, 7]), basis)
    return {
        'id': 'S61_INFINITY_FLAGS_20261031', 'all_passed': True, 'field': 'QQ',
        'model': 'Degree-two homogenizations of the134 generators in45 fixed-first-jet coordinates, BEFORE saturation by the homogenizing variable.',
        'eliminated_reconstruction_max_degree_in_weight_two': max_direction_degree,
        'Bianchi20_coefficients_only_constant_or_weight_one': True,
        'weight_one_order': [f'u{i+1}' for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 1],
        'representative': representative,
        'representative_only_nonzero_leading_equation': 'equation84 (zero-based) = 2*u4',
        'additional_flag_check_all_first_jets': additional,
        'incidence_injectivity': incidence_injectivity(),
        'claim': 'The unsaturated infinity locus contains the closed isotropic-flag incidence T_beta(v,w)=beta(v)(v,w)=0. It has pure dimension3 for every nonzero primitive first jet; for beta=0 the whole4-dimensional flag variety lies at infinity.',
        'proof_note': 'research/notes/S61_INFINITY_FLAGS_20261031.md',
        'limitations': ['This does not prove positive-dimensional affine fibers.',
                        'The true projective closure requires homogenizing-variable saturation; these infinity points may disappear.',
                        'No upper bound on total generic-fiber dimension and no unknown-component exclusion.'],
        'input_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in (Path(__file__), ROOT/'code/s61_model.py', ROOT/'code/weight1_structure.py',
                                      model.source/'symplectic_slice_reduced_equations.json',
                                      model.source/'symplectic_slice_equations.json',
                                      model.source/'symplectic_slice_certificate.json',
                                      model.source/'symplectic_bianchi_identities.json')},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = verify()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    print(json.dumps({'all_passed': True, 'canonical_residual': '2*u4',
                      'incidence_injectivity_minor': '625', 'output': str(args.output)}))


if __name__ == '__main__':
    main()

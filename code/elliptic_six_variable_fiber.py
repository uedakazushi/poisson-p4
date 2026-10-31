#!/usr/bin/env python3
"""Exact Sp4 and scalar elimination for elliptic first-jet fibers on a dense open.

The 30 invariant equalities in six, then five, variables are given as rational
straight-line prescriptions. A certificate run checks exact
rational witnesses and the explicit order-50 parameter redundancy subgroup.
It does NOT determine either elliptic component degree or the full parameter
redundancy. See research/notes/ELLIPTIC_SIX_VARIABLE_FIBER_20261031.md.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path

import sympy as s

from elliptic_weight1_dominance import (elliptic_family, projective_affine,
                                      symplectic_frame, primitive_vector)
from weight1_structure import J, COORDS

PAIRS = list(combinations(range(4), 2))
ROOT = Path(__file__).resolve().parents[1]
FRAME_DEGREES = (3, 5, 7, 9)
SIGNATURE_DEGREES = ([FRAME_DEGREES[i]+FRAME_DEGREES[j] for i, j in PAIRS]
                     + [1+FRAME_DEGREES[d]+FRAME_DEGREES[i]+FRAME_DEGREES[j]
                        for (i, j), d in COORDS])


def lower_tensor(value):
    """T[a][b,c] is the y_a coefficient of beta_bc."""
    tensor = [s.zeros(4) for _ in range(4)]
    for coefficient, ((b, c), a) in zip(value, COORDS, strict=True):
        tensor[a][b, c] = coefficient
        tensor[a][c, b] = -coefficient
    return tensor


def covariant_frame(value):
    """Return degree-(2,2,3;3,5,7,9) covariants M,A,v,R.

    The degree of A is two, M is symmetric, A=J M is Hamiltonian.
    Under the covariant tensor action of S in Sp4, R transforms to S R.
    """
    tensor = lower_tensor(value)
    # J has only four nonzero entries. Keep contractions short and exact.
    nz = [(i, j, J[i, j]) for i in range(4) for j in range(4) if J[i, j]]
    moment = s.Matrix(4, 4, lambda i, a: sum(
        tensor[i][j, k]*tensor[a][b, c]*u*v
        for j, b, u in nz for k, c, v in nz))
    operator = J*moment
    raised = J*moment*J.T
    covector = s.Matrix([sum(tensor[a][b, i]*raised[a, b]
                            for a in range(4) for b in range(4))
                         for i in range(4)])
    vector = J*covector
    columns = [vector]
    for _ in range(3):
        columns.append(operator*columns[-1])
    frame = s.Matrix.hstack(*columns)
    return moment, operator, vector, frame


def invariant_signature(value):
    """6 symplectic Gram entries and 24 tensor entries in the covariant frame.

    This list separates Sp4 orbits wherever det(R) is nonzero. Equalities of
    signatures define the orbit fiber there scheme-theoretically.
    """
    tensor = lower_tensor(value)
    frame = covariant_frame(value)[-1]
    gram = frame.T*J*frame
    transformed = [sum((frame[a, d]*frame.T*tensor[a]*frame
                        for a in range(4)), s.zeros(4)) for d in range(4)]
    return s.Matrix([gram[i, j] for i, j in PAIRS]
                    + [transformed[d][i, j] for (i, j), d in COORDS])


def transform_jet(value, matrix):
    """Pushforward action on a normalized jet, all three indices lowered."""
    tensor = lower_tensor(value)
    inverse = matrix.inv()
    transformed = [sum((inverse[a, d]*inverse.T*tensor[a]*inverse
                        for a in range(4)), s.zeros(4)) for d in range(4)]
    return s.Matrix([transformed[d][i, j] for (i, j), d in COORDS])


@lru_cache(maxsize=2)
def affine_family(kind):
    x, lam, bivector = elliptic_family(kind)
    z = s.symbols('z1:5')
    affine = projective_affine(bivector, x, z)
    return lam, z, affine, tuple(affine.diff(v) for v in z)


def normalized_jet(kind, parameter, point, scalar=1):
    """Rational first jet b_kind(lambda,p,c) before the remaining Sp4 action.

    Domain: lambda*c*a*Delta != 0, where a,Delta are the fixed symplectic
    frame pivots of K=c P_lambda(p). See fixed_jet_equations for the six-variable
    system. Exact numeric arguments keep the default certificates compact.
    """
    lam, z, affine, derivatives = affine_family(kind)
    subst = dict(zip(z, point, strict=True)) | {lam: s.sympify(parameter)}
    constant = scalar*affine.subs(subst, simultaneous=True)
    differential = [scalar*d.subs(subst, simultaneous=True) for d in derivatives]
    frame = symplectic_frame(constant)
    inverse = -constant*frame.T*J
    linear = [frame*sum((inverse[h, j]*differential[h] for h in range(4)),
                        s.zeros(4))*frame.T for j in range(4)]
    return primitive_vector(linear)


def fixed_jet_equations(kind, target, parameter, point, scalar):
    """Exact 30 rational equations in (parameter,point[0:4],scalar).

    Keep this as a straight-line rational prescription; expanding symbolic
    numerator polynomials is unnecessary and can be expensive. Localize at
    lambda*c*a*Delta*det(R_source)*det(R_target). The unique recovered frame
    sending the source jet to target is R_target * R_source^{-1}.
    """
    source = normalized_jet(kind, parameter, point, scalar)
    return invariant_signature(source)-invariant_signature(target)


def scalar_invariants(value):
    """I4, I14, h=I14/I4^3; h scales by c after normalizing c*pi."""
    tensor = lower_tensor(value)
    _, operator, vector, _ = covariant_frame(value)
    moved = operator*vector
    i4 = s.trace(operator*operator)
    i14 = sum(moved[a]*(vector.T*tensor[a]*moved)[0] for a in range(4))
    return i4, i14, i14/i4**3


def fixed_jet_five_variable_equations(kind, target, parameter, point):
    """Return recovered c and the exact 30 rational equations in lambda,p1:4.

    In addition to the six-variable open conditions, localize at I4*I14 of
    target and source at c=1. All30 equalities are retained; this is not a
    claim that five arbitrarily selected invariants separate orbits.
    """
    at_one = normalized_jet(kind, parameter, point, s.Integer(1))
    recovered_scalar = scalar_invariants(target)[2]/scalar_invariants(at_one)[2]
    source_signature = invariant_signature(at_one)
    target_signature = invariant_signature(target)
    equations = s.Matrix([recovered_scalar**(degree//2)*source_signature[i]
                          - target_signature[i]
                          for i, degree in enumerate(SIGNATURE_DEGREES)])
    return recovered_scalar, equations


def family_symmetries(kind):
    """Check Heisenberg generators and reflection for general lambda.

    Cyclic shift C and diagonal D=diag(1,zeta,...,zeta^4) preserve B.
    Reflection H: x_i -> x_-i sends B to -B. The projective group has order50.
    D is checked using monomial weights modulo5, without floating roots.
    """
    x, lam, bivector = elliptic_family(kind)
    cyclic, reflection = s.zeros(5), s.zeros(5)
    for i in range(5):
        cyclic[i, (i+1) % 5] = 1
        reflection[i, (-i) % 5] = 1
    for matrix, character in [(cyclic, 1), (reflection, -1)]:
        inverse = matrix.T
        pushed = matrix*bivector.subs(dict(zip(x, inverse*s.Matrix(x))),
                                     simultaneous=True)*matrix.T
        assert all(s.expand(v) == 0 for v in pushed-character*bivector)
    for i, j in combinations(range(5), 2):
        for monomial, _ in s.Poly(bivector[i, j], *x).terms():
            assert sum(k*monomial[k] for k in range(5)) % 5 == (i+j) % 5
    # All 10 permutation parts (sign,shift) are distinct, and each has five
    # projectively different diagonal phases, hence exactly50 projective maps.
    permutations = {tuple((sign*i+a) % 5 for i in range(5))
                    for sign in [1, -1] for a in range(5)}
    assert len(permutations) == 10
    return {'cyclic_preserves_general_bivector': True,
            'diagonal_weight_preserves_general_bivector': True,
            'reflection_negates_general_bivector': True,
            'projective_subgroup_order': 50,
            'reflection_parameter_action': 'c -> -c; p -> H(p)',
            'full_parameter_redundancy_determined': False}


def serialize(matrix):
    return [[str(v) for v in row] for row in matrix.tolist()]


def certify():
    # Cayley transform of a Hamiltonian matrix: an exact symplectic element.
    symmetric = s.Matrix([[1, 2, 0, 3], [2, -1, 4, 0],
                           [0, 4, 2, 1], [3, 0, 1, -2]])
    hamiltonian = J*symmetric
    transform = (s.eye(4)-hamiltonian).inv()*(s.eye(4)+hamiltonian)
    assert transform*J*transform.T == J
    report = {'id': 'ELLIPTIC_SIX_VARIABLE_FIBER_20261031',
              'generated_utc': datetime.now(timezone.utc).isoformat(),
              'proof_note': 'research/notes/ELLIPTIC_SIX_VARIABLE_FIBER_20261031.md',
              'source': 'https://arxiv.org/pdf/1001.4422v2',
              'input_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                               for path in [Path(__file__).resolve(),
                                            ROOT/'code/elliptic_weight1_dominance.py',
                                            ROOT/'code/weight1_structure.py']},
              'all_passed': True, 'component_degrees_determined': False,
              'generic_model_variables': ['lambda', 'p1', 'p2', 'p3', 'p4', 'c'],
              'signature_equations': 30, 'signature_degrees': SIGNATURE_DEGREES,
              'equations_format': 'nonexpanded rational straight-line prescription',
              'further_reduced_variables': ['lambda', 'p1', 'p2', 'p3', 'p4'],
              'target_open': 'det(R_beta)*I4(beta)*I14(beta) != 0',
              'transformation_witness': serialize(transform), 'families': []}
    for kind in [1, 2]:
        value = normalized_jet(kind, s.Integer(2), [1, 2, 3, 4], s.Integer(1))
        moment, operator, vector, frame = covariant_frame(value)
        assert moment == moment.T
        assert operator*J+J*operator.T == s.zeros(4)
        determinant = s.factor(frame.det())
        assert determinant != 0
        transformed = transform_jet(value, transform)
        m2, a2, v2, r2 = covariant_frame(transformed)
        assert m2 == transform.inv().T*moment*transform.inv()
        assert a2 == transform*operator*transform.inv()
        assert v2 == transform*vector
        assert r2 == transform*frame
        assert invariant_signature(value) == invariant_signature(transformed)
        recovered = r2*frame.inv()
        assert recovered == transform
        assert transform_jet(value, recovered) == transformed
        assert fixed_jet_equations(kind, value, 2, [1, 2, 3, 4], 1) == s.zeros(30, 1)
        i4, i14, h = scalar_invariants(value)
        assert i4 != 0 and i14 != 0
        scaled = normalized_jet(kind, 2, [1, 2, 3, 4], s.Integer(3))
        scaled_i4, scaled_i14, scaled_h = scalar_invariants(scaled)
        assert scaled_i4 == 3**2*i4 and scaled_i14 == 3**7*i14
        assert scaled_h == 3*h
        assert all(degree % 2 == 0 for degree in SIGNATURE_DEGREES)
        assert invariant_signature(scaled) == s.Matrix([
            3**(degree//2)*invariant_signature(value)[i]
            for i, degree in enumerate(SIGNATURE_DEGREES)])
        assert covariant_frame(scaled)[-1].det() == 3**12*determinant
        recovered_c, equations5 = fixed_jet_five_variable_equations(
            kind, scaled, 2, [1, 2, 3, 4])
        assert recovered_c == 3 and equations5 == s.zeros(30, 1)
        signature = invariant_signature(value)
        assert signature == invariant_signature(normalized_jet(kind, 2, [2, 3, 4, 1], 1))
        assert signature == invariant_signature(normalized_jet(kind, 2, [4, 3, 2, 1], -1))
        report['families'].append({
            'kind': kind, 'component_id': 'T39' if kind == 1 else 'T40',
            'point': [1, 2, 3, 4], 'lambda': 2, 'scalar': 1,
            'primitive_jet': serialize(value), 'M': serialize(moment),
            'A': serialize(operator), 'v': serialize(vector),
            'R': serialize(frame), 'det_R': str(determinant),
            'I4': str(i4), 'I14': str(i14), 'h': str(h),
            'scalar_three_exact_recovery_test': True,
            'cyclic_and_reflection_signature_tests': True,
            'covariance_exact_rational_test': True,
            'signature_and_unique_reconstruction_test': True,
            'symmetry_subgroup': family_symmetries(kind)})
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True,
                        help='New JSON certificate path; existing files are never overwritten.')
    args = parser.parse_args()
    report = certify()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print('PASS: exact covariant-frame witnesses, 30-equation five-variable reduction, '
          'order50 subgroup identities; no elliptic degree asserted.')


if __name__ == '__main__':
    main()

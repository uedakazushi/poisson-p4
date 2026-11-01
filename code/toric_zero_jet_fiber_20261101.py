#!/usr/bin/env python3
"""A ten-dimensional actual toric zero-jet boundary and necessary modular bounds.

The upper bound uses the explicitly named historical classification theorems;
this certificate does not independently repeat that classification.
"""
import argparse
import hashlib
import json
from itertools import combinations, combinations_with_replacement
from pathlib import Path

import sympy as s

from s61_model import S61, X
from weight1_structure import J, sp_basis

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    r, u, z, v, t, a, b, c = s.symbols('r u z v t a b c')
    y = s.Matrix([r, u, z, v])
    frame = s.Matrix([[1, t, 0, 0], [-1, t, 0, 0],
                      [0, 0, 1, t], [0, 0, -1, t]])
    ell = s.ones(4, 1) + frame*y
    logarithmic = s.zeros(4)
    logarithmic[0, 1] = logarithmic[2, 3] = 1/(2*t)
    for i in range(2):
        for j in range(2, 4):
            logarithmic[i, j] = (a*(-1)**(i+j)+b)/4
    logarithmic -= logarithmic.T
    assert frame.det() == 4*t*t
    assert s.factor(logarithmic.det()) == (a*b*t*t-1)**2/(16*t**4)
    poisson = (-frame.inv()*s.diag(*ell)*logarithmic.inv()
               *s.diag(*ell)*frame.inv().T).applyfunc(s.cancel)
    # Every entry is regular at t=0, and remains a quadratic in affine y.
    denominator = a*b*t*t-1
    numerator = poisson.applyfunc(lambda f: s.cancel(f*denominator))
    assert all(s.denom(f) == 1 for f in numerator)
    assert all(s.Poly(f, *y).total_degree() <= 2 for f in numerator)
    limit = poisson.applyfunc(lambda f: f.subs(t, 0))
    expected = s.zeros(4)
    expected[0, 1], expected[2, 3] = 1-r*r, 1-z*z
    expected[1, 3] = -a-b*r*z
    expected -= expected.T
    assert limit == expected
    # The constant normalization and a subsequent symplectic shear.
    shift = s.eye(4)
    shift[1, 2] = a+c
    shift[3, 0] = c
    normalized = (shift*limit*shift.T).applyfunc(s.expand)
    g = c*r*r-b*r*z-(a+c)*z*z
    target = s.zeros(4)
    target[0, 1], target[2, 3], target[1, 3] = 1-r*r, 1-z*z, g
    target -= target.T
    assert (normalized-target).applyfunc(s.expand) == s.zeros(4)
    assert normalized.subs(dict.fromkeys(y, 0)) == J
    assert all(s.diff(f, yi).subs(dict.fromkeys(y, 0)) == 0
               for f in normalized for yi in y)
    pairs = list(combinations(range(4), 2))
    for i, j, k in combinations(range(4), 3):
        assert s.expand(sum(target[i, h]*s.diff(target[j, k], y[h])
                            +target[j, h]*s.diff(target[k, i], y[h])
                            +target[k, h]*s.diff(target[i, j], y[h])
                            for h in range(4))) == 0
    pfaffian = target[0, 1]*target[2, 3]-target[0, 2]*target[1, 3]+target[0, 3]*target[1, 2]
    assert s.factor(pfaffian) == (r-1)*(r+1)*(z-1)*(z+1)

    model = S61()
    homogeneous = s.zeros(5)
    for i, j in pairs:
        polynomial = s.Poly(target[i, j], *y)
        value = sum(co*X[0]**(2-sum(mon))*s.prod(X[k+1]**e for k, e in enumerate(mon))
                    for mon, co in polynomial.terms())
        homogeneous[i+1, j+1], homogeneous[j+1, i+1] = value, -value
    normalized_s61 = model.normalize(homogeneous)
    model.verify_point(normalized_s61)
    one = [i for i, weight in enumerate(model.reduced['parameter_weights']) if weight == 1]
    assert all(normalized_s61['parameters'][i] == 0 for i in one)
    assert all(value == 0 for value in normalized_s61['perspective'])

    monomials = [y[i]*y[j] for i, j in combinations_with_replacement(range(4), 2)]
    def coefficients(matrix):
        return s.Matrix([s.Poly(s.expand(matrix[i, j]), *y).coeff_monomial(mon)
                         for i, j in pairs for mon in monomials])
    quadratic = target-J
    sample = {a: -4, b: -2, c: 1}  # g=r^2+2rz+3z^2
    q0 = quadratic.subs(sample)
    orbit = [coefficients(action*q0+q0*action.T-sum(
        ((action*y)[k]*q0.diff(y[k]) for k in range(4)), s.zeros(4)))
        for action in sp_basis()]
    variations = [coefficients(quadratic.diff(parameter)) for parameter in (a, b, c)]
    tangent = s.Matrix.hstack(*orbit, *variations)
    columns = list(tangent.rref()[1])
    rows = list(tangent[:, columns].T.rref()[1])
    determinant = tangent.extract(rows, columns).det()
    assert len(columns) == 10 and determinant != 0

    # A constant symplectic change takes the family to the historical AB core.
    x, yy, zz, ww = s.symbols('x y z w')
    old = s.Matrix([x, yy, zz, ww])
    change = s.Matrix([[s.I/2, 0, s.I/2, 0], [0, -s.I, 0, -s.I],
                       [s.Rational(1, 2), 0, -s.Rational(1, 2), 0], [0, 1, 0, -1]])
    assert change*J*change.T == J
    changed = (change.inv()*quadratic.subs(dict(zip(y, change*old)), simultaneous=True)
               *change.inv().T).applyfunc(s.expand)
    assert changed[0, 1] == changed[2, 3] == x*zz/2
    assert changed[0, 2] == 0
    assert changed[0, 3] == (x*x+zz*zz)/4
    assert changed[1, 2] == -(x*x+zz*zz)/4
    changed_p = s.Matrix([s.Poly(changed[1, 3], x, zz).coeff_monomial(mon)
                          for mon in (x*x, x*zz, zz*zz)])
    assert changed_p.jacobian((a, b, c)).det() != 0

    # The universal factorization/modular-invariance calculation.
    aa = s.Matrix(s.symbols('a0:4'))
    bb = s.Matrix(s.symbols('b0:4'))
    la, lb = (aa.T*y)[0], (bb.T*y)[0]
    trace, quartic = -la*la-lb*lb, la*la*lb*lb
    pairing = (aa.T*J*bb)[0]
    bracket = (s.Matrix([s.diff(trace, yi) for yi in y]).T*J
               *s.Matrix([s.diff(quartic, yi) for yi in y]))[0]
    assert s.expand(bracket+4*pairing*la*lb*(la*la-lb*lb)) == 0
    hamiltonian = J*s.hessian(trace, y)
    assert (hamiltonian**2-4*pairing*J*(aa*bb.T-bb*aa.T)).applyfunc(s.expand) == s.zeros(4)

    sources = [Path(__file__), ROOT/'code/s61_model.py', ROOT/'code/weight1_structure.py',
               ROOT/'legacy/v07/quadratic_classification.tex',
               ROOT/'legacy/v07/proofs/zero_reduced_components_theorem.txt',
               ROOT/'legacy/v07/code/certify_closed_squarezero.py',
               ROOT/'legacy/v06/ambient_families.tex',
               ROOT/'legacy/v06/data/independent_ambient_families_audit.json']
    report = dict(
        id='TORIC_ZERO_JET_FIBER_20261101', all_passed=True,
        field='QQ(a,b,t) and QQ[a,b,c]; one explicitly displayed QQ(i) coordinate change',
        logarithmic_hyperplanes=list(map(str, ell)),
        logarithmic_matrix=[[str(f) for f in row] for row in logarithmic.tolist()],
        affine_frame_determinant=str(frame.det()),
        logarithmic_matrix_determinant=str(s.factor(logarithmic.det())),
        rational_toric_P=[[str(f) for f in row] for row in poisson.tolist()],
        constant_normalization=[[str(f) for f in row] for row in shift.tolist()],
        limit_quadratic=[[str(f) for f in row] for row in quadratic.tolist()],
        arbitrary_quadratic='g(r,z)=c*r^2-b*r*z-(a+c)*z^2',
        S61_parameters=list(map(str, normalized_s61['parameters'])),
        universal_full134_and_reconstruction_verified=True,
        parameter_sample={str(k): str(vv) for k, vv in sample.items()},
        orbit_plus_family_rank=10, rank_minor_rows=rows, rank_minor_columns=columns,
        rank_minor_determinant=str(determinant),
        rank_minor=[[str(f) for f in row] for row in tangent.extract(rows, columns).tolist()],
        historical_identification='The entire reduced nil22 AB component of the compatible quadratic space lies in the actual T01 zero-first-jet fiber.',
        AB_change=[[str(f) for f in row] for row in change.tolist()],
        AB_core_a='1/4', AB_core_b='-1/4', AB_core_product='-1/16',
        AB_free_p=list(map(str, changed_p)),
        necessary_condition='The quadratic trace Hamiltonian matrix is square-zero; only zero, (2,1,1), and (2,2) modular strata can occur.',
        dimension_bounds={'lower': 10, 'upper': 12},
        upper_bound_dependencies=['v07 exhaustive 22-component classification of the compatible quadratic space and its zero-modular fiber',
                                  'v06 global smooth witnesses in both M0 and M1, whose unique ambient components differ from T01'],
        limitations=['No exact reduced zero-first-jet fiber, upper-bound equality, or exhaustive list of its components is proved.',
                     'Complete factorization of the degeneracy quintic is necessary, not sufficient for T01 membership.',
                     'The historical global classification results are cited dependencies, not rerun by this script.',
                     'No new ambient Poisson component and no full S61 generic exhaustion.'],
        proof_note='research/notes/TORIC_ZERO_JET_FIBER_20261101.md',
        source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'zero_jet_dimension_bounds': [10, 12],
                      'entire_AB_component_in_T01': True}))


if __name__ == '__main__':
    main()

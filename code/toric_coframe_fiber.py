#!/usr/bin/env python3
"""Exact coframe reconstruction of transverse toric first-jet fibers.

The map GL(4) -> W16 parametrizes centered, ordered five-hyperplane toric
structures with P(0)=J. Its S5 quotient is the transverse toric locus. This
does not determine the degree of this map or exhaust non-toric S61 fibers.
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path

import sympy as s

from weight1_structure import COORDS, J, constraints, minor_witness

ROOT = Path(__file__).resolve().parents[1]
Y = s.Matrix(s.symbols("y1:5"))
ONES = s.ones(1, 4)


def linear_frame(z):
    """Linear part of (I-z*1^t/5) diag(1+z_i-sum(z)/5)."""
    return s.diag(*z) - sum(z) * s.eye(4) / 5 - z * ONES / 5


def jet(coframe):
    """Return beta=d(P^{-1}) at zero, in weight1_structure coordinates."""
    c = coframe.inv() * linear_frame(coframe * Y) * coframe
    beta = J * c + c.T * J
    return s.Matrix([s.expand(beta[i, j]).coeff(Y[k]) for (i, j), k in COORDS])


def differential(coframe):
    """Differentiate using d(R^-1)=-R^-1(dR)R^-1, entirely over QQ."""
    inverse = coframe.inv()
    h = linear_frame(coframe * Y)
    c = inverse * h * coframe
    columns = []
    for i in range(4):
        for j in range(4):
            e = s.zeros(4)
            e[i, j] = 1
            dc = -inverse * e * c + inverse * linear_frame(e * Y) * coframe + inverse * h * e
            dbeta = J * dc + dc.T * J
            columns.append(s.Matrix([s.expand(dbeta[a, b]).coeff(Y[k])
                                     for (a, b), k in COORDS]))
    return s.Matrix.hstack(*columns)


def affine_bivector(coframe):
    z = coframe * Y
    h0 = 1 - sum(z) / 5
    h = s.diag(*(h0 + x for x in z))
    frame = coframe.inv() * (s.eye(4) - z * ONES / 5) * h * coframe
    return (frame * J * frame.T).applyfunc(s.expand)


def pfaffian(p):
    return s.expand(p[0, 1] * p[2, 3] - p[0, 2] * p[1, 3] + p[0, 3] * p[1, 2])


def permutation_matrix(permutation):
    """R' rows are r_sigma(i)-r_sigma(0), with r_0=0."""
    rows = [s.zeros(1, 4)] + [s.eye(4)[i, :] for i in range(4)]
    return s.Matrix.vstack(*(rows[permutation[i]] - rows[permutation[0]]
                             for i in range(1, 5)))


def universal_check():
    """Check the full polynomial formula before specializing the period matrix."""
    z = Y
    a, b, c, d, e, f = s.symbols("a b c d e f")
    lam = s.Matrix([[0,a,b,c],[-a,0,d,e],[-b,-d,0,f],[-c,-e,-f,0]])
    h0 = 1 - sum(z) / 5
    h = [h0] + [h0 + x for x in z]
    frame = (s.eye(4) - z * ONES / 5) * s.diag(*h[1:])
    p = (frame * lam * frame.T).applyfunc(s.expand)
    assert all(s.Poly(p[i,j],*z).total_degree() <= 3 for i,j in combinations(range(4),2))
    assert s.expand(pfaffian(p) - (a*f-b*e+c*d)*s.prod(h)) == 0
    for i,j,k in combinations(range(4),3):
        jac = sum(p[i,m]*s.diff(p[j,k],z[m]) + p[j,m]*s.diff(p[k,i],z[m])
                  + p[k,m]*s.diff(p[i,j],z[m]) for m in range(4))
        assert s.expand(jac) == 0
    return {"all_parameter_Jacobi": True, "max_affine_degree": 3,
            "pfaffian_identity": "Pf(P_z)=Pf(Lambda)*product(h_0,...,h_4)",
            "h_0": str(h0), "h_i": "h_0+z_i"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    universal = universal_check()
    # Inverse of the known toric Darboux frame at (2,3,5,7,11,13).
    frame = s.Matrix([[1,0,0,0],[0,s.Rational(1,2),0,0],
                      [s.Rational(7,2),s.Rational(-3,2),1,0],
                      [s.Rational(11,28),s.Rational(-5,28),0,s.Rational(1,14)]])
    r = frame.inv()
    value = jet(r)
    derivative = differential(r)
    constraint = s.Matrix.vstack(*constraints())
    assert constraint * value == s.zeros(8,1)
    assert constraint * derivative == s.zeros(8,16)
    witness = minor_witness(derivative)
    assert witness["rank"] == 16
    p = affine_bivector(r)
    assert p.subs(dict(zip(Y,[0]*4))) == J
    linear = p.applyfunc(lambda entry: sum(s.expand(entry).coeff(y,1).subs(
        {v:0 for v in Y if v != y})*y for y in Y))
    beta = -J*linear*J
    assert s.Matrix([s.expand(beta[i,j]).coeff(Y[k]) for (i,j),k in COORDS]) == value
    assert s.expand(linear[0,1]+linear[2,3]) == 0
    # Four adjacent transpositions generate S5. Check both jets and full fields.
    symmetry = []
    for i in range(4):
        permutation = list(range(5))
        permutation[i], permutation[i+1] = permutation[i+1], permutation[i]
        m = permutation_matrix(permutation)
        assert m.det() == -1
        assert jet(m*r) == value
        assert affine_bivector(m*r) == p
        symmetry.append({"permutation":permutation,"matrix":m.tolist()})
    # Independent comparison with the original six-parameter formula.
    from toric_weight1_dominance import toric_primitive_jet
    params, _, old_jet = toric_primitive_jet()
    assert old_jet.subs(dict(zip(params,[2,3,5,7,11,13]))) == value
    report = {
        "id": "TORIC_COFRAME_FIBER_20261031", "all_passed": True, "field": "QQ",
        "universal_checks": universal,
        "coframe_order": "row-major r_11,...,r_44",
        "sample_coframe": [[str(x) for x in row] for row in r.tolist()],
        "sample_first_jet": [str(x) for x in value],
        "differential_witness": witness,
        "S5_generators": [{"permutation":v["permutation"],
                           "matrix":[[str(x) for x in row] for row in v["matrix"]]}
                          for v in symmetry],
        "toric_formula_comparison": True,
        "fiber_equations": "coeff_y(J*adj(R)*H(R*y)*R + transpose(adj(R)*H(R*y)*R)*J) - det(R)*beta = 0; det(R) != 0",
        "independent_equations": 16, "coframe_variables": 16,
        "polynomial_degrees_in_coframe": [4,5],
        "quotient": "free S5 action on GL4 by permuting the ordered five centered hyperplanes",
        "claim": "The transverse toric fiber over beta is the reduced solution set of the displayed coframe equations modulo S5. The coframe map is dominant and generically finite.",
        "limitations": ["The degree of the coframe map and the degree divided by 120 are not computed.",
                        "This parametrizes the toric open locus, not its boundary and not other S61 components.",
                        "No global exhaustion or new component is asserted."],
        "source_sha256": {"code/toric_coframe_fiber.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "proof_note": "research/notes/TORIC_COFRAME_FIBER_20261031.md",
        "reproduction_command": ".venv/bin/python code/toric_coframe_fiber.py --output reports/<new-name>.json",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as out:
        json.dump(report,out,ensure_ascii=False,indent=2)
        out.write("\n")
    print(json.dumps({"all_passed":True,"rank":16,"minor":witness["determinant"],"output":str(args.output)}))


if __name__ == "__main__":
    main()

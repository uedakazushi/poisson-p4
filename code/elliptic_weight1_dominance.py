#!/usr/bin/env python3
"""Exact dominance certificates for the q(5,1) and q(5,2) first-jet maps.

The rational families are those of Ortenzi--Rubtsov--Tagne Pelap,
arXiv:1001.4422v2, equations (16), (18), (19). Equation (21) has an inconsistent
minus sign: the +lambda of (18) is used and Jacobi is checked independently.
No global exhaustion, degree, or exceptional-fiber claim is made here.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path

import sympy as s

from weight1_structure import (COORDS, J, constraints, infinitesimal_action,
                               minor_witness, sp_basis)

ROOT = Path(__file__).resolve().parents[1]
PAIRS = list(combinations(range(4), 2))
SOURCE_URL = "https://arxiv.org/pdf/1001.4422v2"


def elliptic_family(kind):
    """Return (homogeneous coordinates, lambda, homogeneous bivector).

    kind=1 is q(5,1); kind=2 is q(5,2). No floating point coefficients.
    """
    x = s.symbols("x0:5")
    lam = s.Symbol("lam", nonzero=True)
    if kind == 1:
        a1 = -s.Rational(3, 5)*lam**2 + 1/(5*lam**3)
        a2, a3 = -2/lam, 1/lam**2
        b1 = -s.Rational(1, 5)*lam**2 - 3/(5*lam**3)
        b2, b3 = 2, lam
    elif kind == 2:
        a1 = s.Rational(2, 5)*lam**2 + 1/(5*lam**3)
        a2, a3 = lam, -1/lam
        b1 = -s.Rational(1, 5)*lam**2 + 2/(5*lam**3)
        b2, b3 = -1/lam**2, 1
    else:
        raise ValueError("kind must be 1 or 2")
    bivector = s.zeros(5)
    for i in range(5):
        j = (i+1) % 5
        value = (a1*x[i]*x[j] + a2*x[(i+2) % 5]*x[(i+4) % 5]
                 + a3*x[(i+3) % 5]**2)
        bivector[i, j], bivector[j, i] = value, -value
        j = (i+2) % 5
        value = (b1*x[i]*x[j] + b2*x[(i+3) % 5]*x[(i+4) % 5]
                 + b3*x[(i+1) % 5]**2)
        bivector[i, j], bivector[j, i] = value, -value
    return x, lam, bivector


def projective_affine(bivector, x, z):
    """Induced bracket of z_i=x_(i+1)/x_0 in the x_0=1 chart."""
    subst = dict(zip(x, (1, *z), strict=True))
    return s.Matrix(4, 4, lambda i, j: s.expand(
        (bivector[i+1, j+1] - x[i+1]*bivector[0, j+1]
         + x[j+1]*bivector[0, i+1]).subs(subst)))


def symplectic_frame(matrix):
    """Rational F with F matrix F^T=J on a specified pivot chart.

    The fixed reordering puts the coordinate-vertex elliptic values in this
    chart. The denominators are a and delta=Pf(R matrix R^T).
    """
    reorder = s.eye(4)[[0, 3, 1, 2], :]
    changed = reorder*matrix*reorder.T
    a, b, c, d, e, f = [changed[i, j] for i, j in PAIRS]
    delta = a*f-b*e+c*d
    frame = s.Matrix([[1, 0, 0, 0], [0, 1/a, 0, 0],
                      [d/a, -b/a, 1, 0],
                      [e/delta, -c/delta, 0, a/delta]])*reorder
    return frame.applyfunc(s.cancel)


def primitive_vector(matrices):
    """Finite perspective gauge on the linear jet, then lower indices."""
    y = s.Matrix(s.symbols("y1:5"))
    linear = sum((y[k]*matrices[k] for k in range(4)), s.zeros(4))
    tau = linear[0, 1]+linear[2, 3]
    alpha = s.Matrix([s.diff(tau, coordinate)/5 for coordinate in y])
    v = J*alpha
    primitive = linear-2*alpha.dot(y)*J+y*v.T-v*y.T
    beta = -J*primitive*J
    return s.Matrix([s.expand(beta[i, j]).coeff(y[k]) for (i, j), k in COORDS])


def sample_data(kind):
    """Exact family and normalized data at lambda=2, z=0, scalar=1."""
    x, lam, bivector = elliptic_family(kind)
    z = s.symbols("z1:5")
    affine = projective_affine(bivector, x, z)
    point = dict(zip(z, [0]*4, strict=True)) | {lam: s.Integer(2)}
    constant = affine.subs(point)
    derivatives = [affine.diff(coordinate).subs(point) for coordinate in z]
    frame = symplectic_frame(constant)
    assert frame*constant*frame.T == J
    inverse = frame.inv()
    linear = [frame*sum((inverse[h, j]*derivatives[h] for h in range(4)),
                         s.zeros(4))*frame.T for j in range(4)]
    return x, lam, bivector, z, affine, point, constant, derivatives, frame, linear


def frame_derivative(constant, direction):
    entries = s.symbols("k0:6")
    symbolic = s.zeros(4)
    for value, (i, j) in zip(entries, PAIRS, strict=True):
        symbolic[i, j], symbolic[j, i] = value, -value
    frame = symplectic_frame(symbolic)
    point = dict(zip(entries, [constant[i, j] for i, j in PAIRS], strict=True))
    return sum((frame.diff(v).subs(point)*direction[i, j]
                for v, (i, j) in zip(entries, PAIRS, strict=True)), s.zeros(4))


def serialized(matrix):
    return [[str(value) for value in row] for row in matrix.tolist()]


def certify(kind):
    (x, lam, bivector, z, affine, point, constant, derivatives,
     frame, linear) = sample_data(kind)
    # All ten identities over QQ(lambda)[x], not merely at a rational sample.
    for i, j, k in combinations(range(5), 3):
        jacobi = sum(bivector[i, h]*s.diff(bivector[j, k], x[h])
                     + bivector[j, h]*s.diff(bivector[k, i], x[h])
                     + bivector[k, h]*s.diff(bivector[i, j], x[h])
                     for h in range(5))
        assert s.expand(jacobi) == 0
    for i, j, k in combinations(range(4), 3):
        jacobi = sum(affine[i, h]*s.diff(affine[j, k], z[h])
                     + affine[j, h]*s.diff(affine[k, i], z[h])
                     + affine[k, h]*s.diff(affine[i, j], z[h])
                     for h in range(4))
        assert s.expand(jacobi) == 0
    inverse = frame.inv()
    value = primitive_vector(linear)
    horizontal_columns, rational_columns = [], []
    for parameter in [lam, *z, None]:
        # None denotes the c derivative of c*affine at c=1.
        dconstant = affine.diff(parameter).subs(point) if parameter else constant
        dderivatives = [(affine.diff(v).diff(parameter).subs(point)
                         if parameter else derivatives[i])
                        for i, v in enumerate(z)]
        horizontal = frame*dconstant*frame.T*J/2
        actual = frame_derivative(constant, dconstant)*inverse
        for h in (horizontal, actual):
            assert h*J+J*h.T+frame*dconstant*frame.T == s.zeros(4)
        difference = actual-horizontal
        assert difference*J+J*difference.T == s.zeros(4)
        columns = []
        for h in (horizontal, actual):
            dlinear = [h*linear[j]+linear[j]*h.T
                       - sum((h[k, j]*linear[k] for k in range(4)), s.zeros(4))
                       + frame*sum((inverse[k, j]*dderivatives[k]
                                    for k in range(4)), s.zeros(4))*frame.T
                       for j in range(4)]
            columns.append(primitive_vector(dlinear))
        assert columns[1]-columns[0] == -infinitesimal_action(value, difference)
        horizontal_columns.append(columns[0])
        rational_columns.append(columns[1])
    action = [infinitesimal_action(value, generator) for generator in sp_basis()]
    horizontal_differential = s.Matrix.hstack(*horizontal_columns, *action)
    rational_differential = s.Matrix.hstack(*rational_columns, *action)
    constraint = s.Matrix.vstack(*constraints())
    assert constraint*value == s.zeros(8, 1)
    assert constraint*horizontal_differential == s.zeros(8, 16)
    assert constraint*rational_differential == s.zeros(8, 16)
    horizontal_witness = minor_witness(horizontal_differential)
    rational_witness = minor_witness(rational_differential)
    assert horizontal_witness["rank"] == rational_witness["rank"] == 16
    assert horizontal_witness["determinant"] == rational_witness["determinant"]
    expected = (s.Rational(142978515625, 512) if kind == 1 else
                -s.Rational(142978515625, 33554432))
    assert s.Rational(rational_witness["determinant"]) == expected
    # A direct countercheck for the inconsistent printed minus in equation (21).
    wrong_sign_failure = None
    if kind == 1:
        wrong = bivector.copy()
        for i in range(5):
            j = (i+2) % 5
            wrong[i, j] -= 2*lam*x[(i+1) % 5]**2
            wrong[j, i] = -wrong[i, j]
        wrong_sign_failure = s.factor(sum(
            wrong[0, h]*s.diff(wrong[1, 2], x[h])
            + wrong[1, h]*s.diff(wrong[2, 0], x[h])
            + wrong[2, h]*s.diff(wrong[0, 1], x[h]) for h in range(5)))
        assert wrong_sign_failure != 0
    return {
        "component_id": "T39" if kind == 1 else "T40",
        "family": f"q(5,{kind})", "all_passed": True,
        "source": {"url": SOURCE_URL, "equations": ["16", "18" if kind == 1 else "19"],
                   "page_one_based": 10},
        "generic_homogeneous_Jacobi_identities_checked": 10,
        "generic_affine_Jacobi_identities_checked": 4,
        "homogeneous_bivector_upper_triangle": [str(bivector[i, j])
                                                  for i, j in combinations(range(5), 2)],
        "sample_lambda": "2", "sample_projective_point": [1, 0, 0, 0, 0],
        "sample_scalar": "1", "sample_constant": serialized(constant),
        "sample_frame": serialized(frame), "sample_primitive_jet": [str(v) for v in value],
        "column_order": ["lambda", "p1", "p2", "p3", "p4", "scalar_c"]
                        + [f"sp4_{i}" for i in range(10)],
        "rational_frame_differential": serialized(rational_differential),
        "rational_frame_minor_witness": rational_witness,
        "horizontal_frame_minor_witness": horizontal_witness,
        "horizontal_and_rational_frame_differentials_differ_only_by_Sp4_columns": True,
        "eq21_wrong_minus_sign_Jacobi_012": (str(wrong_sign_failure)
                                              if wrong_sign_failure is not None else None),
        "claim": "Dominant primitive first-jet projection to W16; generically finite "
                 "using the existing 25-dimensional projective component and S61 dimension formula.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT/"reports/ELLIPTIC_WEIGHT1_DOMINANCE_20261031.json")
    args = parser.parse_args()
    families = [certify(1), certify(2)]
    report = {
        "id": "ELLIPTIC_WEIGHT1_DOMINANCE_20261031",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "all_passed": True, "field": "QQ", "families": families,
        "domain": "open subset of (lambda != 0) x A4_basepoint x Gm_scalar x Sp4",
        "target": "W16: primitive closed linear two-forms for standard J",
        "coordinate_order": [{"pair_zero_based": list(pair), "linear_variable_zero_based": k}
                             for pair, k in COORDS],
        "limitations": ["No degrees or branch divisors computed.",
                        "No assertion that all jets lift or that lifts are unique.",
                        "Unknown components and exceptional fibers remain unclassified.",
                        "Component labels and dimensions use the existing catalog and cited literature."],
        "source_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in (Path(__file__), ROOT/"code/weight1_structure.py")},
        "proof_note": "research/notes/ELLIPTIC_WEIGHT1_20261031.md",
        "reproduction_command": ".venv/bin/python code/elliptic_weight1_dominance.py --output <new-report.json>",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"all_passed": True, "ranks": [16, 16],
                      "determinants": [v["rational_frame_minor_witness"]["determinant"] for v in families],
                      "output": str(args.output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()

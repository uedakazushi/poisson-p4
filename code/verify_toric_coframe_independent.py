#!/usr/bin/env python3
"""Independent symbolic checks of the toric coframe construction.

Checks the actual projective coordinate change and universal S5 invariance,
with six free toric coefficients. It does not compute a generic degree.
"""
import argparse
import hashlib
from itertools import permutations
import json
from pathlib import Path

import sympy as s

from toric_coframe_fiber import ONES, Y, permutation_matrix

ROOT = Path(__file__).resolve().parents[1]


def verify():
    z = Y
    a, b, c, d, e, f = s.symbols("a b c d e f")
    lam = s.Matrix([[0, a, b, c], [-a, 0, d, e],
                    [-b, -d, 0, f], [-c, -e, -f, 0]])

    def frame(y):
        h0 = 1-sum(y)/5
        return (s.eye(4)-y*ONES/5)*s.diag(*(h0+t for t in y))

    bivector = (frame(z)*lam*frame(z).T).applyfunc(s.expand)
    for i in range(4):
        permutation = list(range(5))
        permutation[i], permutation[i+1] = permutation[i+1], permutation[i]
        matrix = permutation_matrix(permutation)
        transformed = matrix.inv()*frame(matrix*z)*matrix
        assert (transformed*lam*transformed.T-bivector).applyfunc(s.expand) == s.zeros(4)

    matrices = [permutation_matrix(p) for p in permutations(range(5))]
    assert len({tuple(matrix) for matrix in matrices}) == 120
    assert all(matrix.det() in (-1, 1) for matrix in matrices)
    assert sum(matrix == s.eye(4) for matrix in matrices) == 1

    # Independent differentiation of the inverse projective coordinate change.
    w = s.Matrix(s.symbols("w1:5"))
    h0 = 5/(1+sum(w))
    inverse_coordinates = h0*(w-s.ones(4, 1))
    derivative = inverse_coordinates.jacobian(w)
    assert (derivative-h0*(s.eye(4)-inverse_coordinates*ONES/5)).applyfunc(s.cancel) == s.zeros(4)
    pushforward = derivative*s.diag(*w)*lam*s.diag(*w)*derivative.T
    substituted = bivector.subs(dict(zip(z, inverse_coordinates)), simultaneous=True)
    assert (pushforward-substituted).applyfunc(s.cancel) == s.zeros(4)
    return {
        "all_passed": True, "field": "QQ(a,b,c,d,e,f)",
        "universal_S5_generator_identities": 4,
        "distinct_integral_permutation_matrices": 120,
        "only_identity_permutation_matrix_is_identity": True,
        "independent_projective_change_checked": True,
        "projective_change": "w_i=h_i/h_0; z_i=5*(w_i-1)/(1+sum(w))",
        "scope": "Symbolic formulas, projective extension, and permutation action; "
                 "no scheme-theoretic quotient or generic-degree computation.",
        "source_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in (Path(__file__), ROOT/"code/toric_coframe_fiber.py")},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = verify()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"all_passed": True, "output": str(args.output)}))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Independent finite exact checks of the S61 adapter's coordinate conventions.

These are regression checks of the implementation, not component-exhaustion
certificates. The generic six-parameter symplectic frame is checked symbolically;
the group-action and gauge tests use explicit rational inputs.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import sympy as s

from s61_model import (BASIS, J, PAIRS, ROOT, S61, X,
                       bivector_from_coefficients, canonical, coefficients,
                       jacobi_projective, pushforward, quotient_coefficients,
                       toric_family)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "reports/S61_ADAPTER_20261030.json")
    args = parser.parse_args()
    checks = []

    def check(name, condition, **details):
        assert condition, name
        checks.append({"name": name, "passed": True, **details})

    def same_matrix(a, b):
        return a.shape == b.shape and all(canonical(x - y) == 0
                                        for x, y in zip(a, b, strict=True))

    def same_vector(a, b):
        return len(a) == len(b) and all(canonical(x - y) == 0
                                      for x, y in zip(a, b, strict=True))

    def rejects(name, function, exception=ValueError):
        try:
            function()
        except exception:
            check(name, True, exception=exception.__name__)
        else:
            raise AssertionError(name + ": invalid input was accepted")

    sys.path.insert(0, str(ROOT / "legacy/v02/code"))
    import poisson_p4_components_check as old
    basis = [(pair, s.Poly(mon, *X).monoms()[0]) for pair, mon in BASIS]
    check("all_150_basis_entries_equal_archived_ordering",
          basis == old.BIVECTOR_BASIS, entries=150)
    values = list(range(-75, 75))
    check("all_150_coefficients_roundtrip",
          coefficients(bivector_from_coefficients(values)) == values,
          entries=150)

    model = S61()
    a, b, c, d, e, f = s.symbols("a b c d e f")
    delta = a*f - b*e + c*d
    lam = s.Matrix([[0,a,b,c],[-a,0,d,e],[-b,-d,0,f],[-c,-e,-f,0]])
    # Independent evaluation of the exact frame specified by the public family.
    symbolic, frame, actual_delta = toric_family(a,b,c,d,e,f)
    check("six_parameter_toric_frame_congruence",
          same_matrix(frame[1:5,1:5] * lam * frame[1:5,1:5].T, J[1:5,1:5]),
          field="QQ(a,b,c,d,e,f)", open_condition="a*(a*f-b*e+c*d) != 0")
    check("six_parameter_frame_determinant",
          canonical(frame.det() - 1/delta) == 0 and actual_delta == delta)
    check("six_parameter_value_normalized_at_e0",
          same_matrix(symbolic.subs(dict(zip(X,[1,0,0,0,0]))), J))

    point, _, _ = toric_family(2,3,5,7,11,13)
    normalized = model.normalize(point)
    verification = model.verify_point(normalized)
    check("rational_toric_point_reconstruction_and_equations", True,
          input_parameters=[2,3,5,7,11,13], **verification)
    normal_bivector = bivector_from_coefficients(normalized["coefficients"])
    check("normalization_idempotent",
          same_vector(model.normalize(normal_bivector)["coefficients"],
                      normalized["coefficients"]))

    # Every one of the 25 standard vertical generators, with symbolic weights.
    gauge_parameters = s.symbols("g0:25")
    field = s.Matrix(5,5,gauge_parameters) * s.Matrix(X)
    vertical = s.Matrix(X)*field.T - field*s.Matrix(X).T
    check("generic_E_wedge_EndV_gauge_invariance",
          same_vector(quotient_coefficients(point + vertical),
                      quotient_coefficients(point)), gauge_parameters=25)

    perspective = s.eye(5)
    perspective[0,1:5] = s.Matrix([[2,-3,5,-7]])
    shifted = model.normalize(pushforward(normal_bivector,perspective))
    check("finite_perspective_orbit_normalization",
          same_vector(shifted["parameters"],normalized["parameters"]) and
          same_vector(shifted["perspective"],[-2,3,-5,7]),
          perspective=[2,-3,5,-7])

    matrix = s.eye(5)
    matrix[0,2]=2
    matrix[3,1]=-3
    other = s.diag(1,2,3,5,7)
    check("pushforward_composition_order",
          same_matrix(pushforward(pushforward(point,other),matrix),
                      pushforward(point,matrix*other)))
    check("pushforward_inverse",
          same_matrix(pushforward(pushforward(point,matrix),matrix.inv()),point))

    wrong_value = X[0]**2 * J * 2
    rejects("reject_wrong_evaluation",lambda:model.normalize(wrong_value))
    rejects("reject_toric_zero_first_pivot",lambda:toric_family(0,1,0,0,1,1))
    rejects("reject_toric_zero_pfaffian",lambda:toric_family(1,1,0,0,1,1))
    rejects("reject_wrong_parameter_count",lambda:model.reconstruct([0]*60))
    wrong_degree = X[0]**2 * J
    wrong_degree[1,2] += X[0]**3
    wrong_degree[2,1] -= X[0]**3
    rejects("reject_cubic_contamination",lambda:model.normalize(wrong_degree))
    wrong_constant = X[0]**2 * J
    wrong_constant[1,2] += 1
    wrong_constant[2,1] -= 1
    rejects("reject_constant_contamination",lambda:coefficients(wrong_constant))
    nonskew = X[0]**2 * J
    nonskew[2,1] = 0
    rejects("reject_non_skew_matrix",lambda:coefficients(nonskew))
    rejects("reject_wrong_matrix_size",lambda:coefficients(s.zeros(4)))

    nonpoisson = X[0]**2*J
    nonpoisson[3,4] += X[1]**2
    nonpoisson[4,3] -= X[1]**2
    check("negative_control_has_nonzero_projective_jacobi",
          any(v != 0 for v in jacobi_projective(nonpoisson)))
    rejects("reject_nonpoisson_point_verification",
            lambda:model.verify_point(model.normalize(nonpoisson)),AssertionError)

    # Exercise the optimized common-denominator checker over QQ(t), including
    # a genuine pole and an independently nonzero symbolic negative control.
    t = s.symbols("t")
    rank_one_parameters = [s.S.Zero]*61
    rank_one_parameters[33] = t/(1+t)
    model.verify_equations(rank_one_parameters)
    check("optimized_symbolic_rank_one_equations",
          all(model.evaluate(eq,rank_one_parameters) == 0
              for eq in model.reduced["equations"]),
          parameter="u34=t/(1+t), all other ui=0", field="QQ(t)",
          open_condition="t != -1", equations_compared=134)
    symbolic_nonpoisson = X[0]**2*J
    symbolic_nonpoisson[3,4] += t*X[1]**2
    symbolic_nonpoisson[4,3] -= t*X[1]**2
    negative_parameters = model.normalize(symbolic_nonpoisson)["parameters"]
    negative_values = [model.evaluate(eq,negative_parameters)
                       for eq in model.reduced["equations"]]
    check("scalar_symbolic_negative_control_nonzero",
          any(value != 0 for value in negative_values),
          nonzero_equations_zero_based=[i for i,v in enumerate(negative_values) if v != 0])
    rejects("optimized_symbolic_negative_control_rejected",
            lambda:model.verify_equations(negative_parameters),AssertionError)

    report = {
        "id":"S61_ADAPTER_20261030",
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "all_passed":all(x["passed"] for x in checks),
        "check_groups":len(checks),
        "checks":checks,
        "implementation_sha256":hashlib.sha256((ROOT/"code/s61_model.py").read_bytes()).hexdigest(),
        "verifier_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope":"Exact coordinate and finite group-action regression checks; no prime decomposition, family dominance, or global exhaustion claim.",
        "reproduction_command":"PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/verify_s61_adapter.py",
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"all_passed":report["all_passed"],
                      "check_groups":len(checks),"report":str(args.output)}))


if __name__ == "__main__":
    main()

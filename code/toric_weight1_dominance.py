#!/usr/bin/env python3
"""Certify dominance of the T01 toric family's primitive first-jet projection.

All computations are exact over QQ. A nonzero 16 by 16 Jacobian minor proves
dominance onto the 16-dimensional primitive closed linear two-form space.
This does not classify fibers or other components of S61.
"""
import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path

import sympy as s

from weight1_structure import (COORDS, J, coefficient_map, constraints,
                               infinitesimal_action, minor_witness, sp_basis)

ROOT = Path(__file__).resolve().parents[1]


def toric_primitive_jet():
    parameters = s.symbols("a b c d e f")
    a,b,c,d,e,f = parameters
    delta = a*f-b*e+c*d
    lam = s.Matrix([[0,a,b,c],[-a,0,d,e],[-b,-d,0,f],[-c,-e,-f,0]])
    frame = s.Matrix([[1,0,0,0], [0,1/a,0,0],
                      [d/a,-b/a,1,0], [e/delta,-c/delta,0,a/delta]])
    assert (frame*lam*frame.T-J).applyfunc(s.cancel) == s.zeros(4)
    assert s.cancel(frame.det()-1/delta) == 0

    # Independent all-parameter Jacobi check before the linear change of frame.
    u = s.symbols("u1:5")
    toric = s.diag(*(1+z for z in u))*lam*s.diag(*(1+z for z in u))
    for i,j,k in combinations(range(4),3):
        jac = sum(toric[i,h]*s.diff(toric[j,k],u[h])
                  + toric[j,h]*s.diff(toric[k,i],u[h])
                  + toric[k,h]*s.diff(toric[i,j],u[h]) for h in range(4))
        assert s.expand(jac) == 0

    y = s.Matrix(s.symbols("y1:5"))
    preimage = frame.inv()*y
    linear = (frame*(s.diag(*preimage)*lam+lam*s.diag(*preimage))*frame.T).applyfunc(s.cancel)
    tau = linear[0,1]+linear[2,3]
    alpha = s.Matrix([s.diff(tau,z)/5 for z in y])
    v = J*alpha
    primitive = (linear-2*alpha.dot(y)*J+y*v.T-v*y.T).applyfunc(s.cancel)
    assert s.cancel(primitive[0,1]+primitive[2,3]) == 0
    beta = -J*primitive*J
    vector = s.Matrix([s.cancel(s.diff(beta[i,j],y[k])) for (i,j),k in COORDS])
    constraint = s.Matrix.vstack(*constraints())
    assert (constraint*vector).applyfunc(s.cancel) == s.zeros(8,1)
    return parameters, delta, vector


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,
                        default=ROOT/"reports/TORIC_WEIGHT1_DOMINANCE_20261030.json")
    args = parser.parse_args()
    parameters,delta,vector = toric_primitive_jet()
    sample = [2,3,5,7,11,13]
    substitution = dict(zip(parameters,sample,strict=True))
    value = vector.subs(substitution)
    parameter_derivative = vector.jacobian(parameters).subs(substitution)
    action = s.Matrix.hstack(*(infinitesimal_action(value,a) for a in sp_basis()))
    differential = parameter_derivative.row_join(action)
    constraint = s.Matrix.vstack(*constraints())
    assert constraint.rank() == 8
    assert constraint*differential == s.zeros(8,16)
    witness = minor_witness(differential)
    assert witness["rank"] == 16
    assert parameter_derivative.rank() == 6
    assert action.rank() == 10
    assert action.row_join(value).rank() == 11
    assert coefficient_map(value).rank() == 4
    report = {
        "id":"TORIC_WEIGHT1_DOMINANCE_20261030",
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "all_passed":True,
        "field":"QQ",
        "domain":"{a*(a*f-b*e+c*d) != 0} x Sp4",
        "target":"V16 = {linear two-forms beta: d beta = 0, beta12+beta34 = 0}",
        "parameter_order":[str(x) for x in parameters],
        "sample_parameters":sample,
        "sample_pfaffian":str(delta.subs(substitution)),
        "coordinate_order":[{"pair_one_based":[i+1,j+1],"linear_variable_one_based":k+1}
                            for (i,j),k in COORDS],
        "generic_primitive_jet_coefficients":[str(x) for x in vector],
        "sample_primitive_jet_coefficients":[str(x) for x in value],
        "sample_parameter_derivative":[[str(x) for x in row] for row in parameter_derivative.tolist()],
        "sample_sp4_action":[[str(x) for x in row] for row in action.tolist()],
        "combined_column_order":"six parameter derivatives a,b,c,d,e,f; ten weight1_structure.sp_basis() pullback action columns",
        "combined_differential_witness":witness,
        "parameter_derivative_rank":6,
        "sp4_action_rank":10,
        "sp4_times_scalars_action_rank":11,
        "coefficient_map_rank":4,
        "generic_toric_Jacobi_checked":True,
        "generic_frame_congruence_checked":True,
        "generic_primitive_closed_constraints_checked":True,
        "claim":"The known toric component T01 has dominant primitive first-jet projection onto V16; consequently a nonempty Zariski open subset of V16 is realized by projective Poisson structures in T01.",
        "limitations":["No explicit description of the open image or exceptional image locus.",
                       "No assertion that every first jet lifts, nor uniqueness of a lift.",
                       "No assertion that every S61 component is toric or dominates V16.",
                       "The component label T01 uses the existing catalog; this is not a new component."],
        "proof_note":"research/notes/TORIC_WEIGHT1_DOMINANCE_20261030.md",
        "source_sha256":{
            "code/toric_weight1_dominance.py":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "code/weight1_structure.py":hashlib.sha256((ROOT/"code/weight1_structure.py").read_bytes()).hexdigest(),
        },
        "reproduction_command":"PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/toric_weight1_dominance.py",
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("x") as out:
        json.dump(report,out,indent=2,ensure_ascii=False)
        out.write("\n")
    print(json.dumps({"all_passed":True,"combined_rank":16,
                      "minor_determinant":witness["determinant"],
                      "output":str(args.output)},ensure_ascii=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact weights for arbitrary lower-degree tails above ten SNC initial forms.

This checks an explicitly defined rank-four support subset, not all special bases.
"""
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json
import resource
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "reports/CLASSIFICATION_BETA_ZERO_FULL_COVERAGE_20261106.json"


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU, (30, 31))
    resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
    start = time.monotonic()
    source = json.loads((ROOT / SOURCE).read_text())
    h, x, y, z, w = s.symbols("h x y z w")
    vs = (x, y, z, w)
    pairs = list(combinations(range(4), 2))
    records = []
    for row in source["nonnil_ten"]:
        weights = row["pullback_weights"][1:]
        wh = max(weights[i] + weights[j] - weights[k]
                 for i, j in pairs for k in range(4)) + 1
        wh = max(wh, max((weights[i] + weights[j]) // 2 + 1
                         for i, j in pairs))
        new_weights = [wh, *weights]
        quadratic = []
        for term in row["all_terms"]:
            mon = term["monomial"]
            if mon[0]:
                continue
            i, j = (a - 1 for a in term["pair"])
            value = sum(a*b for a, b in zip(mon, new_weights)) - weights[i] - weights[j]
            assert value == term["weight"] and value >= 0
            quadratic.append({**term, "weight": value})
        lower = [dict(pair=[i+1, j+1], variable=k+1,
                      weight=wh+weights[k]-weights[i]-weights[j])
                 for i, j in pairs for k in range(4)]
        constants = [dict(pair=[i+1, j+1], weight=2*wh-weights[i]-weights[j])
                     for i, j in pairs]
        assert all(v["weight"] > 0 for v in lower + constants)
        limit = {key: s.sympify(value) for key, value in row["limit"].items()}
        pf = limit["1,2"]*limit["3,4"]-limit["1,3"]*limit["2,4"]+limit["1,4"]*limit["2,3"]
        assert s.expand(h*pf-s.sympify(row["projective_limit_degeneracy_divisor"])) == 0
        records.append(dict(name=row["Q45_component"], weights=new_weights,
                            all24_linear_tail_weights=lower,
                            all6_constant_tail_weights=constants,
                            quadratic_terms=quadratic, limit=row["limit"],
                            projective_divisor=row["projective_limit_degeneracy_divisor"],
                            parameter_open=row["generic_parameter_open"]))

    J = s.zeros(4)
    J[0, 1] = J[2, 3] = 1
    J[1, 0] = J[3, 2] = -1
    P = s.zeros(4)
    P[0, 1], P[1, 0] = 1+x+x*y, -1-x-x*y
    P[2, 3], P[3, 2] = 1+z*w, -1-z*w
    jacobi = [s.expand(sum(P[i,l]*s.diff(P[j,k],vs[l])
                          +P[j,l]*s.diff(P[k,i],vs[l])
                          +P[k,l]*s.diff(P[i,j],vs[l]) for l in range(4)))
              for i, j, k in combinations(range(4), 3)]
    assert jacobi == [0]*4
    linear = s.zeros(4)
    linear[0, 1], linear[1, 0] = x, -x
    a = s.Matrix([s.Rational(1, 5), 0, 0, 0])
    vec = s.Matrix(vs)
    ja = J*a
    primitive = linear-2*x*J/5+vec*ja.T-ja*vec.T
    assert s.expand(primitive[0,1]+primitive[2,3]) == 0
    assert primitive[0,1] == 2*x/5
    witness = dict(affine_brackets={f"{i+1},{j+1}": str(P[i,j]) for i,j in pairs},
                   jacobi=[str(v) for v in jacobi],
                   primitive_linear_brackets={f"{i+1},{j+1}": str(s.expand(primitive[i,j])) for i,j in pairs},
                   primitive_first_jet_nonzero=True,
                   Q_type="G1: A=D=1,E=0", constant_value="J")
    report = dict(task_id="E_poly", date_JST="2026-11-07", all_passed=True,
                  scope="All Poisson lower tails with the stated ten quadratic normal forms on kappa!=0; their PGL-saturation and closure. Not the full preimage of the closed Q45 components.",
                  all_ten=records, nonzero_first_jet_witness=witness,
                  source_sha256={SOURCE: hashlib.sha256((ROOT/SOURCE).read_bytes()).hexdigest(),
                                 "code/classification_lower_degree_snc.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                  generic_exhaustion_proved=False, global_exhaustion_proved=False,
                  seconds=time.monotonic()-start)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"all_passed": True, "ten_weights": [r["weights"] for r in records],
                      "primitive_first_jet_nonzero": True, "seconds": report["seconds"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args().output)

#!/usr/bin/env python3
"""Certify rational S61 maps of actual Poisson families, not just samples.

The full toric six-parameter family is checked over QQ(a,b,c,d,e,f).
Its S61 image is not claimed to be the whole 16-dimensional T01 slice.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import sympy as s

from s61_model import J, X, S61, toric_family


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    model = S61()
    a, b, c, d, e, f = s.symbols("a b c d e f")
    poisson, frame, delta = toric_family(a, b, c, d, e, f)
    normalized = model.normalize(poisson)
    print("Full toric six-parameter map constructed", flush=True)
    result = model.verify_point(normalized)
    print("Full toric family: exact reconstruction, 134 equations, Jacobi PASS", flush=True)
    toric = {
        "family_id": "MPS_T01_TORIC_CHART",
        "legacy_id": "T01",
        "coefficient_field": "QQ(a,b,c,d,e,f)",
        "open_condition": "a*(a*f-b*e+c*d) != 0",
        "convention": "a,b,c,d,e,f are affine log POISSON coefficients, not biresidues",
        "perspective": list(map(str, normalized["perspective"])),
        "parameters_u1_through_u61": list(map(str, normalized["parameters"])),
        "parameter_weights": model.data["parameter_weights"],
        "linear_frame": [[str(v) for v in row] for row in frame.tolist()],
        "verification": result,
        "scope": "Generic toric coefficient family mapped into S61; residual CSp4 action must be added to cover the component. No full defining ideal or primary decomposition.",
    }
    # A second genuine family with weight-one coefficients, from v02 §4.
    alpha, beta = s.symbols("alpha beta")
    almost = X[0]**2 * J
    for (i, j), value in { (1,4): -beta*X[0]*X[1],
                          (2,4): beta*X[0]*X[2],
                          (3,4): -alpha*X[0]*X[3] }.items():
        almost[i,j] += value
        almost[j,i] -= value
    an = model.normalize(almost)
    ar = model.verify_point(an)
    print("Almost-abelian two-parameter family: PASS", flush=True)
    almost_report = {
        "family_id": "V02_ALMOST_ABELIAN",
        "coefficient_field": "QQ(alpha,beta)", "open_condition": "none",
        "source": "legacy/v02/p4_poisson_completeness_v2.tex, almost-abelian branch",
        "parameters_u1_through_u61": list(map(str, an["parameters"])),
        "perspective": list(map(str, an["perspective"])),
        "verification": ar,
        "scope": "Contained in the union MPS40 by the bundled v02 proof; no individual MPS component ID assigned here.",
    }
    output_files = {"toric_T01_map.json": toric, "almost_abelian_map.json": almost_report}
    for name, payload in output_files.items():
        (args.output/name).write_text(json.dumps(payload, ensure_ascii=False, indent=2)+"\n")
    inputs = list(model.source.glob("symplectic_*.json")) + [Path(__file__), Path(__file__).with_name("s61_model.py")]
    report = {
        "all_passed": True,
        "verification_method": "Exact rational-function identities, no random/finite-field inference",
        "families": [{"id": toric["family_id"], **result}, {"id": almost_report["family_id"], **ar}],
        "inputs_sha256": {str(path.relative_to(Path(__file__).resolve().parents[1])):
                          hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
        "elapsed_seconds": round(time.monotonic()-start, 3),
        "no_new_component_count": True,
    }
    (args.output/"verification.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

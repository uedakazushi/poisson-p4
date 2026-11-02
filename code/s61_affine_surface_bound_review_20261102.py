"""Independent intrinsic trace construction of the marked binary factor."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import sympy as sym

from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import J

ROOT = Path(__file__).resolve().parents[1]
STEM = "S61_AFFINE_SURFACE_BOUND_REVIEW_20261102"
NOTE = f"research/notes/{STEM}.md"
PRIMARY = "reports/S61_AFFINE_SURFACE_BOUND_20261102.json"


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def verify():
    start = time.monotonic()
    report = json.loads((ROOT/PRIMARY).read_text())
    assert report["all_passed"] and report["generic_affine_dimension_upper"] == 2
    paths = [str(Path(__file__).relative_to(ROOT)),NOTE,PRIMARY]
    dependencies = [PRIMARY, "reports/S61_AFFINE_BOUNDARY_DIMENSION_20261102.json",
                    "reports/S61_BINARY_H11_RANK_DROP_RESIDUAL_FINITE_20261102.json"]
    for source in dependencies:
        data = json.loads((ROOT/source).read_text())
        assert data["all_passed"]
        for p,h in data["input_sha256"].items():
            assert sha(p) == h,p
        paths += [source,*data["input_sha256"]]

    z,u,beta = primitive_beta(S61())
    by_name = {str(x):x for x in u}
    plane = sym.eye(4)[:,[1,2]]
    # Verify both defining plane equations directly, before using the marked trace.
    restricted = [sym.factor((plane.T*beta.subs(dict(zip(z,plane[:,i])))*plane)[0,1]) for i in range(2)]
    equations = [by_name["u12"],by_name["u13"]]
    system = sym.Matrix(restricted).jacobian(equations)
    assert system.det() and sym.Matrix(restricted)==system*sym.Matrix(equations)
    beta = beta.subs(dict.fromkeys(equations,0))
    l0,l1 = sym.symbols("l0 l1")
    ell = sym.Matrix([l0,0,0,l1])
    matrix = sym.Matrix.hstack(*[J*beta.subs(dict(zip(z,plane[:,i])))*J*ell for i in range(2)])
    assert matrix[[0,3],:] == sym.zeros(2)
    trace_functional = sym.expand(-5*matrix[[1,2],:].trace())
    expected = l0*(5*by_name["u14"]-by_name["u49"]) + l1*(5*by_name["u35"]+by_name["u50"])
    assert trace_functional == sym.expand(expected)

    # A product in the marked relative line has the asserted order-two discriminant.
    a,b,c,d = sym.symbols("a b c d")
    A,B,C = a*c, a*d+b*c, b*d
    assert sym.expand(B**2-4*A*C-(a*d-b*c)**2) == 0
    m = sym.Symbol("m",positive=True,integer=True)
    assert sym.expand((2*m-1)*20+2*m*(-20)) == -20
    paths += ["research/notes/S61_BINARY_DEEPER_REDUCTION_20261101.md",
              "research/notes/PRIMITIVE_JET_ELLIPTIC_CURVE_20261101.md",
              "research/notes/S61_AFFINE_BOUNDARY_DIMENSION_20261102.md",
              "research/notes/S61_BINARY_H11_RANK_DROP_RESIDUAL_FINITE_20261102.md"]
    return dict(id=STEM,all_passed=True,review_kind="Independent intrinsic trace and symbolic local identities, plus proof and source audit",
                intrinsic_functional="-5 trace on Pi of v -> J beta_v J ell, ell in ann(Pi)",
                plane_conditions=list(map(str,restricted)),trace_in_standard_frame=str(trace_functional),
                full_local_parabolic_matrix_rerun=False,line_factor_discriminant_verified=True,
                filtration_degree_after_etale_pullback=-20,generic_true_infinity_dimension_upper=1,
                generic_affine_dimension_upper=2,finite_projection_coordinates=2,affine_exhaustion_claimed=False,
                proof_note=NOTE,input_sha256={p:sha(p) for p in dict.fromkeys(paths)},seconds=time.monotonic()-start)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    target=parser.parse_args().output
    result=verify()
    with target.open("x") as f:
        json.dump(result,f,indent=2)
        f.write("\n")
    print(f"PASS {result['seconds']:.3f}s, {len(result['input_sha256'])} hashes")

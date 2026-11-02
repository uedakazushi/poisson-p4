"""Ambient normal length at the nonzero T7 point above the sparse flag beta."""
import argparse
import hashlib
import json
import time
from functools import lru_cache
from pathlib import Path

import sympy as sym

from s61_model import S61, X, J
from weight1_structure import sp_basis

ROOT = Path(__file__).resolve().parents[1]
STEM = "S61_T7_CRITICAL_AMBIENT_TRANSVERSE_20261102"
NOTE = f"research/notes/{STEM}.md"
FORMAL = "reports/S61_R0_MIXED_SURFACE_T7_FORMAL_20261101.json"
INFINITY = "reports/S61_FLAG_T7_BOUNDARY_20261101.json"
EXTRA = "reports/S61_UNIVERSAL_WEIGHT6_20261031.json"


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def verify():
    start = time.monotonic()
    model = S61()
    boundary = json.loads((ROOT/INFINITY).read_text())
    generic = json.loads((ROOT/FORMAL).read_text())
    for report in (boundary,generic):
        assert report["all_passed"]
        for path,digest in report["input_sha256"].items():
            assert sha(path) == digest
    assert generic["free_rank"] == 4 and generic["homogenizer_zero"] and generic["relative_free_rank"] == 4
    h11 = sym.Matrix(boundary["fixed_beta_kernel"]).applyfunc(sym.Rational)
    a = sym.zeros(11,1)
    a[0] = a[8] = 1
    beta = h11*a
    e = sym.Symbol("e")
    q = sym.Matrix(boundary["projective_q"]).applyfunc(sym.sympify).subs(e,1)
    assert q[24] == 3 and q[25] == -1
    one = [i for i,w in enumerate(model.reduced["parameter_weights"]) if w==1]
    two = [i for i,w in enumerate(model.reduced["parameter_weights"]) if w==2]
    x,y,z,w = X[1:]
    spatial = sym.Matrix(X[1:])
    tensor = sym.zeros(4)
    tensor[1,2],tensor[1,3] = 3*w**2,-x**2
    tensor -= tensor.T
    orbit_columns = []
    for generator in sp_basis():
        velocity = generator*spatial
        dq = generator*tensor + tensor*generator.T - sum((velocity[i]*tensor.diff(spatial[i]) for i in range(4)),sym.zeros(4))
        ambient = X[0]**2*J
        ambient[1:,1:] += dq
        variation = sym.Matrix(model.normalize(ambient)["parameters"])[two,:]
        orbit_columns.append(variation+variation[25]*q)
    orbit = sym.Matrix.hstack(*orbit_columns)
    beta_rows = list(h11.T.to_DM().rref()[1])
    q_rows = list(orbit.T.to_DM().rref()[1])
    group_cols = list(orbit[q_rows,:].to_DM().rref()[1])
    assert len(beta_rows) == 11 and len(q_rows) == len(group_cols) == 7
    # The incidence minor is block triangular since varying H11 fixes Q.
    incidence_minor = h11[beta_rows,:].det()*orbit[q_rows,group_cols].det()
    assert incidence_minor != 0
    point = dict(zip(one,beta)) | dict(zip(two,q)) | {61:sym.S.Zero}
    frozen = {two[25],*[one[i] for i in beta_rows],*[two[i] for i in q_rows]}
    variables = [i for i in range(62) if i not in frozen]
    assert len(variables) == 43 and 61 in variables
    location = {j:i for i,j in enumerate(variables)}
    equations = model.reduced["equations"]+[json.loads((ROOT/EXTRA).read_text())["universal_new_equation"]]
    assert len(equations)==135
    jacobian = sym.zeros(135,43)
    higher = []
    for i,terms in enumerate(equations):
        row = {}
        for mon,numerator,denominator in terms:
            q_degree = sum(k in two for k in mon)
            expansion = {():sym.Rational(numerator,denominator)}
            for k in [*mon,*([61]*(2-q_degree))]:
                new = {}
                for indices,coef in expansion.items():
                    if point[k]:
                        new[indices] = new.get(indices,0)+point[k]*coef
                    if k in location:
                        key = tuple(sorted((*indices,location[k])))
                        new[key] = new.get(key,0)+coef
                expansion = new
            for key,coef in expansion.items():
                row[key] = row.get(key,0)+coef
        row = {key:coef for key,coef in row.items() if coef}
        assert () not in row
        for key,coef in row.items():
            if len(key)==1:
                jacobian[i,key[0]]=coef
        higher.append({key:coef for key,coef in row.items() if len(key)>=2})
    pivots = list(jacobian.to_DM().rref()[1])
    rows = list(jacobian[:,pivots].T.to_DM().rref()[1])
    assert len(pivots)==42
    block = jacobian[rows,pivots]
    inverse = block.inv()
    right = sym.zeros(43,42)
    for i,j in enumerate(pivots):
        right[j,:]=inverse[i,:]
    assert jacobian[rows,:]*right == sym.eye(42)
    tangent = jacobian.nullspace()[0]
    free = next(i for i in range(43) if i not in pivots)
    assert tangent[free]==1
    lifts = {1:tangent}
    residuals = {}
    for degree in range(2,5):
        @lru_cache(None)
        def product(indices,order):
            if not indices:
                return int(order==0)
            return sum(lifts[k][indices[0]]*product(indices[1:],order-k) for k in range(1,min(order,degree-1)+1))
        raw = sym.Matrix([sum(coef*product(indices,degree) for indices,coef in row.items() if len(indices)<=degree) for row in higher])
        lifts[degree] = -right*raw[rows,:]
        residual = raw + jacobian*lifts[degree]
        residuals[degree] = {str(i):str(v) for i,v in enumerate(residual) if v}
        assert not any(residual[i] for i in rows)
    assert residuals[2] == residuals[3] == {}
    assert residuals[4] == {"103":"-267070968/4654064375","123":"801212904/23270321875"}
    paths = [str(Path(__file__).relative_to(ROOT)),NOTE,FORMAL,INFINITY,EXTRA,
             "research/notes/S61_R0_MIXED_SURFACE_T7_FORMAL_20261101.md",
             "research/notes/S61_R0_MIXED_SURFACE_T7_EXHAUSTION_20261101.md",
             "reports/S61_FLAG_T7_WEIGHTED_PRIMARY_20261102.json",*boundary["input_sha256"],*generic["input_sha256"]]
    return dict(id=STEM,all_passed=True,field="QQ",H11_coordinates=list(map(int,a)),beta=list(map(str,beta)),
                q_support={str(i):str(v) for i,v in enumerate(q) if v},
                incidence_dimension=18,beta_cut_indices=beta_rows,q_cut_indices=q_rows,orbit_source_columns=group_cols,
                incidence_minor=str(incidence_minor),transverse_variables=variables,rank=len(pivots),
                jacobian_minor_rows=rows,jacobian_minor_columns=pivots,jacobian_minor=str(block.det(method="domain-ge")),
                normal_raw_index=variables[free],t_index=61,formal_corrections={str(k):list(map(str,col)) for k,col in lifts.items()},
                full_residuals=residuals,transverse_length=4,ambient_formal_parameters=18,
                full_ambient_algebra="QQ[[z1,...,z18]][u]/(u^4)",homogenizer_zero_in_ambient=True,
                flag_endpoint_claim_in_this_report=False,proof_note=NOTE,
                input_sha256={p:sha(p) for p in dict.fromkeys(paths)},seconds=time.monotonic()-start)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    output=parser.parse_args().output
    report=verify()
    with output.open("x") as stream:
        json.dump(report,stream,indent=2)
        stream.write("\n")
    print(f"PASS {report['seconds']:.3f}s, {len(report['input_sha256'])} hashes")

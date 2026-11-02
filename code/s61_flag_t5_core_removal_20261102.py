"""The full ordinary T5 flag core is removed over sufficiently general beta."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import sympy as sym

from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_flag_t5_core_geometry_20261102 import CORE,curve_certificate
from s61_t5_critical_ambient_length_20261102 import ambient_certificate
from s61_model import S61
from weight1_structure import COORDS,infinitesimal_action,sp_basis

ROOT=Path(__file__).resolve().parents[1]
STEM="S61_FLAG_T5_CORE_REMOVAL_20261102"
NOTE=f"research/notes/{STEM}.md"
PRIOR=["S61_R0_T5_FORMAL_RANK23_20261101","S61_R0_T5_ZERO_L_EXHAUSTION_20261101",
       "S61_FLAG_REMAINING_CORE_WEIGHTED_POINTS_20261102","S61_FLAG_T7_CORE_REMOVAL_20261102",
       "S61_FLAG_T7_WEIGHTED_PRIMARY_20261102","S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101",
       "S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101","S61_FLAG_UNIFORM_OPEN_20261101"]


def sha(p):
    return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def verify():
    start=time.monotonic()
    paths=[str(Path(__file__).relative_to(ROOT)),NOTE,"code/s61_flag_t5_core_geometry_20261102.py",
           "code/s61_t5_critical_ambient_length_20261102.py","code/s61_flag_t7_weighted_primary_20261102.py",
           "code/s61_flag_linear_artifact_weighted_engine_20261102.py","code/s61_flag_ordinary_weighted_rows_20261102.cpp",
           "code/s61_nil22_infinity_20261101.py"]
    previous={}
    for stem in PRIOR:
        p=f"reports/{stem}.json"
        report=json.loads((ROOT/p).read_text())
        assert report["all_passed"]
        for source,digest in report["input_sha256"].items():
            assert sha(source)==digest,source
        paths.extend([p,report["proof_note"],*report["input_sha256"]])
        previous[stem]=report
    special=curve_certificate()
    ambient=ambient_certificate()
    assert ambient["actual_transverse_length"]==23
    known=previous["S61_R0_T5_FORMAL_RANK23_20261101"]
    assert known["relative"]["finite_free_rank"]==known["relative"]["homogenizer_power_zero"]==23
    assert known["relative"]["smooth_equivariant_descent_to_ambient"]
    m=S61()
    one=[i for i,w in enumerate(m.reduced["parameter_weights"]) if w==1]
    point=sym.zeros(61,1)
    point[4]=point[48]=1
    point[12]=CORE
    nuisance=[24,33,34,35,36,47,49]
    for k,i in enumerate(nuisance):
        point[i]=k+1
    z,u,beta=primitive_beta(m)
    vector=sym.Matrix([beta[i,j].diff(z[k]) for (i,j),k in COORDS])
    representation=vector.jacobian(u)
    rows=list(representation.T.to_DM().rref()[1])
    inverse=representation[rows,:].inv()
    action=sym.Matrix.hstack(*[inverse*infinitesimal_action(vector,g)[rows,:] for g in sp_basis()])
    action=action.subs(dict(zip(u,[point[i] for i in one])))
    source=sym.eye(16)[:,[one.index(j) for j in nuisance]].row_join(action)
    columns=list(source.to_DM().rref()[1])
    assert len(columns)==16
    determinant=source[:,columns].det(method="domain-ge")
    assert determinant and 2-CORE==sym.Rational(56,25)
    return dict(id=STEM,all_passed=True,core=str(CORE),fixed_beta_geometry=special,ambient_nonzero_point=ambient,
                actual_power_module="t^23 S",finite_over_full_H14=True,relative_nilpotence_exponent_determined=False,
                fixed_beta_orbit_minor=2,relative_compensator_minor="56/25",nuisance_weights=[3,3,1,2,2,2,1],
                source_rank=16,source_minor_columns=columns,source_minor=str(determinant),
                source_beta=[str(point[i]) for i in one],all_sufficiently_general_beta_T5_core_removed=True,
                core_two_or_full_S61_exhaustion_claimed=False,proof_note=NOTE,
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

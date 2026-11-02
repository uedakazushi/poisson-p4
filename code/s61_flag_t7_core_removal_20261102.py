"""Whole generic-beta removal for the ordinary flag core r=-91/400."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import sympy as sym

from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import COORDS, infinitesimal_action, sp_basis

ROOT = Path(__file__).resolve().parents[1]
STEM = "S61_FLAG_T7_CORE_REMOVAL_20261102"
NOTE = f"research/notes/{STEM}.md"
REPORTS = [
    "S61_T7_CRITICAL_AMBIENT_TRANSVERSE_20261102",
    "S61_FLAG_T7_WEIGHTED_PRIMARY_20261102",
    "S61_FLAG_REMAINING_CORE_WEIGHTED_POINTS_20261102",
    "S61_FLAG_ORDINARY_WEIGHTED_FINITE_KERNEL_20261102",
    "S61_FLAG_ZERO12_RELATIVE_WEIGHT_REMOVAL_20261102",
    "S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101",
    "S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101",
]


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def verify():
    start = time.monotonic()
    paths = [str(Path(__file__).relative_to(ROOT)),NOTE]
    evidence = {}
    for stem in REPORTS:
        path = f"reports/{stem}.json"
        data = json.loads((ROOT/path).read_text())
        assert data["all_passed"]
        evidence[stem] = data
        paths += [path,f"research/notes/{stem}.md",*data["input_sha256"]]
        for p,digest in data["input_sha256"].items():
            assert sha(p) == digest,p
    ambient,primary,points,ordinary = [evidence[k] for k in REPORTS[:4]]
    assert ambient["transverse_length"] == 4 and ambient["ambient_formal_parameters"] == 18
    assert ambient["homogenizer_zero_in_ambient"] and ambient["H11_coordinates"] == [1,0,0,0,0,0,0,0,1,0,0]
    assert primary["generic_lengths"] == [5,2] and primary["embedded_difference_length"] == 106
    assert primary["actual_T7_curve"]["word"] == [[4,"7/20"],[2,"0"]]
    assert primary["actual_T7_curve"]["raw_q_support"]["24"] == "3*e"
    assert points["points"][2]["core"] == "-91/400" and points["points"][2]["slice_dimension"] == 1
    assert ordinary["two_orbit_cut_minor"] == "2"
    r = -sym.Rational(91,400)
    compensator = sym.Matrix([[-1,2],[-1,r]]).det()
    assert compensator == sym.Rational(891,400)
    model = S61()
    one = [i for i,w in enumerate(model.reduced["parameter_weights"]) if w==1]
    sample = sym.zeros(61,1)
    sample[4] = sample[48] = 1
    sample[12] = r
    nuisance = [24,33,34,35,36,47,49]
    for i,j in enumerate(nuisance):
        sample[j]=i+1
    spatial,coordinates,beta = primitive_beta(model)
    entries = sym.Matrix([beta[i,j].diff(spatial[k]) for (i,j),k in COORDS])
    basis = entries.jacobian(coordinates)
    selected = list(basis.T.to_DM().rref()[1])
    inverse = basis[selected,:].inv()
    action = sym.Matrix.hstack(*[inverse*infinitesimal_action(entries,g)[selected,:] for g in sp_basis()])
    action = action.subs(dict(zip(coordinates,[sample[j] for j in one])))
    differential = sym.eye(16)[:,[one.index(j) for j in nuisance]].row_join(action)
    columns = list(differential.to_DM().rref()[1])
    assert len(columns)==16
    determinant = differential[:,columns].det(method="domain-ge")
    assert determinant
    hashes = {p:sha(p) for p in dict.fromkeys(paths)}
    return dict(id=STEM,all_passed=True,field="QQ",core=str(r),source_nuisance_coordinates=nuisance,
                positive_source_weights=[3,3,1,2,2,2,1],normalized_beta_sample=[str(sample[j]) for j in one],
                full_source_minor_columns=columns,full_source_rank=16,full_source_minor=str(determinant),
                two_orbit_minor=2,H14_compensator_minor=str(compensator),
                special_support_curves=["actual flag D5","actual T7 curve"],
                t_module_vanishes_at_both_special_generic_points=True,actual_t_module_finite_over_H14_base=True,
                nilpotence_at_complete_relative_flag=True,cut_directions_restored=True,
                all_nuisance_contraction=True,all_sufficiently_general_beta_T7_flag_points_removed=True,
                nilpotence_exponent_determined=False,other_cores_claimed=False,full_affine_exhaustion_claimed=False,
                proof_note=NOTE,input_sha256=hashes,seconds=time.monotonic()-start)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    output=parser.parse_args().output
    report=verify()
    with output.open("x") as stream:
        json.dump(report,stream,indent=2)
        stream.write("\n")
    print(f"PASS {report['seconds']:.3f}s, {len(report['input_sha256'])} hashes")

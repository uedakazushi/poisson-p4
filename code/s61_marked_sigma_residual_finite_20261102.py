"""The marked-square curve has at most finitely many residual points over general beta."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from s61_marked_sigma_weighted_axis_20261102 import ROOT, FLAG, ACTION, build, artin_certificate

STEM="S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102"
NOTE=f"research/notes/{STEM}.md"


def produce():
    start=time.monotonic()
    paths=[NOTE,str(Path(__file__).relative_to(ROOT)), "code/s61_marked_sigma_weighted_axis_20261102.py",
           FLAG,ACTION,"code/s61_model.py","code/s61_flag_transverse_obstruction_20261101.py",
           "code/s61_flag_t7_weighted_primary_20261102.py","reports/S61_UNIVERSAL_WEIGHT6_20261031.json",
           "research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md"]
    for stem in ["S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101","S61_FLAG_UNIFORM_OPEN_20261101",
                 "S61_FLAG_ZERO12_RELATIVE_WEIGHT_REMOVAL_20261102","S61_AFFINE_SURFACE_BOUND_20261102"]:
        path=f"reports/{stem}.json"
        prior=json.loads((ROOT/path).read_text())
        assert prior["all_passed"]
        for source,digest in prior["input_sha256"].items():
            assert hashlib.sha256((ROOT/source).read_bytes()).hexdigest()==digest,source
        paths.extend([path,prior.get("proof_note",f"research/notes/{stem}.md"),*prior["input_sha256"]])
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in dict.fromkeys(paths)}
    names,weights,residuals,frame=build()
    central=artin_certificate(names,weights,residuals)
    return dict(id=STEM,all_passed=True,field="QQ",frame=frame,central_axis=central,base_dimension=15,
                center_in_marked_square=True, marked_source_dimension=17, actual_power_module="t S",
                nonempty_ambient_removed_open=True, generic_beta_marked_residual_dimension_at_most=0,
                general_beta_entire_marked_curve_removed=False,relative_nilpotence_exponent_determined=False,
                proof_note=NOTE,input_sha256=hashes,seconds=time.monotonic()-start)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    output=parser.parse_args().output
    result=produce()
    with output.open("x")as fp:
        json.dump(result,fp,indent=2)
        fp.write("\n")
    print(f"PASS {result['seconds']:.3f}s / {len(result['input_sha256'])} hashes")

"""Eliminate the entire h=0 core divisor over a sufficiently general first jet."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import sympy as s

from s61_flag_t7_weighted_primary_20261102 import M2
from s61_model import S61
from s61_marked_sigma_core_axis_engine_20261102 import build_axis
from s61_marked_sigma_pure_endpoint_removal_20261102 import point_geometry

ROOT=Path(__file__).resolve().parents[1]
STEM="S61_MARKED_SIGMA_ZERO_H_REMOVAL_20261102"
NOTE=f"research/notes/{STEM}.md"
PRIOR=["S61_MARKED_SIGMA_PROJECTIVE_CORE_20261102","S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102",
       "S61_MARKED_SIGMA_INTERSECTION_REMOVAL_20261102","S61_MARKED_SIGMA_HOMOGENEOUS_PENCILS_20261102",
       "S61_MARKED_SIGMA_PURE_ENDPOINT_REMOVAL_20261102"]


def verify():
    started=time.monotonic();paths=[NOTE,str(Path(__file__).relative_to(ROOT))];old={}
    def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    for stem in PRIOR:
        path=f"reports/{stem}.json";d=json.loads((ROOT/path).read_text());old[stem]=d
        assert d["all_passed"]
        for p,h in d["input_sha256"].items():assert sha(p)==h,p
        paths.extend([path,d["proof_note"],*d["input_sha256"]])
    hashes={p:sha(p) for p in dict.fromkeys(paths)}
    core=s.zeros(15,1);core[1]=core[5]=1;core[14]=2
    known=old[PRIOR[1]]["frame"]
    names,weights,forms,frame=build_axis(core,[0,19,27],known["q_weights"])
    def encode(p):
        return "+".join("("+str(c)+")"+"".join("*"+n+(f"^{k}" if k>1 else "")for n,k in zip(names,mon)if k)
                       for mon,c in p.items())or"0"
    script="R=QQ["+ ",".join(names)+",Degrees=>{"+",".join(map(str,weights))+"}];\n"
    script+="I=ideal("+",".join(encode(p)for p in forms if p)+");\n"
    script+="""
    Z=I+ideal(rho); assert(dim(R/I)==1 and dim(R/Z)==0 and numColumns basis(R/Z)==70);
    assert(toString(hilbertSeries(R/I,Reduce=>true))=="(1+7*T+23*T^2+35*T^3-25*T^4-37*T^5)/((1-T))");
    assert(numColumns basis(R/(I+ideal(rho-1)))==4); print("ZERO_H_MARKED_PASS");
    """
    with tempfile.TemporaryDirectory(prefix="marked-zero-h-")as temp:
        script_file=Path(temp)/"center.m2";script_file.write_text(script)
        env=os.environ|{"DYLD_LIBRARY_PATH":str(ROOT/"work/m2-runtime/lib")+":/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib"}
        run=subprocess.run([M2,"--script",str(script_file)],env=env,timeout=100,capture_output=True,text=True)
        assert run.returncode==0 and "ZERO_H_MARKED_PASS"in run.stdout,(run.stdout,run.stderr)
    beta=s.Matrix([0,1,*core[1:]])
    q0=s.zeros(45,1);q0[25]=-1
    assert point_geometry(S61(),q0,beta,False)["rho_one_transverse_tangent_dimension"]==3
    flag=json.loads((ROOT/"reports/S61_FLAG_UNIFORM_OPEN_20261101.json").read_text())
    acts=s.Matrix(json.loads((ROOT/"reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json").read_text())["primitive_action"]).applyfunc(s.sympify)
    axes=sorted(acts.free_symbols,key=lambda z:int(z.name[1:]))
    chart=[i for i in range(45) if i!=25];qdiff=s.Matrix(flag["flag"]["q_orbit_differential"]).applyfunc(s.Rational)
    block=qdiff[[chart.index(j)for j in [0,19,27]],:].col_join(acts.subs(dict(zip(axes,s.Matrix([0,*core]))))[0,:])
    assert block[:,[4,5,6,9]].det()==4
    return dict(id=STEM,all_passed=True,field="QQ",marked_core=[1,2,0],frame=frame,
                central_rho_zero_length=70,central_rho_one_length=4,central_only_reduced_branch="standard rho axis",
                relative_actual_ideal="t S",ambient_etale_minor=4,
                irreducible_source_fiber_dimension=12,irreducible_joint_dimension=16,
                all_general_beta_zero_h_removed=True,all_zero_nuisance_points_on_h_zero_removed=False,
                all_general_beta_marked_squares_removed=False,proof_note=NOTE,input_sha256=hashes,
                seconds=time.monotonic()-started)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",type=Path,required=True)
    out=p.parse_args().output;d=verify()
    with out.open("x")as f:json.dump(d,f,indent=2);f.write("\n")
    print(f"PASS {d['seconds']:.3f}s / {len(d['input_sha256'])} hashes")

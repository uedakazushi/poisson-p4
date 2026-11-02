"""The point [0:0:1] completing c=2a has an ambient removed neighborhood."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_flag_t7_weighted_primary_20261102 import M2

ROOT=Path(__file__).resolve().parents[1]
STEM="S61_MARKED_SIGMA_PENCIL_INFINITY_20261102"
NOTE=f"research/notes/{STEM}.md"


def verify():
    t0=time.monotonic();paths=[NOTE,str(Path(__file__).relative_to(ROOT))]
    def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
    old={}
    for stem in ["S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102","S61_MARKED_SIGMA_PROJECTIVE_CORE_20261102",
                 "S61_MARKED_SIGMA_PENCIL_MODEL_20261102"]:
        p=f"reports/{stem}.json";d=json.loads((ROOT/p).read_text());old[stem]=d;assert d["all_passed"]
        for k,v in d["input_sha256"].items():assert sha(k)==v,k
        paths.extend([p,d["proof_note"],*d["input_sha256"]])
    refs={p:sha(p)for p in dict.fromkeys(paths)}
    model=S61();one=[i for i,v in enumerate(model.reduced["parameter_weights"])if v==1]
    chart=[j for j in range(45)if j!=25];cuts=[0,9,27,32]
    cols=[j for j in chart if j not in cuts]
    qw=old["S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102"]["frame"]["q_weights"]
    rho=s.Symbol("rho")
    beta=dict.fromkeys(one,s.S.Zero);beta.update({4:rho,11:1,26:1})
    _,forms=specialize(model,beta,True)
    variables=[f"q{j}"for j in cols]+["t","rho"]
    weights=[qw[chart.index(j)]for j in cols]+[3,1]
    assert len(variables)==42 and min(weights)>0
    polys=[]
    for form in forms:
        terms=[];degrees=set()
        for mon,coeff in form.items():
            if any(j in cuts for j in mon):continue
            coeff*=(-1)**mon.count(25)
            subs=[f"q{j}"for j in mon if j!=25]+["t"]*(2-len(mon))
            for (d,),c in s.Poly(coeff,rho).terms():
                if c:degrees.add(d+sum(weights[variables.index(x)]for x in subs))
            if coeff:terms.append("("+str(coeff).replace("**","^")+")"+"".join("*"+x for x in subs))
        assert len(degrees)<=1
        if terms:polys.append("+".join(terms))
    script="R=QQ["+ ",".join(variables)+",Degrees=>{"+",".join(map(str,weights))+"}];\n"
    script+="I=ideal("+",".join(polys)+");\n"
    script+="""
    U=I+ideal(rho); assert(dim(R/U)==0 and numColumns basis(R/U)==79);
    assert(toString(hilbertSeries(R/I,Reduce=>true))=="(1+7*T+23*T^2+35*T^3-24*T^4-36*T^5-2*T^6)/((1-T))");
    assert(numColumns basis(R/(I+ideal(rho-1)))==4 and t%I!=0 and t^2%I==0);
    print("INFINITY_AXIS_PASS");
    """
    with tempfile.TemporaryDirectory(prefix="marked-pencil-infinity-")as tmp:
        p=Path(tmp)/"raw.m2";p.write_text(script)
        env=os.environ|{"DYLD_LIBRARY_PATH":str(ROOT/"work/m2-runtime/lib")+":/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib"}
        run=subprocess.run([M2,"--script",str(p)],env=env,capture_output=True,text=True,timeout=90)
        assert run.returncode==0 and "INFINITY_AXIS_PASS" in run.stdout,(run.stdout,run.stderr)
    J=s.Matrix([[-f.get(tuple(sorted((25,j))),0) for j in cols]+[-f.get((25,),0)] for f in forms]).subs(rho,1)
    assert J.cols-J.to_DM().rank()==3
    flag=json.loads((ROOT/"reports/S61_FLAG_UNIFORM_OPEN_20261101.json").read_text())["flag"]
    act=s.Matrix(json.loads((ROOT/"reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json").read_text())["primitive_action"]).applyfunc(s.sympify)
    us=sorted(act.free_symbols,key=lambda z:int(z.name[1:]))
    at=act.subs(dict(zip(us,[s.sympify(beta[j]).subs(rho,0) for j in one])))
    orbit=s.Matrix(flag["q_orbit_differential"]).applyfunc(s.Rational)
    block=orbit[[chart.index(j)for j in [0,9,27]],:].col_join(at[0,:])
    assert block[:,[4,5,6,9]].det()==-2
    return dict(id=STEM,all_passed=True,field="QQ",projective_core=[0,0,1],raw_weighted_variable_count=42,
                rho_zero_length=79,rho_one_length=4,standard_flag_rho_one_tangent_dimension=3,
                only_reduced_central_branch="standard flag rho axis",relative_actual_ideal="t S",
                full_ambient_etale_minor=-2,zero_nuisance_full_ambient_removed=True,
                all_nuisance_at_this_core_removed=True,all_finite_h_removed=False,
                proof_note=NOTE,input_sha256=refs,seconds=time.monotonic()-t0)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",type=Path,required=True)
    out=p.parse_args().output;d=verify()
    with out.open("x")as f:json.dump(d,f,indent=2);f.write("\n")
    print(f"PASS {d['seconds']:.3f}s / {len(d['input_sha256'])} hashes")

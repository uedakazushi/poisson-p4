"""Ambient removal of every b1=0 marked square over sufficiently general beta."""
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
from s61_marked_sigma_core_axis_engine_20261102 import build_axis

ROOT=Path(__file__).resolve().parents[1]
STEM="S61_MARKED_SIGMA_INTERSECTION_REMOVAL_20261102"
NOTE=f"research/notes/{STEM}.md"
PRIOR=["S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102","S61_FLAG_TRUE_CORE_WEIGHTED_REMOVAL_20261102",
       "S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101","S61_FLAG_UNIFORM_OPEN_20261101",
       "S61_AFFINE_SURFACE_BOUND_20261102"]


def verify():
    start=time.monotonic()
    paths=[str(Path(__file__).relative_to(ROOT)),NOTE,"code/s61_marked_sigma_core_axis_engine_20261102.py",
           "reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json"]
    primary={}
    def sha(path):
        return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    for stem in PRIOR:
        p=f"reports/{stem}.json";primary[stem]=data=json.loads((ROOT/p).read_text())
        assert data["all_passed"]
        for path,digest in data["input_sha256"].items(): assert sha(path)==digest,path
        paths.extend([p,data.get("proof_note",f"research/notes/{stem}.md"),*data["input_sha256"]])
    hashes={p:sha(p) for p in dict.fromkeys(paths)}
    known=primary[PRIOR[0]]
    weights=known["frame"]["q_weights"]
    core=s.zeros(15,1)
    core[3],core[12],core[5],core[6],core[14]=s.Rational(4,5),1,1,3,2
    names,ws,polynomials,frame=build_axis(core,[9,19,27],weights)
    assert len(names)==14 and frame["Jacobian_rank"]==28
    def encode(p):
        return "+".join(f"({c})"+"".join("*"+n+(f"^{i}" if i>1 else "")
                   for n,i in zip(names,mon) if i) for mon,c in p.items())or"0"
    program="R=QQ["+ ",".join(names)+",Degrees=>{"+",".join(map(str,ws))+"}];\n"
    program+="I=ideal("+",".join(encode(p)for p in polynomials if p)+");\n"
    program+="""
P=ideal apply(0..12,i->R_i); A=saturate(I,ideal(rho)); U=I+ideal(rho);
assert(dim(R/I)==1 and dim(R/U)==0 and numColumns basis(R/U)==76);
assert(radical A==P and A==saturate(A,ideal(rho)));
assert(numColumns basis(R/(A+ideal(rho)))==31 and v12^2%A==0);
assert(toString(hilbertSeries(R/I,Reduce=>true))=="(1+7*T+23*T^2+35*T^3-8*T^4-22*T^5-5*T^6)/((1-T))");
assert(v12%I!=0 and v12^2%I==0); print("BOTH_MARKED_PASS");
"""
    env=os.environ|{"DYLD_LIBRARY_PATH":str(ROOT/"work/m2-runtime/lib")+":/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib"}
    with tempfile.TemporaryDirectory(prefix="marked-intersection-axis-")as temp:
        p=Path(temp)/"section.m2";p.write_text(program)
        run=subprocess.run([M2,"--script",str(p)],capture_output=True,text=True,env=env,timeout=110)
        assert run.returncode==0 and "BOTH_MARKED_PASS"in run.stdout,(run.stdout,run.stderr)
    flag=primary["S61_FLAG_UNIFORM_OPEN_20261101"]
    H=s.Matrix(flag["flag"]["beta_kernel"]).applyfunc(s.Rational)
    d=s.Matrix(flag["frame_change"]["beta_derivative"]).applyfunc(s.sympify)
    b=s.symbols("b0:15");moving=core+s.eye(15)[:,0];stages=[]
    for j,parameter in [(4,lambda y:-y[3]/y[0]),(2,lambda y:(y[14]-y[5])/(3*y[0])),
                         (3,lambda y:y[5]/y[12]),(0,lambda y:-y[6]/y[0])]:
        n=d[:,j].jacobian(b);assert n**4==s.zeros(15)
        t=parameter(moving)
        moving=sum((t**k/s.factorial(k)*n**k for k in range(4)),s.zeros(15))*moving
        stages.append([j,str(t)])
    assert moving[0]==moving[12]==1 and all(moving[j]==0 for j in [1,3,5,6,14])
    assert moving[2]==-s.Rational(4,25)
    acts=s.Matrix(json.loads((ROOT/"reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json").read_text())["primitive_action"]).applyfunc(s.sympify)
    symbols=sorted(acts.free_symbols,key=lambda z:int(z.name[1:]))
    act=acts.subs(dict(zip(symbols,H*core)))
    chart=[i for i in range(45)if i!=25]
    orbit=s.Matrix(flag["flag"]["q_orbit_differential"]).applyfunc(s.Rational)
    four=orbit[[chart.index(i) for i in [9,19,27]],:].col_join(act[0,:])
    cols=list(four.to_DM().rref()[1]);assert len(cols)==4
    minor=four[:,cols].det()
    origin=s.zeros(15,1);origin[0]=origin[12]=1;origin[2]=-s.Rational(4,25)
    at_zero=acts.subs(dict(zip(symbols,H*origin)))
    compensators=at_zero[[0,2],[5,6]]
    assert abs(compensators.det())==s.Rational(54,25)
    dependent=s.Matrix(b);dependent[0]=dependent[1]=0;dependent[3]=s.Rational(4,5)*b[12]
    parameters=sorted(dependent.free_symbols,key=str)
    incid=H*dependent
    source=incid.jacobian(parameters).row_join(acts.subs(dict(zip(symbols,incid))))
    sm=source.subs(dict(zip(parameters,range(2,2+len(parameters)))))
    chosen=list(sm.to_DM().rref()[1]);assert len(chosen)==16
    return dict(id=STEM,all_passed=True,field="QQ",core=list(map(str,core)),frame=frame,
                central_axis_unique_reduced_prime=True,central_axis_generic_thickness=31,
                central_rho_zero_length=76,central_nonzero_t_square_zero=True,
                rho_one_finite_flag_transform=stages,rho_one_normalized_b=list(map(str,moving)),
                trueD_full_H14_compensator_minor=str(compensators.det()),relative_actual_module="t^2 S",
                ambient_orbit_minor=str(minor),orbit_columns=cols,
                restricted_source_dimension=16,restricted_source_beta_minor=str(sm[:,chosen].det(method="domain-ge")),
                all_general_beta_marked_b1_zero_removed=True,all_marked_squares_removed=False,
                proof_note=NOTE,input_sha256=hashes,seconds=time.monotonic()-start)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--output",type=Path,required=True)
    out=parser.parse_args().output;result=verify()
    with out.open("x")as f:json.dump(result,f,indent=2);f.write("\n")
    print(f"PASS {result['seconds']:.3f}s / {len(result['input_sha256'])} hashes")

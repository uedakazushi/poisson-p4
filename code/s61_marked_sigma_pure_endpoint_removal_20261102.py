"""Eliminate the four pure-cubic minimal flag endpoint presentations over general beta."""
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
from s61_model import S61, X, J as J5
from s61_infinity_flags import leading_equations_matrix
from s61_marked_sigma_core_axis_engine_20261102 import build_axis
from weight1_structure import sp_basis

ROOT=Path(__file__).resolve().parents[1]
STEM="S61_MARKED_SIGMA_PURE_ENDPOINT_REMOVAL_20261102"
NOTE=f"research/notes/{STEM}.md"
CANDIDATE="code/certificates/s61_marked_sigma_pure_endpoint_points_20261102.json"
EXTRA="reports/S61_UNIVERSAL_WEIGHT6_20261031.json"
PRIOR=["S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102","S61_MARKED_SIGMA_INTERSECTION_REMOVAL_20261102",
       "S61_MINIMAL_FLAG_ENDPOINT_INVARIANTS_20261101","S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101",
       "S61_MARKED_SIGMA_PROJECTIVE_CORE_20261102"]


def point_geometry(model,q,beta,need_incidence):
    one=[i for i,w in enumerate(model.reduced["parameter_weights"])if w==1]
    two=[i for i,w in enumerate(model.reduced["parameter_weights"])if w==2]
    chart=[i for i in range(45)if i!=25]
    point=s.zeros(61,1)
    for i,z in zip(two,q):point[i]=z
    Q=model.reconstruct(list(point))-X[0]**2*J5
    for i,z in zip(one,beta):point[i]=z
    assert q[25]==-1 and all(q[i]==0 for i in [0,19,27,32])
    variables=one+[two[j] for j in chart]+[61]
    values=list(point)+[s.S.Zero]
    rows=model.reduced["equations"]+[json.loads((ROOT/EXTRA).read_text())["universal_new_equation"]]
    jac=s.zeros(135,len(variables))
    for i,row in enumerate(rows):
        total=0
        for mon,num,den in row:
            indices=[*mon,*([61]*(2-sum(j in two for j in mon)))]
            c=s.Rational(num,den)
            total+=c*s.prod(values[j]for j in indices)
            for k,j in enumerate(indices):
                if j in variables:
                    jac[i,variables.index(j)]+=c*s.prod(values[z] for z in indices[:k]+indices[k+1:])
        assert total==0
    selected=[16+chart.index(j) for j in chart if j not in [0,19,27,32]]+[60]
    fiber_rank=jac[:,selected].to_DM().rank()
    if not need_incidence:
        assert 41-fiber_rank==3
        return dict(rho_one_transverse_tangent_dimension=3)
    # The orbit and the entire linear kernel, together, give a smooth
    # t=0 incidence of the same dimension as the total tangent space.
    leading=leading_equations_matrix(model,point)
    assert leading[:,0]==s.zeros(leading.rows,1)
    H=s.Matrix.hstack(*leading[:,1:].nullspace())
    assert H.cols==9 and leading[:,1:]*H==s.zeros(leading.rows,9)
    assert leading[:,1:]*beta==s.zeros(leading.rows,1)
    form=s.zeros(16)
    for mon,num,den in rows[-1]:
        qq=[two.index(j) for j in mon if j in two]
        if len(qq)!=2:continue
        bb=[one.index(j)for j in mon if j in one]
        assert len(bb)==2
        form[bb[0],bb[1]]+=s.Rational(num,den)*q[qq[0]]*q[qq[1]]
    assert H.T*(form+form.T)*H==s.zeros(9)
    velocity=[]
    for g in sp_basis():
        A=s.zeros(5);A[1:,1:]=g
        d=A*s.Matrix(X)
        derivative=A*Q+Q*A.T-sum((d[i]*Q.diff(X[i])for i in range(5)),s.zeros(5))
        coords=s.Matrix(model.normalize(X[0]**2*J5+derivative)["parameters"])[two,:]
        velocity.append((coords+coords[25]*q)[chart,:])
    orbit=s.Matrix.hstack(*velocity)
    def witness(M):
        pc=list(M.to_DM().rref()[1]);rs=list(M[:,pc].T.to_DM().rref()[1])
        return dict(columns=pc,rows=rs,determinant=str(M[rs,pc].det(method="domain-ge")))
    orbit_w=witness(orbit)
    ambient_w=witness(jac)
    assert len(orbit_w["columns"])==8 and len(ambient_w["columns"])==44 and fiber_rank==41
    return dict(H_dimension=9,orbit_minor=orbit_w,total_jacobian_minor=ambient_w,
                full_smooth_dimension=17,full_ambient_t_zero=True,rho_one_simple_point=True)


def verify():
    start=time.monotonic()
    paths=[NOTE,str(Path(__file__).relative_to(ROOT)),CANDIDATE,EXTRA,"code/s61_marked_sigma_core_axis_engine_20261102.py",
           "code/s61_infinity_flags.py","code/s61_model.py","code/weight1_structure.py"]
    def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
    primary={}
    for stem in PRIOR:
        p=f"reports/{stem}.json";d=json.loads((ROOT/p).read_text());primary[stem]=d
        assert d["all_passed"]
        for k,v in d["input_sha256"].items():assert sha(k)==v,k
        paths.extend([p,d.get("proof_note",f"research/notes/{stem}.md"),*d["input_sha256"]])
    h={p:sha(p) for p in dict.fromkeys(paths)}
    core=s.zeros(15,1);core[1]=core[5]=1;core[14]=2;core[6]=4
    assert core[6]**2+16*core[5]**2*(core[14]-3*core[5])==0
    ref=primary[PRIOR[0]]["frame"]
    names,weights,polys,frame=build_axis(core,[0,19,27],ref["q_weights"])
    def encode(p):
        return "+".join("("+str(c)+")"+"".join("*"+n+(f"^{k}" if k>1 else "")for n,k in zip(names,mon)if k)
                       for mon,c in p.items())or"0"
    source="R=QQ["+ ",".join(names)+",Degrees=>{"+",".join(map(str,weights))+"}];\n"
    source+="I=ideal("+",".join(encode(p)for p in polys if p)+");\n"
    source+="""
    assert(dim(R/I)==1 and dim(R/(I+ideal(rho)))==0 and numColumns basis(R/(I+ideal(rho)))==70);
    assert(toString(hilbertSeries(R/I,Reduce=>true))=="(1+7*T+23*T^2+35*T^3-30*T^4-30*T^5)/((1-T))");
    assert(numColumns basis(R/(I+ideal(rho-1)))==6);
    print("PURE_MARKED_CENTRAL_PASS");
    """
    with tempfile.TemporaryDirectory(prefix="marked-pure-flag-")as temp:
        p=Path(temp)/"model.m2";p.write_text(source)
        env=os.environ|{"DYLD_LIBRARY_PATH":str(ROOT/"work/m2-runtime/lib")+":/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib"}
        run=subprocess.run([M2,"--script",str(p)],env=env,timeout=115,capture_output=True,text=True)
        assert run.returncode==0 and "PURE_MARKED_CENTRAL_PASS"in run.stdout,(run.stdout,run.stderr)
    model=S61()
    beta=s.Matrix([0,1,*core[1:]])
    q0=s.zeros(45,1);q0[25]=-1
    qs=[s.Matrix(p).applyfunc(s.Rational) for p in json.loads((ROOT/CANDIDATE).read_text())["q_points"]]
    assert len(qs)==2 and qs[0]!=qs[1] and all(q!=q0 for q in qs)
    geo=[point_geometry(model,q,beta,i!=0)for i,q in enumerate([q0,*qs])]
    # For the four ambient orbit directions the determinant equals 4*a;
    # verify at this core, not just on the preceding neighboring section.
    flag=json.loads((ROOT/"reports/S61_FLAG_UNIFORM_OPEN_20261101.json").read_text())
    acts=s.Matrix(json.loads((ROOT/"reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json").read_text())["primitive_action"]).applyfunc(s.sympify)
    bs=sorted(acts.free_symbols,key=lambda z:int(z.name[1:]))
    act=acts.subs(dict(zip(bs,s.Matrix([0,*core]))))
    qchart=[i for i in range(45)if i!=25]
    QD=s.Matrix(flag["flag"]["q_orbit_differential"]).applyfunc(s.Rational)
    block=QD[[qchart.index(j)for j in [0,19,27]],:].col_join(act[0,:])
    assert block[:,[4,5,6,9]].det()==4
    return dict(id=STEM,all_passed=True,field="QQ",central_core=[1,2,4],frame=frame,
                central_dimension=1,rho_zero_length=70,rho_one_length=6,three_point_lengths=[4,1,1],
                finite_fiber_point_geometry=geo,base_dimension=15,actual_module="t S",ambient_orbit_minor=4,
                irreducible_joint_dimension=16,all_generic_beta_pure_cubic_flag_endpoints_removed=True,
                all_generic_beta_marked_squares_removed=False,proof_note=NOTE,input_sha256=h,seconds=time.monotonic()-start)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--output",type=Path,required=True)
    out=parser.parse_args().output;r=verify()
    with out.open("x")as f:json.dump(r,f,indent=2);f.write("\n")
    print(f"PASS {r['seconds']:.3f}s / {len(r['input_sha256'])} hashes")

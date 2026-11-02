"""Recover the two proper primary factors and actual T5 line at the flag."""
import json
import os
import subprocess
import tempfile
from pathlib import Path

import sympy as sym

from s61_flag_linear_artifact_weighted_engine_20261102 import build
from s61_flag_t7_weighted_primary_20261102 import M2, residuals
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_model import X,J
from s61_nil22_infinity_20261101 import beta_matrix
from weight1_structure import sp_basis

ROOT=Path(__file__).resolve().parents[1]
CORE=-sym.Rational(6,25)
POINTS="reports/S61_FLAG_REMAINING_CORE_WEIGHTED_POINTS_20261102.json"
FLAG="reports/S61_FLAG_UNIFORM_OPEN_20261101.json"


def fixed_representative(model):
    e=sym.Symbol("e")
    x,y,z,w=X[1:]
    q=sym.zeros(4)
    q[1,2],q[1,3]=2*e*x*w,-x**2
    q-=q.T
    linear=beta_matrix(model,q)
    reduced,pivots=linear.to_DM().to_field().rref()
    equations=reduced.to_Matrix()[:len(pivots),:]
    assert len(pivots)==4 and not any(a.has(e) for a in equations)
    kernel=sym.Matrix.hstack(*equations.nullspace())
    assert linear*kernel==sym.zeros(linear.rows,12)
    one=[i for i,w in enumerate(model.reduced["parameter_weights"]) if w==1]
    beta=sym.zeros(16,1)
    beta[one.index(4)]=beta[one.index(48)]=1
    beta[one.index(13)]=sym.Rational(2,5)
    assert equations*beta==sym.zeros(4,1)
    # The displayed nonzero line really has the single T5 projective orbit.
    scale=sym.Symbol("a", nonzero=True)
    g=sym.diag(scale,1/scale,1,1)
    spatial=sym.Matrix(X[1:])
    moving=g.inv()*q.subs(e,1).subs(dict(zip(spatial,g*spatial)),simultaneous=True)*g.inv().T
    assert (moving-scale**3*q.subs(e,1/scale)).applyfunc(sym.expand)==sym.zeros(4)
    return q,kernel,beta


def identify_line(model,frame,tangent,inverse):
    data=json.loads((ROOT/FLAG).read_text())
    stabilizer=sym.Matrix(data["frame_change"]["stabilizer"]).applyfunc(sym.sympify)
    beta_action=sym.Matrix(data["frame_change"]["beta_derivative"]).applyfunc(sym.sympify)
    b=sym.symbols("b0:15")
    values=sym.zeros(15,1)
    values[0]=values[12]=1
    values[2]=CORE
    t=sym.Rational(2,5)
    deriv=beta_action[:,4].jacobian(b)
    normalized=sum((t**j/sym.factorial(j)*deriv**j for j in range(4)),sym.zeros(15))*values
    assert deriv**4==sym.zeros(15)
    space=sum((stabilizer[k,4]*g for k,g in enumerate(sp_basis())),sym.zeros(4))
    change=sym.eye(4)+t*space
    q,kernel,beta=fixed_representative(model)
    flag_kernel=sym.Matrix(data["flag"]["beta_kernel"]).applyfunc(sym.Rational)
    assert flag_kernel*normalized==beta
    vec=sym.Matrix(X[1:])
    moved=change*q.subs(dict(zip(vec,change.inv()*vec)),simultaneous=True)*change.T
    ambient=X[0]**2*J
    ambient[1:,1:]+=moved
    one=[i for i,w in enumerate(model.reduced["parameter_weights"]) if w==1]
    two=[i for i,w in enumerate(model.reduced["parameter_weights"]) if w==2]
    raw=sym.Matrix(model.normalize(ambient)["parameters"])[two,:]
    assert raw[25]==-1 and raw[19]==raw[26]==0
    _,equations=specialize(model,dict(zip(one,flag_kernel*values)),True)
    for row in equations:
        assert sym.expand(sum(c*sym.prod(raw[j] for j in mon) for mon,c in row.items() if len(mon)==2))==0
    free=[]
    for i in range(7):
        at=next(j for j in range(45) if list(tangent.row(j))==list(sym.eye(7).row(i)) and not any(inverse.row(j)))
        assert at<44
        free.append(frame["chart"][at])
    normal=[raw[i] for i in free]
    e=sym.Symbol("e")
    assert normal==[-sym.Rational(8,25)*e,-sym.Rational(4,5)*e,0,0,0,0,0]
    return dict(stabilizer_column=4,finite_parameter=str(t),H12_rank=kernel.cols,
                raw_q_support={str(i):str(x) for i,x in enumerate(raw) if x},free_q_indices=free,
                seven_free_values=list(map(str,normal)),full135_t_zero_identity=True)


M2_CERTIFICATE = """
I=trim I;
C=saturate(I,ideal(x4));E=saturate(I,ideal(x1+2*x4));U=intersect(C,E);K=U/I;
Pc=ideal(x1+2*x4,x0+4*x4,x5-2*x6,2*x4^2-x6,2*x3+x6,x2+2*x6);
Pe=ideal(x4,5*x0-2*x1,x6,x5,x3,x2);
assert(radical C==Pc and radical E==Pe);
assert(not isSubset(C,E) and not isSubset(E,C));
assert(C==saturate(C,ideal(x4)) and E==saturate(E,ideal(x1)));
assert(dim(R/(C+ideal x4))==0 and numColumns basis(R/(C+ideal x4))==5);
assert(dim(R/(E+ideal x1))==0 and numColumns basis(R/(E+ideal x1))==7);
assert(apply(toList(0..3),i->hilbertFunction(i,R/(C+ideal x4)))=={1,2,2,0});
assert(apply(toList(0..3),i->hilbertFunction(i,R/(E+ideal x1)))=={1,2,4,0});
assert(dim(R/(C+E))==0 and numColumns basis(R/(C+E))==6);
assert(apply(toList(0..3),i->hilbertFunction(i,R/(C+E)))=={1,3,2,0});
assert(dim K==0 and numColumns basis K==78);
assert(apply(toList(0..7),i->hilbertFunction(i,K))=={0,0,0,10,28,37,3,0});
assert(tval%C==0 and tval%E==0 and tval%I!=0);
print("T5_PRIMARY_PASS");
"""


def curve_certificate():
    with tempfile.TemporaryDirectory(prefix="flag-t5-primary-") as path:
        temp=Path(path)
        frame,model,_,_,tangent,inverse=build(temp/"raw",CORE)
        previous=json.loads((ROOT/POINTS).read_text())
        assert frame==previous["points"][1]["frame"]
        line=identify_line(model,frame,tangent,inverse)
        environment=os.environ|{"DYLD_LIBRARY_PATH":str(ROOT/"work/m2-runtime/lib")+":/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib"}
        expressions=residuals(temp,frame,environment)
        program=temp/"verify.m2"
        program.write_text("R=QQ[x0,x1,x2,x3,x4,x5,x6,Degrees=>{1,1,2,2,1,2,2}];\nI=ideal("+",".join(expressions[:135])+");\ntval="+expressions[135]+";\n"+M2_CERTIFICATE)
        result=subprocess.run([M2,"--script",str(program)],env=environment,check=False,capture_output=True,text=True,timeout=70)
        assert result.returncode==0 and "T5_PRIMARY_PASS" in result.stdout,(result.stdout,result.stderr)
    return dict(actual_T5_curve=line,actual_flag_length=5,actual_T5_primary_length=7,
                flag_parameter="x4",T5_parameter="x1",both_parameter_weights=1,
                T5_prime="(x4,5*x0-2*x1,x6,x5,x3,x2)",T5_primary="I:(x1+2*x4)^infinity",
                flag_primary="I:x4^infinity",intersection_Hilbert=[1,3,2],intersection_length=6,
                embedded_difference_Hilbert=[0,0,0,10,28,37,3],embedded_difference_length=78,
                vertex_primary="I+(x0,...,x6)^7",t_in_both_isolated_primary_ideals=True)

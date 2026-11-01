#!/usr/bin/env python3
"""The core r=-91/400 is an actual fixed-beta T7 flag connection."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_nil22_infinity_20261101 import beta_matrix
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)
enc=lambda A:[[str(v)for v in row]for row in A.tolist()]


def run():
    start=time.time();model=S61();e=s.Symbol('e');x,y,z,w=X[1:];spatial=s.Matrix([x,y,z,w]);Q=s.zeros(4);Q[1,2]=3*e*w*w;Q[1,3]=-x*x;Q-=Q.T
    A=beta_matrix(model,Q);rr,piv=A.to_DM().to_field().rref();rr=rr.to_Matrix();L=rr[:len(piv),:]
    assert len(piv)==5 and not any(v.has(e)for v in L)
    H=s.Matrix.hstack(*L.nullspace());assert A*H==s.zeros(A.rows,11)
    free=s.symbols('a0:11');beta=H*s.Matrix(free)
    one=[i for i,wg in enumerate(model.reduced['parameter_weights'])if wg==1];two=[i for i,wg in enumerate(model.reduced['parameter_weights'])if wg==2]
    P=X[0]**2*J;P[1:,1:]+=Q;q=s.Matrix(model.normalize(P)['parameters'])[two,:]
    _,F=specialize(model,dict(zip(one,beta)),True)
    for f in F:
        value=sum(c*s.prod(q[j]for j in mon)for mon,c in f.items()if len(mon)==2)
        assert s.expand(value)==0
    assert all(not v.has(e)for v in beta)
    a=s.Symbol('a',nonzero=True);g=s.diag(a,1/a,1,1);Q1=Q.subs(e,1)
    assert(g.inv()*Q1.subs(dict(zip(spatial,g*spatial)),simultaneous=True)*g.inv().T-a**3*Q.subs(e,1/a**2)).applyfunc(s.expand)==s.zeros(4)

    old=json.loads((ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json').read_text());b=s.symbols('b0:15');Hflag=M(old['flag']['beta_kernel']);db=M(old['frame_change']['beta_derivative']);stab=M(old['frame_change']['stabilizer'])
    conditions=s.Matrix([b[1],b[2],b[3]-s.Rational(7,20)*b[12],b[13]+5*b[8]])
    sub={b[1]:0,b[2]:0,b[3]:s.Rational(7,20)*b[12],b[8]:-b[13]/5}
    assert(L*Hflag*s.Matrix(b)).subs(sub)==s.zeros(5,1)
    N=(conditions.jacobian(b)*db).subs(sub).applyfunc(s.factor)
    assert N.to_DM().rank()==2 and N[[1,3],[4,2]].det()==s.Rational(9,40)*b[12]**2
    zz,u,B=primitive_beta(model);vec=s.Matrix([B[i,j].diff(zz[k])for(i,j),k in COORDS]);basis=vec.jacobian(u);rows=list(basis.T.rref()[1]);iv=basis[rows,:].inv();action=s.Matrix.hstack(*[iv*infinitesimal_action(vec,S)[rows,:]for S in sp_basis()]);action=action.subs(dict(zip(u,Hflag*s.Matrix(b))))
    target=(Hflag*s.Matrix(b)).subs(sub);parameters=[v for v in b if v not in sub];D=target.jacobian(parameters).row_join(action.subs(sub));O=M(old['flag']['q_orbit_differential']);sample={v:j+2 for j,v in enumerate(b)}
    assert D.subs(sample).rank()==16 and D.col_join(s.zeros(44,11).row_join(O)).subs(sample).rank()==17

    # Exact finite word sends every point of the regular core into H11.
    I=b[0]*b[2]-b[3]*b[12]+b[3]**2
    core={b[1]:0,b[2]:(-s.Rational(91,400)*b[12]**2+b[3]*b[12]-b[3]**2)/b[0]}
    values=s.Matrix(b).subs(core);word=[]
    for column,parameter in [(4,(s.Rational(7,20)*b[12]-b[3])/b[0]),(2,None)]:
        if parameter is None:parameter=s.cancel(-4*(values[13]+5*values[8])/(3*values[12]))
        Dflow=db[:,column].jacobian(b);assert Dflow**4==s.zeros(15)
        flow=sum((parameter**j/s.factorial(j)*Dflow**j for j in range(4)),s.zeros(15));values=(flow*values).applyfunc(s.cancel)
        space=sum((stab[j,column]*sp_basis()[j]for j in range(10)),s.zeros(4));assert space*space==s.zeros(4)
        word.append(dict(stabilizer_column=column,parameter=str(parameter)))
    assert (L*Hflag*values).applyfunc(s.cancel)==s.zeros(5,1)
    assert s.cancel(I.subs(sub)+s.Rational(91,400)*b[12]**2)==0
    paths=[Path(__file__),ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json',ROOT/'reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json',ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_FLAG_T7_BOUNDARY_20261101',all_passed=True,potential='x^2*z+e*w^3',fixed_beta_constraints=enc(L),fixed_beta_kernel=enc(H),beta_independent_of_e=True,all135_highest_identities=True,projective_q=list(map(str,q)),nonzero_e_orbit_normalization='diag(a,a^-1,1,1), e=a^-2, projective multiplier a^3',flag_limit_constraints=list(map(str,conditions)),stabilizer_normal_matrix=enc(N),normal_rank=2,normal_minor='9*b12^2/40',beta_rank=16,joint_rank=17,core='(b0*b2-b3*b12+b3^2)/b12^2=-91/400',finite_word=word,finite_word_open='b0*b12 !=0, b1=0',limitations=['Provides actual fixed-beta reduced T7 curves through every point of the stated regular core.','No flat multiplicity4 thickening at the flag endpoint is asserted.','No t-saturation removal or entire T7 flag-boundary exhaustion outside this core chart.'],proof_note='research/notes/S61_FLAG_T7_BOUNDARY_20261101.md',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))

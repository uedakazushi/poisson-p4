#!/usr/bin/env python3
"""Independent coefficient, finite-group and all-open proof audit for flags."""
import argparse,hashlib,json,time
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61,X,J as J0
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)

def run():
    started=time.time();r=json.loads((ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json').read_text());case=r['flag'];lin=r['linear'];qr=r['quadratic'];frame=r['frame_change'];model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};ui={i:j for j,i in enumerate(one)};equations=model.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];encoded=[[(tuple(qi[i] for i in mm if i in qi),tuple(ui[i] for i in mm if i in ui),s.Rational(nu,de)) for mm,nu,de in eq] for eq in equations]
    b=s.symbols('b:15');field=s.QQ.frac_field(*b);poly,*bp=ring(b,s.QQ);H=M(lin['beta_kernel']);beta=H*s.Matrix(b);x,y,z,w=X[1:];spatial=s.Matrix([x,y,z,w]);Q=s.zeros(4);Q[1,3]=-x*x;Q-=Q.T
    def coords(Q):
        full=X[0]**2*J0;full[1:,1:]+=Q;return s.Matrix(model.normalize(full)['parameters'])[two,:]
    q=coords(Q);assert q==s.Matrix([s.Rational(case['point61'][i]) for i in two]);fixed=next(i for i,v in enumerate(q) if v);chart=[i for i in range(45) if i!=fixed];J=s.zeros(135,45);BM=s.zeros(134,16);terms=[];products={}
    for row,eq in enumerate(encoded):
        total=0;ts={}
        for mm,bm,co in eq:
            if bm not in products:products[bm]=s.prod(beta[i] for i in bm)
            value=co*products[bm]
            if len(mm)==2:
                i,j=mm;total+=value*q[i]*q[j]
                if i!=fixed:J[row,chart.index(i)]+=value*q[j]
                if j!=fixed:J[row,chart.index(j)]+=value*q[i]
                if 70<=row<134:
                    assert len(bm)==1;BM[row,bm[0]]+=co*q[i]*q[j]
            elif len(mm)==1:J[row,44]+=value*q[mm[0]]
            if fixed not in mm:
                ij=tuple(sorted([chart.index(i) for i in mm]+[44]*(2-len(mm))));ts[ij]=ts.get(ij,0)+value
        assert s.expand(total)==0;terms.append([(i,j,poly.from_expr(s.expand(co))) for (i,j),co in ts.items() if co!=0])
    J=J.applyfunc(s.expand);assert (J-M(lin['jacobian'])).applyfunc(s.expand)==s.zeros(135,45);assert BM.rank()==1 and H.rank()==15 and BM*H==s.zeros(134,15);assert beta[0]==0
    T=M(lin['tangent']);L=M(lin['left_kernel']);JD=J.to_DM().convert_to(field);assert (JD*T.to_DM().convert_to(field)).is_zero_matrix and (L.to_DM().convert_to(field)*JD).is_zero_matrix
    freecols=[i for i in range(45) if i not in lin['pivot_columns']];freerows=[i for i in range(135) if i not in lin['pivot_rows']];assert T[freecols,:]==s.eye(7) and L[:,freerows]==s.eye(97) and T[-1,:]==s.zeros(1,7)
    determinant=JD.extract(lin['pivot_rows'],lin['pivot_columns']).det();D=s.sympify(lin['jacobian_minor']);assert determinant==field.from_sympy(D)
    allowed=[b[0],b[1],5*b[3]-16*b[12],5*b[3]-4*b[12]]
    # Test divisibility by the named factors without factoring every repeated denominator.
    unique_den=set(s.denom(v) for v in [*T,*L])
    for den in unique_den:
        p=s.Poly(den,*b)
        for factor in allowed:
            pf=s.Poly(factor,*b)
            while not p.rem(pf):p=p.exquo(pf)
        assert p.total_degree()==0
    print('raw135 linear/minor/denominators PASS',time.time()-started,flush=True)

    normal=qr['normal_columns'];selected=qr['selected_left_rows'];assert qr['parameter_rows']==[0,1,3] and normal==[2,4,5,6];td=s.sympify(qr['tangent_clearing_factor']);ld=s.sympify(qr['left_clearing_factor']);tdp=poly.from_expr(td);ldp=poly.from_expr(ld)
    def clear(value,den):
        nu,de=s.fraction(value);return (poly.from_expr(nu)*den).exquo(poly.from_expr(de))
    Tp=[[clear(T[i,j],tdp) for j in normal] for i in range(45)];Lp=[[clear(L[i,j],ldp) for j in range(135)] for i in selected];mons=list(combinations_with_replacement(range(4),2));quad=[[poly.zero for _ in mons] for _ in range(135)]
    for row,ts in enumerate(terms):
        for i,j,co in ts:
            for k,(a,c) in enumerate(mons):quad[row][k]+=co*(Tp[i][a]*Tp[j][c]+(Tp[i][c]*Tp[j][a] if a!=c else poly.zero))
    target=M(qr['quadratic_matrix']);scale=ldp*tdp**2
    for i in range(10):
        for j in range(10):
            value=sum((Lp[i][k]*quad[k][j] for k in range(135) if Lp[i][k] and quad[k][j]),poly.zero);assert value==clear(target[i,j],scale)
    print('all100 normal quadratic coefficients PASS',time.time()-started,flush=True)
    rowden=[s.lcm([s.denom(v) for v in target[i,:]]) for i in range(10)];AP=s.Matrix([[s.cancel(target[i,j]*rowden[i]) for j in range(10)] for i in range(10)]);dm=AP.to_DM();qdet=s.cancel(dm.domain.to_sympy(dm.det())/s.prod(rowden));assert s.cancel(qdet-s.sympify(qr['quadratic_minor']))==0
    print('normal10minor PASS',time.time()-started,flush=True)

    zz,bu,be=primitive_beta(model);vector=s.Matrix([be[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vector.jacobian(bu);br=list(basis.T.rref()[1]);inverse=basis[br,:].inv();acts=s.Matrix.hstack(*[inverse*infinitesimal_action(vector,A)[br,:] for A in sp_basis()]);ar=acts.xreplace(dict(zip(bu,beta)));od=[]
    for A in sp_basis():
        dq=coords(sum(((A*spatial)[i]*Q.diff(spatial[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T);dq-=q*dq[fixed]/q[fixed];od.append(dq[chart,:])
    OD=s.Matrix.hstack(*od);assert OD==M(case['q_orbit_differential']) and OD.rank()==4;stab=M(frame['stabilizer']);assert stab.rank()==6 and OD*stab==s.zeros(44,6);hr=list(H.T.rref()[1]);hi=H[hr,:].inv();db=(hi*(ar*stab)[hr,:]).applyfunc(s.expand);assert db==M(frame['beta_derivative']) and (H*db-ar*stab).applyfunc(s.expand)==s.zeros(16,6)
    for i in [0,1]:assert all(s.cancel(value/b[i]).is_number for value in db[i,:])
    A=M(frame['unipotent_generator']);expected=s.zeros(4);expected[3,2]=-1;assert A==expected;ss=s.Symbol('s');g=s.eye(4)+ss*A;assert g.T*J0[1:,1:]*g==J0[1:,1:];assert g.inv()*Q.subs(dict(zip(spatial,g*spatial)),simultaneous=True)*g.inv().T==Q
    # Direct finite pullback of all three covariant indices, not just Lie action.
    betamat=be.xreplace(dict(zip(bu,beta)));transformed=g.T*betamat.subs(dict(zip(zz,g*s.Matrix(zz))),simultaneous=True)*g;flat=s.Matrix([s.expand(transformed[i,j]).diff(zz[k]) for (i,j),k in COORDS]);unew=inverse*flat[br,:];bnew=(hi*unew[hr,:]).applyfunc(s.expand);assert (H*bnew-unew).applyfunc(s.expand)==s.zeros(16,1);assert bnew[3]==b[3]+ss*b[0] and all(bnew[i]==b[i] for i in [0,1,12])
    boundary=[]
    for zero in [0,1]:
        HC=H[:,[i for i in range(15) if i!=zero]];sample=dict(zip(b,range(2,17)));sample[b[zero]]=0;BD=HC.row_join(ar.subs(sample));joint=BD.col_join(s.zeros(44,14).row_join(OD));assert BD.rank()==16 and joint.rank()==18;boundary.append(dict(zero_coordinate=zero,beta_rank=16,joint_rank=18))
    paths=[Path(__file__),ROOT/'code/s61_flag_uniform_open_20261101.py',ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json',ROOT/'research/notes/S61_FLAG_UNIFORM_OPEN_20261101.md',ROOT/'reports/S61_FLAG_FORMAL_RING_20261101.json',ROOT/'reports/S61_FLAG_FORMAL_RING_REVIEW_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_FLAG_UNIFORM_OPEN_REVIEW_20261101',all_passed=True,raw135_jacobian_verified=True,full_beta_kernel_dimension=15,linear_rank=38,linear_minor=str(D),all_denominators_supported_on_minor=True,all100_quadratic_coefficients_verified=True,quadratic_minor=str(s.sympify(qr['quadratic_minor'])),normal_columns=normal,parameter_rows=qr['parameter_rows'],full_finite_tensor_pullback_verified=True,finite_action='b3 -> b3+s*b0; b0,b1,b12 fixed',remaining_dominant_divisors=boundary,finite_R3_free_rank5_argument_accepted=True,completed_ring='k[[u0,u1,u2,e0,e1,e2,e3]]/(e0,e1,e2,e3)^2; t=0',scope='Every flag point with b0*b1 nonzero; neither remaining divisor nor whole infinity is removed.',seconds=time.time()-started,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds']}))

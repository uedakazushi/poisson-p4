#!/usr/bin/env python3
"""Independent universal-coefficient and exceptional-divisor audit for xzw infinity."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)
def run():
    start=time.time();path=ROOT/'reports/S61_NONBINARY_CANONICAL_EXHAUSTION_20261101.json';report=json.loads(path.read_text());assert report['all_passed'];un=report['universal_open'];artifact=report['artifact_divisors'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};u=s.symbols('u1:62')
    beta={i:u[i] for i in one};beta.update({3:0,4:0,11:0,12:0,49:-5*u[34],50:5*u[25]});source=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];polys=[]
    for eq in source:
        f={}
        for mon,a,b in eq:
            mm=tuple(qi[i] for i in mon if i in qi);cc=s.Rational(a,b)*s.prod(beta[i] for i in mon if i in beta);f[mm]=f.get(mm,s.S.Zero)+cc
        polys.append({mm:s.expand(cc) for mm,cc in f.items() if cc})
    fixed=qi[40];cols=[i for i in range(45) if i!=fixed];q=s.zeros(45,1);q[fixed]=-1;Jac=s.zeros(135,45)
    for row,f in enumerate(polys):
        assert sum(cc*s.prod(q[i] for i in mm) for mm,cc in f.items() if len(mm)==2)==0
        for mm,cc in f.items():
            if len(mm)==2:
                for j,i in enumerate(cols):Jac[row,j]+=cc*sum(q[mm[1-h]] for h in range(2) if mm[h]==i)
            elif len(mm)==1:Jac[row,-1]+=cc*q[mm[0]]
    assert Jac==M(un['jacobian']);T=M(un['tangent']);L=M(un['left_kernel']);assert (Jac.to_DM()*T.to_DM()).is_zero_matrix and (L.to_DM()*Jac.to_DM()).is_zero_matrix
    free=[i for i in range(45) if i not in un['pivot_columns']];leftfree=[i for i in range(135) if i not in un['pivot_rows']];assert T[free,:]==s.eye(2) and L[:,leftfree]==s.eye(92) and T[-1,:]==s.zeros(1,2)
    dmat=Jac.extract(un['pivot_rows'],un['pivot_columns']).to_DM();det=s.factor(dmat.domain.to_sympy(dmat.det()));assert s.expand(det-s.sympify(un['jacobian_minor']))==0
    denoms=set()
    for value in list(T)+list(L):
        den=s.denom(s.cancel(value))
        for fac,power in s.factor_list(den)[1]:
            assert s.rem(s.Poly(det,*[u[i] for i in one]),s.Poly(fac,*[u[i] for i in one])).is_zero
            denoms.add(str(fac))
    def quadratic(tangent,F):
        E=s.zeros(46,tangent.cols)
        for j,i in enumerate(cols):E[i,:]=tangent[j,:]
        E[45,:]=tangent[-1,:];H=s.zeros(135,3)
        for row,f in enumerate(F):
            for mm,cc in f.items():
                if fixed in mm:continue
                i,j=(*mm,*([45]*(2-len(mm))))
                H[row,0]+=cc*E[i,0]*E[j,0];H[row,1]+=cc*(E[i,0]*E[j,1]+E[i,1]*E[j,0]);H[row,2]+=cc*E[i,1]*E[j,1]
        return H
    field=s.QQ.frac_field(*[u[i] for i in one]);zero=field.zero;cache={}
    def cv(value):
        if value not in cache:cache[value]=field.from_sympy(value)
        return cache[value]
    ET=[[zero,zero] for _ in range(46)]
    for j,i in enumerate(cols):ET[i]=[cv(T[j,h]) for h in range(2)]
    ET[45]=[cv(T[-1,h]) for h in range(2)];HQ=[[zero,zero,zero] for _ in range(135)]
    for row,f in enumerate(polys):
        for mm,value in f.items():
            if fixed in mm:continue
            i,j=(*mm,*([45]*(2-len(mm))))
            if not any(ET[i]) or not any(ET[j]):continue
            cc=cv(value)
            HQ[row][0]+=cc*ET[i][0]*ET[j][0];HQ[row][1]+=cc*(ET[i][0]*ET[j][1]+ET[i][1]*ET[j][0]);HQ[row][2]+=cc*ET[i][1]*ET[j][1]
    obs=[]
    for row in L.tolist():
        val=[zero,zero,zero]
        for j,c in enumerate(row):
            if c:
                cc=cv(c)
                for h in range(3):val[h]+=cc*HQ[j][h]
        obs.append(val)
    assert obs==[[cv(s.sympify(c)) for c in row] for row in un['quadratic']]
    assert [(i,list(row)) for i,row in enumerate(obs) if any(row)]==[(40,[cv(-4),zero,zero]),(72,[cv(-2*u[48]),zero,zero])]
    # Direct fixed-Q orbit and beta action; use the explicit six coordinate constraints.
    x,y,z,w=X[1:];Y=s.Matrix([x,y,z,w]);Q=s.zeros(4);Q[1,2]=x*z;Q[1,3]=-x*w;Q-=Q.T;orbit=[]
    for HH in sp_basis():
        delta=sum(((HH*Y)[i]*Q.diff(Y[i]) for i in range(4)),s.zeros(4))-HH*Q-Q*HH.T;P=X[0]**2*J;P[1:,1:]+=delta;pp=m.normalize(P)['parameters'];dq=s.Matrix([pp[i] for i in two]);dq-=q*dq[fixed]/q[fixed];orbit.append(dq[cols,:])
    orbit=s.Matrix.hstack(*orbit);assert orbit==M(artifact['projective_Q_orbit_differential']) and orbit.rank()==7;stab=s.Matrix.hstack(*orbit.nullspace());assert stab.cols==3
    zz,bu,B=primitive_beta(m);vv=s.Matrix([B[i,j].diff(zz[k]) for (i,j),k in COORDS]);bas=vv.jacobian(bu);rr=list(bas.T.rref()[1]);left=bas[rr,:].inv();acts=s.Matrix.hstack(*[left*infinitesimal_action(vv,HH)[rr,:] for HH in sp_basis()])
    C=s.zeros(6,16)
    for row,i in enumerate((3,4,11,12)):C[row,one.index(i)]=1
    C[4,one.index(49)]=1;C[4,one.index(34)]=5;C[5,one.index(50)]=1;C[5,one.index(25)]=-5;K=s.Matrix.hstack(*C.nullspace());assert K.shape==(16,10)
    assert C*M(artifact['canonical_beta_kernel'])==s.zeros(6,10)
    cases=[]
    for record in artifact['cases']:
        condition=s.sympify(record['linear_condition']);cc=s.Matrix([condition]).jacobian(bu);HH=K*s.Matrix.hstack(*(cc*K).nullspace());assert HH.cols==9
        r=s.symbols('r0:9');beta_generic=HH*s.Matrix(r);stab_beta=(acts*stab).subs(dict(zip(bu,beta_generic)),simultaneous=True).applyfunc(s.expand);assert (C*stab_beta).applyfunc(s.expand)==s.zeros(6,3) and (cc*stab_beta).applyfunc(s.expand)==s.zeros(1,3)
        ld=record['local_data'];point=list(map(s.Rational,ld['point61']));sub=dict(zip(u,point));JP=Jac.xreplace(sub);assert JP==M(ld['jacobian']) and JP.to_DM().rank()==43
        TP=M(ld['tangent_basis']);assert JP*TP==s.zeros(135,2) and TP.rank()==2 and TP[-1,:]==s.zeros(1,2)
        LF=s.Matrix.hstack(*JP.T.nullspace()).T;FS=[{mm:cc.xreplace(sub) for mm,cc in f.items()} for f in polys];OP=LF*quadratic(TP,FS);assert OP.to_DM().rank()==1
        first=next(list(row) for row in OP.tolist() if any(row));hess=s.Matrix([[2*first[0],first[1]],[first[1],2*first[2]]]);assert hess.rank()==1
        BP=acts.subs(dict(zip(bu,[point[i] for i in one])));BD=K.row_join(BP);QD=s.zeros(44,10).row_join(orbit);assert BD.rank()==16;fiber=QD*s.Matrix.hstack(*BD.nullspace());assert fiber.rank()==1
        tangent=M(record['fixed_beta_curve_tangent']);actual=(TP*tangent)[:44,:];assert actual!=s.zeros(44,1) and fiber.row_join(actual).rank()==1 and hess*tangent==s.zeros(2,1)
        restricted=HH.row_join(BP);joint=restricted.col_join(s.zeros(44,9).row_join(orbit));assert restricted.rank()==joint.rank()==16
        reverse=list(reversed(range(19)));DR=restricted[:,reverse];piv=list(DR.to_DM().rref()[1]);minor=DR[:,piv].det(method='domain-ge');assert minor!=0
        cases.append(dict(label=record['label'],full_rank=43,quadratic_rank=1,hessian_rank=1,actual_curve_verified=True,connected_stabilizer_preserves_slice=True,alternative_beta_minor=str(minor),stratum_dimension=16))
    paths=[Path(__file__),path,ROOT/'code/s61_nonbinary_canonical_exhaustion_20261101.py',ROOT/'research/notes/S61_NONBINARY_CANONICAL_EXHAUSTION_20261101.md',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',m.source/'symplectic_slice_reduced_equations.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'research/notes/NULL_CORRELATION_DUALITY_20261101.md',ROOT/'research/notes/S61_R0_UNIFORM_VERTEX_20261101.md',ROOT/'reports/S61_NONBINARY_FORMAL_RING_20261101.json',ROOT/'reports/S61_NONBINARY_RESONANCE_REVIEW_20261101.json',ROOT/'reports/S61_NONBINARY_ZERO49_ENDPOINT_REVIEW_20261101.json',ROOT/'reports/S61_NONBINARY_RESONANCE_EXCLUSIONS_REVIEW_20261101.json']
    return dict(id='S61_NONBINARY_CANONICAL_REVIEW_20261101',all_passed=True,universal_rank=43,universal_minor=str(det),denominator_factors=sorted(denoms),universal_obstruction='-4*v0^2 and -2B*v0^2',artifact_cases=cases,accepted_entire_generic_beta_canonical_orbit_removed=True,scope='Only the canonical xzw orbit, including all its exceptional divisor strata; not other R0 normal forms or whole infinity.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'artifact_cases':len(r['artifact_cases']),'seconds':r['seconds']}))

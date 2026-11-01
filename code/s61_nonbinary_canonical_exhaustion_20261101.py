#!/usr/bin/env python3
"""Remove the whole canonical xzw nonbinary orbit over a general first jet."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_nil22_infinity_20261101 import beta_matrix
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import sp_basis,COORDS,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
def enc(M):return [[str(v) for v in row] for row in M.tolist()]
def universal():
    start=time.time();m=S61();u=s.symbols('u1:62');b={i:u[i] for i,w in enumerate(m.reduced['parameter_weights']) if w==1};b.update({3:0,4:0,11:0,12:0,49:-5*u[34],50:5*u[25]})
    two,F=specialize(m,b,True);fixed=two.index(40);cols=[j for j in range(45) if j!=fixed];J=s.zeros(135,45)
    for i,f in enumerate(F):
     for mon,c in f.items():
      if len(mon)==2:
       a,z=mon
       if a!=fixed and z==fixed:J[i,cols.index(a)]-=c
       if z!=fixed and a==fixed:J[i,cols.index(z)]-=c
      elif len(mon)==1 and mon[0]==fixed:J[i,44]-=c

    # Fraction-field RREF directly uses sparse low-degree structure.
    rr,piv=J.to_DM().rref();rr=rr.to_Matrix();assert len(piv)==43
    free=[i for i in range(45) if i not in piv];T=s.zeros(45,len(free))
    for k,j in enumerate(free):
     T[j,k]=1
     for i,c in enumerate(piv):T[c,k]=-rr[i,j]
    assert (J.to_DM()*T.to_DM()).is_zero_matrix and T[-1,:]==s.zeros(1,2)
    # Independent left projection using pivot rows; calculate only three quadratic columns.
    lr,lp=J.T.to_DM().rref();lr=lr.to_Matrix();lf=[j for j in range(135) if j not in lp];L=s.zeros(len(lf),135)
    for k,j in enumerate(lf):
     L[k,j]=1
     for i,c in enumerate(lp):L[k,c]=-lr[i,j]
    assert (L.to_DM()*J.to_DM()).is_zero_matrix
    minor=J[list(lp),list(piv)];dm=minor.to_DM();det=s.factor(dm.domain.to_sympy(dm.det()))
    Q=s.zeros(135,3);mons=[(0,0),(0,1),(1,1)]
    for r,f in enumerate(F):
     for mon,c in f.items():
      if fixed in mon:continue
      a,z=[cols.index(v) for v in mon]+[44]*(2-len(mon))
      Q[r,0]+=c*T[a,0]*T[z,0];Q[r,1]+=c*(T[a,0]*T[z,1]+T[a,1]*T[z,0]);Q[r,2]+=c*T[a,1]*T[z,1]
    A=(L.to_DM()*Q.to_DM()).to_Matrix().applyfunc(s.factor)
    assert [(i,list(A[i,:])) for i in range(A.rows) if any(A[i,:])]==[(40,[-4,0,0]),(72,[-2*u[48],0,0])]


    expected=-s.Rational(1594323,2441406250)*u[25]*u[34]*u[48]**7*(u[13]-2*u[48])*(u[13]+u[48])*(2*u[13]-u[48])**2*(5*u[13]-4*u[48])**2*(5*u[13]-3*u[48])*(5*u[13]-2*u[48])*(5*u[13]-u[48])**3*(10*u[13]+u[48])
    assert s.expand(det-expected)==0
    denominator_factors=set()
    for value in set(T).union(set(L)):
        denominator=s.fraction(s.cancel(value))[1]
        denominator_factors.update(str(s.factor(f)) for f,e in s.factor_list(denominator)[1])
    allowed={str(s.factor(f)) for f,e in s.factor_list(det)[1]}
    assert denominator_factors.issubset(allowed)
    return dict(beta={str(u[i]):str(v) for i,v in b.items()},jacobian=enc(J),tangent=enc(T),left_kernel=enc(L),pivot_rows=list(lp),pivot_columns=list(piv),jacobian_minor=str(det),quadratic=enc(A),normal_quadratic='-4*v0**2',denominator_factors=sorted(denominator_factors),full_rank=43,tangent_t_zero=True)

def artifact_strata():
    m=S61();x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);Q=s.zeros(4);Q[1,2]=x*z;Q[1,3]=-x*w;Q=Q-Q.T
    one=[i for i,a in enumerate(m.reduced['parameter_weights']) if a==1];two=[i for i,a in enumerate(m.reduced['parameter_weights']) if a==2]
    M=beta_matrix(m,Q);H=s.Matrix.hstack(*M.nullspace());assert H.cols==10
    zz,bu,B=primitive_beta(m);vector=s.Matrix([B[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vector.jacobian(bu);rows=list(basis.T.to_DM().rref()[1]);left=basis[rows,:].inv()
    acts=s.Matrix.hstack(*[(left*infinitesimal_action(vector,A)[rows,:]).applyfunc(s.expand) for A in sp_basis()])
    def coords(Q):
        P=X[0]**2*J0;P[1:,1:]+=Q;pp=m.normalize(P)['parameters'];return s.Matrix([pp[i] for i in two])
    q0=coords(Q);fixed=next(i for i,v in enumerate(q0) if v);chart=[i for i in range(45) if i!=fixed];orbit=[]
    for A in sp_basis():
        dQ=sum(((A*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T;dq=coords(dQ);dq-=q0*dq[fixed]/q0[fixed];orbit.append(dq[chart,:])
    QD=s.Matrix.hstack(*orbit);stab=s.Matrix.hstack(*QD.nullspace());assert QD.rank()==7 and stab.cols==3
    A,U,VV,Bb=[bu[one.index(i)] for i in [13,25,34,48]]
    specs=[('U0',U,{25:0}),('V0',VV,{34:0}),('AplusB',A+Bb,{13:-14}),('A2B',A-2*Bb,{13:28}),('10AplusB',10*A+Bb,{13:-s.Rational(7,5)})]
    out=[]
    for name,condition,values in specs:
        C=s.Matrix([condition]).jacobian(bu);HH=H*s.Matrix.hstack(*(C*H).nullspace());assert HH.cols==9
        b=s.symbols('b0:9');sa=(acts*stab).subs(dict(zip(bu,HH*s.Matrix(b))),simultaneous=True).applyfunc(s.expand)
        assert (C*sa).applyfunc(s.expand)==s.zeros(1,3) and (M*sa).applyfunc(s.expand)==s.zeros(134,3)
        point=s.zeros(61,1)
        for j,i in enumerate(one):point[i]=j+2
        for i in [3,4,11,12]:point[i]=0
        point[48]=14
        for i,v in values.items():point[i]=v
        point[49]=-5*point[34];point[50]=5*point[25];point[40]=-1
        assert condition.subs(dict(zip(bu,[point[i] for i in one])))==0
        local=calculate(m,True,point)
        assert (local['homogenized_jacobian_rank'],local['quadratic_span_rank'],local['common_hessian_rank'],local['homogenizer_zero_on_tangent'])==(43,1,1,True)
        sub=dict(zip(bu,[point[i] for i in one]));BD=H.row_join(acts.subs(sub));sourceQ=s.zeros(44,10).row_join(QD)
        kernel=s.Matrix.hstack(*BD.nullspace());fiber=sourceQ*kernel;fp=list(fiber.to_DM().rref()[1]);assert BD.rank()==16 and len(fp)==1
        tangent=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);td=fiber[:,fp[0]].col_join(s.zeros(1,1));jr=s.Matrix(local['jacobian']).applyfunc(s.Rational);assert jr*td==s.zeros(135,1)
        tr=list(tangent.T.to_DM().rref()[1]);tv=tangent[tr,:].inv()*td[tr,:];assert tv!=s.zeros(2,1)
        mon=local['quadratic_monomials'];quad=[sum(s.Rational(c)*tv[i]*tv[j] for c,(i,j) in zip(row,mon)) for row in local['quadratic_obstruction_basis']];assert quad==[0]
        restricted=HH.row_join(acts.subs(sub));joint=restricted.col_join(s.zeros(44,9).row_join(QD));cp=list(restricted.to_DM().rref()[1]);assert len(cp)==16 and joint.rank()==16
        dm=restricted[:,cp].to_DM();det=dm.domain.to_sympy(dm.det());assert det
        out.append(dict(label=name,linear_condition=str(condition),restricted_beta_kernel=enc(HH),stabilizer_preserves_slice=True,beta_rank=16,joint_rank=16,restricted_beta_minor_columns=cp,restricted_beta_minor=str(det),local_data=local,full_source_beta_differential=enc(BD),fixed_beta_curve_tangent=enc(tv),transverse_length_upper_bound=2,generic_curve_length_lower_bound=2,completed_ring='C[[u,e]]/(e^2)',homogenizer_zero_at_sample=True,stratum_dimension=16,whole_generic_beta_stratum_removed=True))
    return dict(canonical_beta_kernel=enc(H),projective_Q_orbit_differential=enc(QD),projective_Q_orbit_dimension=7,stabilizer=enc(stab),cases=out)
def verify():
    start=time.time();univ=universal();art=artifact_strata()
    deps=['S61_NONBINARY_FORMAL_RING','S61_NONBINARY_RESONANCE_NODES','S61_NONBINARY_ZERO49_ENDPOINT','S61_NONBINARY_RESONANCE_EXCLUSIONS']
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
    for stem in deps:
        path=ROOT/f'reports/{stem}_20261101.json';r=json.loads(path.read_text());assert r['all_passed'];paths.append(path)
    paths += [ROOT/'research/notes/S61_NONBINARY_ZERO49_ENDPOINT_20261101.md',ROOT/'research/notes/S61_R0_GENERIC_EXHAUSTION_20261101.md',ROOT/'research/notes/S61_R0_UNIFORM_VERTEX_20261101.md',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_NONBINARY_CANONICAL_EXHAUSTION_20261101',all_passed=True,field='QQ',universal_open=univ,artifact_divisors=art,conclusion='For sufficiently general beta, every infinity point in the Sp4 orbit of d_y wedge(xz d_z-xw d_w) is removed by t-saturation.',proof_note='research/notes/S61_NONBINARY_CANONICAL_EXHAUSTION_20261101.md',limitations=['Only the canonical xzw nonbinary orbit and its proved sixteen-dimensional exceptional strata.','Other nonbinary R0 normal forms z^2*w and z*w^2, their tails, flag and binary endpoints require separate analysis.','Not full infinity removal, affine finiteness, or component exhaustion.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,artifact_divisors=len(r['artifact_divisors']['cases']),seconds=r['seconds'])))

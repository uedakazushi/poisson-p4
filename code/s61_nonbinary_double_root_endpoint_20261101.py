#!/usr/bin/env python3
"""A triple nonbinary/minimal endpoint of generic multiplicity five."""
import argparse,json,time,hashlib
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X,J as J0
from s61_nil22_infinity_20261101 import beta_matrix
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
from s61_local_cubic_certificate_20261101 import profile
ROOT=Path(__file__).resolve().parents[1]
def enc(M):return [[str(v) for v in row] for row in M.tolist()]
def geometry(local):
    m=S61();x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);Q=s.zeros(4);Q[1,2]=z*z;Q[1,3]=-2*z*w;Q=Q-Q.T
    one=[i for i,a in enumerate(m.reduced['parameter_weights']) if a==1];two=[i for i,a in enumerate(m.reduced['parameter_weights']) if a==2];M=beta_matrix(m,Q);H=s.Matrix.hstack(*M.nullspace())
    zz,bu,B=primitive_beta(m);vector=s.Matrix([B[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vector.jacobian(bu);rows=list(basis.T.rref()[1]);left=basis[rows,:].inv();acts=s.Matrix.hstack(*[(left*infinitesimal_action(vector,A)[rows,:]).applyfunc(s.expand) for A in sp_basis()])
    def coords(Q):
     P=X[0]**2*J0;P[1:,1:]+=Q
     return s.Matrix([m.normalize(P)['parameters'][i] for i in two])
    q0=coords(Q);fixed=next(i for i,v in enumerate(q0) if v);chart=[i for i in range(45) if i!=fixed];orbit=[]
    for A in sp_basis():
     dQ=sum(((A*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T;dq=coords(dQ);dq-=q0*dq[fixed]/q0[fixed];orbit.append(dq[chart,:])
    QD=s.Matrix.hstack(*orbit);stab=s.Matrix.hstack(*QD.nullspace());print('orbit',QD.rank(),'stab',stab.cols)
    r=local['local'];point=s.Matrix(r['point61']).applyfunc(s.Rational);sub=dict(zip(bu,[point[i] for i in one]));BD=H.row_join(acts.subs(sub));sourceQ=s.zeros(44,10).row_join(QD);C=sourceQ*s.Matrix.hstack(*BD.nullspace());cp=list(C.to_DM().rref()[1]);assert len(cp)==1;print('beta',BD.rank(),'curve',len(cp));T=s.Matrix(r['tangent_basis']).applyfunc(s.Rational);tr=list(T.T.to_DM().rref()[1]);tv=T[tr,:].inv()*C[:,cp[0]].col_join(s.zeros(1,1))[tr,:];print('main tangent',tv)
    constraint=s.zeros(1,16);constraint[0,one.index(13)]=1;constraint[0,one.index(48)]=-2;HH=H*s.Matrix.hstack(*(constraint*H).nullspace());b=s.symbols('b:9');sa=(acts*stab).subs(dict(zip(bu,HH*s.Matrix(b))),simultaneous=True).applyfunc(s.expand);print('stabconstraint',constraint*sa);restricted=HH.row_join(acts.subs(sub));print('restricted',restricted.rank(),restricted.col_join(s.zeros(44,9).row_join(QD)).rank());cp=list(restricted.to_DM().rref()[1]);dm=restricted[:,cp].to_DM();det=dm.domain.to_sympy(dm.det());print('minor',det)
    assert BD.rank()==16 and len(cp)==16 and QD.rank()==7 and stab.cols==3
    assert (constraint*sa).applyfunc(s.expand)==s.zeros(1,3)
    assert (M*sa).applyfunc(s.expand)==s.zeros(134,3)
    assert tv==s.Matrix([-s.Rational(39,4),0,-s.Rational(1,4),0])
    return dict(beta_kernel=enc(H),restricted_beta_kernel=enc(HH),source_beta_differential=enc(BD),restricted_beta_differential=enc(restricted),q_orbit_differential=enc(QD),stabilizer=enc(stab),main_tangent=enc(tv),beta_minor_columns=cp,beta_minor=str(det),stratum_dimension=16)

def verify():
    start=time.time();point=s.zeros(61,1)
    for idx,value in {14:-s.Rational(2,3),49:-s.Rational(1,3),26:2,35:3,50:-15,51:10,25:2,27:3,34:5,36:7,37:11,48:13,46:-2}.items():point[idx-1]=value
    local=profile(point);assert local['sections']['0']['zero_dimensional'] and local['sections']['0']['length']==6
    assert local['sections']['1']['zero_dimensional'] and local['sections']['1']['length']==5
    geo=geometry(local)
    branch_path=ROOT/'reports/S61_MINIMAL_DOUBLE_ROOT_ENDPOINT_BRANCHES_20261101.json';branches=json.loads(branch_path.read_text());assert branches['all_passed']
    m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two,F=specialize(m,{i:point[i] for i in one},True)
    old=s.Matrix(branches['point61']).applyfunc(s.Rational);assert all(old[i]==point[i] for i in one) and all(old[i]==-point[i] for i in two)
    T=s.Matrix(local['local']['tangent_basis']).applyfunc(s.Rational);assert T==s.Matrix(branches['tangent_basis']).applyfunc(s.Rational)
    tau=s.Symbol('tau');tails=[];normalized_curves={}
    for label,c in branches['cases'].items():
        curve=-s.Matrix(c['q_curve']).applyfunc(s.sympify);assert curve.subs(tau,0)==s.Matrix([point[i] for i in two])
        assert all(s.expand(sum(coef*s.prod(curve[i] for i in mon) for mon,coef in f.items() if len(mon)==2))==0 for f in F)
        tangent=-s.Matrix(c['local_tangent_coordinates']).applyfunc(s.Rational)
        assert T*tangent==-s.Matrix(c['chart_tangent']).applyfunc(s.Rational)
        tails.append(tangent);normalized_curves[label]=list(map(str,curve))
    main=s.Matrix(geo['main_tangent']).applyfunc(s.Rational)
    assert s.Matrix.hstack(main,*tails).rank()==3 and all(v[0]!=0 for v in [main]+tails)
    # Exact nonhomogeneous transverse quotient (m^3=0 follows from the initial section).
    n=4;ms2=[tuple(v) for v in local['quadratic_monomials']];ms3=[tuple(v) for v in local['cubic_monomials']];ii2=[i for i,v in enumerate(ms2) if 0 not in v];ii3=[i for i,v in enumerate(ms3) if 0 not in v]
    mon=[(i,) for i in range(1,n)]+[ms2[i] for i in ii2]+[ms3[i] for i in ii3];ix={v:i for i,v in enumerate(mon)};A2=s.Matrix(local['A2']).applyfunc(s.Rational);A3=s.Matrix(local['A3']).applyfunc(s.Rational);rows=[]
    for k in range(135):
        rows.append([0]*3+list(A2[k,ii2])+list(A3[k,ii3]))
        for a in range(1,4):
            out=[s.S.Zero]*len(mon)
            for i in ii2:out[ix[tuple(sorted((a,)+ms2[i]))]]+=A2[k,i]
            if any(out):rows.append(out)
    R=s.Matrix(rows);rank=len(R.to_DM().rref()[1]);length=1+len(mon)-rank
    tq=list(map(s.Rational,local['homogenizer_quadratic']));tc=list(map(s.Rational,local['homogenizer_cubic']));t=[0]*3+[tq[i] for i in ii2]+[tc[i] for i in ii3];tnonzero=len(R.col_join(s.Matrix([t])).to_DM().rref()[1])>rank
    assert length==6 and not tnonzero
    paths=[Path(__file__),ROOT/'code/s61_local_cubic_certificate_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',branch_path,ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101.json',ROOT/'research/notes/S61_NONBINARY_ZERO49_ENDPOINT_20261101.md',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_NONBINARY_DOUBLE_ROOT_ENDPOINT_20261101',all_passed=True,field='QQ',local=local,geometry=geo,branch_generic_length_lower_bounds=[3,1,1],normalized_tail_curves=normalized_curves,branch_tangents=enc(s.Matrix.hstack(main,*tails)),generic_rank=5,transverse_length=length,finite_length_torsion_minimal_generators=1,nonflat=True,non_Cohen_Macaulay=True,transverse_t_nonzero=tnonzero,transverse_relation_rank=rank,homogenizer_nilpotent=True,whole_generic_beta_stratum_removed=True,proof_note='research/notes/S61_NONBINARY_DOUBLE_ROOT_ENDPOINT_20261101.md',limitations=['No minimum nilpotence exponent or full multiplication table.','Only the A=2B double-root16D stratum; other factors remain separate.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({k:r[k] for k in ['all_passed','generic_rank','transverse_length','transverse_t_nonzero','seconds']}))

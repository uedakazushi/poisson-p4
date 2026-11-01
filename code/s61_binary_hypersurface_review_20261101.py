#!/usr/bin/env python3
"""Parent audit of the eighteen-dimensional binary rank-drop stratum."""
import argparse,json,hashlib,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J as J0
from s61_binary_endpoints_review_20261101 import review_case
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
def verify():
    start=time.time();path=ROOT/'reports/S61_BINARY_HYPERSURFACE_RANK_DROP_20261101.json';r=json.loads(path.read_text());assert r['all_passed'];point=review_case(r['point'])
    m=S61();u=s.symbols('u1:62');one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qindex={k:i for i,k in enumerate(two)}
    extra=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json';equations=m.reduced['equations']+[json.loads(extra.read_text())['universal_new_equation']];D=s.zeros(135,45)
    for row,f in enumerate(equations):
        for mon,num,den in f:
            qs=[i for i in mon if i in qindex]
            if len(qs)!=2:continue
            b=s.prod((u[i] for i in mon if i not in qindex),start=s.S.One).subs({u[11]:0,u[12]:0});c=s.Rational(num,den)*b
            for j in range(2):
                if qs[1-j]==29:D[row,qindex[qs[j]]]+=c
    K=s.Matrix(r['constant_first70_kernel']).applyfunc(s.Rational);assert D[:70,:].rank()==16 and K.rank()==29 and D[:70,:]*K==s.zeros(70,29)
    E=(D[70:,:]*K).applyfunc(s.expand);rr=r['residual_block_rows'];cc=r['residual_block_columns'];assert sorted(sum(rr,[]))==list(range(65)) and sorted(sum(cc,[]))==list(range(29))
    for i in range(4):
        for j in range(4):
            if i>j or i==3 or j==3 or (i==0 and j==2):assert E[rr[i],cc[j]]==s.zeros(len(rr[i]),len(cc[j]))
    for i in range(3):assert E[rr[i],cc[i]]==s.Matrix(r['diagonal_blocks'][i]).applyfunc(s.sympify)
    B=E[rr[1],cc[1]].subs(u[49],-5*u[34]);kernel=s.Matrix(r['B_kernel']).applyfunc(s.sympify);assert kernel!=s.zeros(9,1) and (B*kernel).applyfunc(s.expand)==s.zeros(12,1)
    # Fresh primitive action and projective-Q stabilizer.
    z,bu,beta=primitive_beta(m);vec=s.Matrix([beta[i,j].diff(z[k]) for (i,j),k in COORDS]);basis=vec.jacobian(bu);br=list(basis.T.to_DM().rref()[1]);left=basis[br,:].inv();acts=s.Matrix.hstack(*[(left*infinitesimal_action(vec,A)[br,:]).applyfunc(s.expand) for A in sp_basis()])
    x,y,zz,w=X[1:];V=s.Matrix([x,y,zz,w]);Q=s.zeros(4);Q[1,2]=x*w;Q=Q-Q.T
    def coords(Q):
        P=X[0]**2*J0;P[1:,1:]+=Q;v=m.normalize(P)['parameters'];return s.Matrix([v[i] for i in two])
    q0=coords(Q);fixed=qindex[29];chart=[i for i in range(45) if i!=fixed];cols=[]
    for A in sp_basis():
        dQ=sum(((A*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T;dq=coords(dQ);dq-=q0*dq[fixed];cols.append(dq[chart,:])
    orbit=s.Matrix.hstack(*cols);st=s.Matrix.hstack(*orbit.nullspace());assert orbit.rank()==5 and st.cols==5
    constraint=s.Matrix([u[11],u[12],5*u[34]+u[49]]).jacobian(bu);sub={u[11]:0,u[12]:0,u[49]:-5*u[34]};assert (constraint*acts*st).subs(sub).applyfunc(s.expand)==s.zeros(3,5)
    sample={s.Symbol(k):s.Rational(v) for k,v in r['point']['sample'].items()};normal=(constraint*acts).subs(sample);reverse=normal[:,list(reversed(range(10)))];cp=list(reverse.to_DM().rref()[1]);assert len(cp)==3;minor=reverse[:,cp].det();assert minor
    S=s.Matrix(r['symmetry_matrix']).applyfunc(s.Rational);assert S.T*J0[1:,1:]*S==J0[1:,1:]
    assert (S.inv()*Q.subs(dict(zip(V,S*V)),simultaneous=True)*S.inv().T+Q).applyfunc(s.expand)==s.zeros(4)
    changed=S.T*beta.subs(dict(zip(z,S*s.Matrix(z))),simultaneous=True)*S;cv=s.Matrix([s.expand(changed[i,j]).coeff(z[k]) for (i,j),k in COORDS]);new=left*cv[br,:]
    assert s.expand(5*new[one.index(34)]+new[one.index(49)]-(5*u[13]-u[48]))==0
    paths=[Path(__file__),path,ROOT/'code/s61_binary_endpoints_review_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',extra,ROOT/'research/notes/S61_BINARY_HYPERSURFACE_RANK_DROP_20261101.md']
    return dict(id='S61_BINARY_HYPERSURFACE_REVIEW_20261101',all_passed=True,point_review=point,raw_universal_first_rank=16,upper_block_sizes=[[4,12],[12,9],[17,5]],B_kernel_verified=True,universal_residual_rank_upper=17,projective_orbit_dimension=5,projective_stabilizer_dimension=5,stabilizer_preserves_Hplus=True,alternative_beta_minor=str(minor),joint_image_dimension=18,generic_beta_fiber_dimension=2,finite_symmetry_verified=True,accepted=['Universal highest/full ranks<=33/34 on Hplus, with equality at the specified point.','Specified completed ring free65, I^4=It=t^2=0, t nonzero and quotient free64.','Relative16D smooth binary support and smooth equivariant base change give an ambient nilpotence open.'],limitations=['No whole18D-stratum generic-beta removal: its bad complement can be17D.','Quartic coefficients input from the main certificate; independent relation algebra and QQ local data reconstructed.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({k:r[k] for k in ['all_passed','alternative_beta_minor','joint_image_dimension','seconds']}))

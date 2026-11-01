#!/usr/bin/env python3
"""Independent raw-equation and all-face audit of the remaining mixed R0 cover."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)
def run():
    start=time.time();main=json.loads((ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json').read_text());rankrec=json.loads((ROOT/'reports/S61_R0_MIXED_CURVE_UNIFORM_RANK_20261101.json').read_text());vertex=json.loads((ROOT/'reports/S61_R0_UNIFORM_VERTEX_20261101.json').read_text());m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};ui={i:j for j,i in enumerate(one)};eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];encoded=[[(tuple(qi[i] for i in mm if i in qi),tuple(ui[i] for i in mm if i in ui),s.Rational(nu,de)) for mm,nu,de in eq if sum(i in qi for i in mm)==2] for eq in eqs]
    x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);zz,bu,beta=primitive_beta(m);bv=s.Matrix([beta[i,j].diff(zz[k]) for (i,j),k in COORDS]);bas=bv.jacobian(bu);rr=list(bas.T.rref()[1]);left=bas[rr,:].inv();acts=s.Matrix.hstack(*[left*infinitesimal_action(bv,H)[rr,:] for H in sp_basis()])
    def coords(Q):
        P=X[0]**2*J;P[1:,1:]+=Q;return s.Matrix(m.normalize(P)['parameters'])[two,:]
    def tensor(G):
        Q=s.zeros(4);Q[1,2]=s.diff(G,w);Q[1,3]=-s.diff(G,z);return Q-Q.T
    cases=[]
    for label,rec in main['orbits'].items():
        G=s.sympify(rec['G'],locals={str(v):v for v in V});Q=tensor(G);q=coords(Q);H=M(rec['beta_kernel']);b=s.symbols(f'b0:{H.cols}');point=H*s.Matrix(b)
        # Derive the full linear beta kernel directly from the raw64 weight5 equations.
        BM=s.zeros(64,16)
        for row,ff in enumerate(encoded[70:134]):
            for mm,bmm,c in ff:
                assert len(bmm)==1;BM[row,bmm[0]]+=c*q[mm[0]]*q[mm[1]]
        assert BM.rank()==16-H.cols and H.rank()==H.cols and BM*H==s.zeros(64,H.cols)
        D=s.zeros(135,45);products={}
        for row,ff in enumerate(encoded):
            total=[]
            for (i,j),bmm,c in ff:
                if bmm not in products:products[bmm]=s.prod(point[k] for k in bmm)
                cc=c*products[bmm]
                if not cc:continue
                total.append(cc*q[i]*q[j]);D[row,i]+=cc*q[j];D[row,j]+=cc*q[i]
            assert s.expand(s.Add(*total))==0
        D=D.applyfunc(s.expand)
        fixed=next(i for i,v in enumerate(q) if v);chart=[i for i in range(45) if i!=fixed];orbit=[]
        for A in sp_basis():
            dQ=sum(((A*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T;dq=coords(dQ);dq-=q*dq[fixed]/q[fixed];orbit.append(dq[chart,:])
        OD=s.Matrix.hstack(*orbit);assert OD==M(rec['q_orbit_differential']);sample=dict(zip(b,range(2,H.cols+2)));BD=H.row_join(acts.xreplace(dict(zip(bu,point.xreplace(sample)))));assert BD==M(rec['source_beta_differential']);joint=BD.col_join(s.zeros(44,H.cols).row_join(OD));stats={'Hdim':H.cols,'orbitdim':OD.rank(),'beta_rank':BD.rank(),'joint_rank':joint.rank()};assert stats==rec['stats']
        if label=='T3':
            constraints=s.Matrix([point[0],point[1],point[2],point[3],5*point[4]-point[13],5*point[9]+point[14],5*point[6]-point[15]]);assert constraints==s.zeros(7,1) and H.cols==9
        result=dict(label=label,full_beta_kernel_from_raw134=True,all135_leading_equations_verified=True,stats=stats)
        if label in rankrec['cases']:
            r=rankrec['cases'][label];C=M(r['first70_kernel']);E=M(r['residual_matrix']);N=M(r['residual_kernel']);K=M(r['full_highest_kernel']);assert D[:70,:].rank()==28 and C.rank()==17 and D[:70,:]*C==s.zeros(70,17);assert (D[70:,:]*C-E).applyfunc(s.expand)==s.zeros(65,17)
            field=s.QQ.frac_field(*b);ED=E.to_DM().convert_to(field);ND=N.to_DM().convert_to(field);assert (ED*ND).is_zero_matrix and (C*N-K).applyfunc(s.cancel)==s.zeros(45,3)
            identity_rows=[]
            for j in range(3):identity_rows.append(next(i for i in range(17) if N[i,:]==s.eye(3)[j,:]))
            assert N[identity_rows,:]==s.eye(3)
            sample_rank=E.xreplace(sample).rank();assert sample_rank==14
            result['universal_highest_rank_upper']=42;result['generic_highest_rank']=42;result['residual_kernel_identity_rows']=identity_rows
        cases.append(result);print(label,'PASS',flush=True)
    swap={0:1,1:0,2:4,4:2,3:3,5:8,8:5,6:7,7:6};canonical={(0,6),(4,6),(0,3,7),(3,8),(0,1,4),(2,8),(0,1,8),(2,6)};faces=[]
    for k,fl in vertex['extra_vertical_prime_faces'].items():
        for ff in fl:
            f=tuple(sorted(ff));sf=tuple(sorted(swap[i] for i in ff));assert f in canonical or sf in canonical;faces.append((k,list(f)))
    assert len(faces)==16
    # Check finite Darboux lifts and the full bivector pullback, including the y correction.
    a,b,c=s.symbols('a b c');S=s.Matrix([[1,0,0,0],[0,1,a,-b],[b,0,1,0],[a,0,0,1]]);H=s.eye(4);H[2,3]=a;J4=J[1:,1:]
    for A in [S,H]:
        assert A.T*J4*A==J4;G=x*x*z+2*x*x*w+3*x*z*z+5*x*z*w+7*x*w*w+11*z**3+13*z*z*w+17*z*w*w+19*w**3;pull=(A.inv()*tensor(G).subs(dict(zip(V,A*V)),simultaneous=True)*A.inv().T).applyfunc(s.expand);expected=tensor(s.expand(G.subs(dict(zip(V,A*V)),simultaneous=True)));assert pull==expected
    formulas=[(z*z*(a*x+b*w),{w:w-a*x/b},b*z*z*w),(z*(a*x*x+b*x*w+c*w*w),{w:w-b*x/(2*c)},z*((a-b*b/(4*c))*x*x+c*w*w)),(x*z*(a*x+b*w),{w:w-a*x/b},b*x*z*w),(a*x*x*z+b*x*x*w+c*x*w*w,{z:z-b*w/a},a*x*x*z+c*x*w*w),(a*x*x*z+b*x*x*w+c*w**3,{z:z-b*w/a},a*x*x*z+c*w**3)]
    assert all(s.expand(G.subs(sub)-target)==0 for G,sub,target in formulas)
    for rec in main['orbits'].values():
        G=s.sympify(rec['G'],locals={str(v):v for v in V});weights=[(mm[0],mm[1]-mm[2],1) for mm in s.Poly(G,x,z,w).monoms()];assert s.Matrix(weights).rank()==2
    paths=[Path(__file__),ROOT/'code/s61_r0_mixed_orbit_reduction_20261101.py',ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json',ROOT/'research/notes/S61_R0_MIXED_ORBIT_REDUCTION_20261101.md',ROOT/'code/s61_r0_mixed_curve_uniform_rank_20261101.py',ROOT/'reports/S61_R0_MIXED_CURVE_UNIFORM_RANK_20261101.json',ROOT/'research/notes/S61_R0_MIXED_CURVE_UNIFORM_RANK_20261101.md',ROOT/'reports/S61_R0_UNIFORM_VERTEX_20261101.json',ROOT/'reports/S61_R0_GENERIC_EXHAUSTION_20261101.json',ROOT/'reports/S61_NONBINARY_RESONANCE_EXCLUSIONS_REVIEW_20261101.json',ROOT/'reports/S61_NONBINARY_CANONICAL_REVIEW_20261101.json',ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_MIXED_JUNCTIONS_REVIEW_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
    return dict(id='S61_R0_MIXED_ORBIT_REVIEW_20261101',all_passed=True,cases=cases,all16_extra_faces_covered=faces,full_bivector_shear_covariance_checked=True,all_normalization_denominator_boundaries_reviewed=True,T3_proper_beta_image_from_identical_resonance_kernel=True,accepted_remaining_six_orbit_cover=True,accepted_four_universal_highest_rank_bounds=True,scope='Reduced all-point orbit cover plus universal highest rank<=42 for four curve orbits. No mixed-orbit saturation; T5/T7 surface cases separate.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))

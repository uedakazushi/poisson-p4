#!/usr/bin/env python3
"""Independent Jordan main-curve coverage, same-beta branch and cubic bound audit."""
import argparse,hashlib,json,time
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)

def length(gb,n):
    assert gb.is_zero_dimensional;leads=[p.LM(order=gb.order).exponents for p in gb.polys];counts=[]
    for degree in range(40):
        count=0
        for mon in combinations_with_replacement(range(n),degree):
            exp=tuple(mon.count(i) for i in range(n))
            if not any(all(a<=b for a,b in zip(lead,exp)) for lead in leads):count+=1
        counts.append(count)
        if count==0:return sum(counts),counts
    raise AssertionError('standard-monomial enumeration did not terminate')

def run():
    started=time.time();r=json.loads((ROOT/'reports/S61_MINIMAL_MAIN_JORDAN_REMOVAL_20261101.json').read_text());model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};ui={i:j for j,i in enumerate(one)};eqs=model.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];enc=[[(tuple(qi[i] for i in mm if i in qi),tuple(ui[i] for i in mm if i in ui),s.Rational(nu,de)) for mm,nu,de in eq]for eq in eqs];x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);a,bb=s.symbols('a b')
    def tensor(a,bb,c):
        Q=s.zeros(4);Q[0,1]=x*x;Q[1,2]=-bb*x*x-c*x*z;Q[1,3]=a*x*x+c*x*w;return Q-Q.T
    def coords(Q):
        P=X[0]**2*J;P[1:,1:]+=Q;return s.Matrix(model.normalize(P)['parameters'])[two,:]
    Q=tensor(0,1,1);q=coords(Q);H=M(r['beta_kernel']);b=s.symbols('b:9');beta=H*s.Matrix(b);fixed=next(i for i,v in enumerate(q)if v);cols=[i for i in range(45)if i!=fixed];A=s.zeros(135,45);BM=s.zeros(134,16);terms=[];products={}
    for row,eq in enumerate(enc):
        total=0;ts={}
        for mm,bm,co in eq:
            if bm not in products:products[bm]=s.prod(beta[i]for i in bm)
            val=co*products[bm]
            if len(mm)==2:
                i,j=mm;total+=val*q[i]*q[j]
                if i!=fixed:A[row,cols.index(i)]+=val*q[j]
                if j!=fixed:A[row,cols.index(j)]+=val*q[i]
                if 70<=row<134:assert len(bm)==1;BM[row,bm[0]]+=co*q[i]*q[j]
            elif len(mm)==1:A[row,44]+=val*q[mm[0]]
            if fixed not in mm:
                ij=tuple([cols.index(i)for i in mm]+[44]*(2-len(mm)));ts[ij]=ts.get(ij,0)+val
        assert s.expand(total)==0;terms.append(ts)
    assert BM.rank()==7 and H.rank()==9 and BM*H==s.zeros(134,9)
    sample=dict(zip(b,range(1,10)));A0=A.subs(sample);local=r['local_data'];assert A0==M(local['jacobian']) and A0.rank()==42 and A0[:,:44].rank()==42;T=M(local['tangent_basis']);assert A0*T==s.zeros(135,3) and T.rank()==3
    zz,bu,be=primitive_beta(model);vec=s.Matrix([be[i,j].diff(zz[k])for(i,j),k in COORDS]);basis=vec.jacobian(bu);rows=list(basis.T.rref()[1]);iv=basis[rows,:].inv();act=s.Matrix.hstack(*[iv*infinitesimal_action(vec,B)[rows,:]for B in sp_basis()]);ar=act.xreplace(dict(zip(bu,beta)));orbit=[]
    for B in sp_basis():
        dq=coords(sum(((B*V)[i]*Q.diff(V[i])for i in range(4)),s.zeros(4))-B*Q-Q*B.T);dq-=q*dq[fixed]/q[fixed];orbit.append(dq)
    O=s.Matrix.hstack(*orbit);stab=M(r['stabilizer']);assert O.rank()==7 and stab.rank()==3 and O*stab==s.zeros(45,3);actions=s.Matrix.hstack(*[M(v)for v in r['stabilizer_actions']]);assert (H*actions-ar*stab).applyfunc(s.expand)==s.zeros(16,3);D=H.row_join(ar);joint=D.col_join(s.zeros(45,9).row_join(O));assert D.subs(sample).rank()==16 and joint.subs(sample).rank()==16
    source=json.loads((ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json').read_text());pars=list(map(s.sympify,source['source_variables']));sb=M(source['source_beta']);sq=M(source['cases']['M0']['q_coordinates']);universal={s.Symbol(k):s.sympify(v)for k,v in r['universal_main_source_point'].items()};point={s.Symbol(k):s.Rational(v)for k,v in r['main_source_point'].items()};assert (sb.subs(universal)-beta).applyfunc(s.cancel)==s.zeros(16,1) and (sq.subs(universal)-q).applyfunc(s.cancel)==s.zeros(45,1)
    BD=sb.jacobian(pars).subs(point).row_join(ar.subs(sample));dQ=sq.jacobian(pars).subs(point);dQ-=q*dQ[fixed,:]/q[fixed];QD=dQ[cols,:].row_join(O[cols,:]);assert BD==M(r['main_beta_differential']) and QD==M(r['main_q_differential']);assert BD.rank()==16 and BD.col_join(QD).rank()==17;direction=M(r['main_source_direction']);tangent=M(r['main_branch_tangent']);assert BD*direction==s.zeros(16,1) and (QD*direction).col_join(s.zeros(1,1))==T*tangent and tangent[0]!=0
    # Rebuild the complete quadratic and cubic residuals from the raw135 coefficients.
    prows=list(A0.T.rref()[1]);pcols=list(A0.rref()[1]);U=s.zeros(45,42);inverse=A0[prows,pcols].inv()
    for i,col in enumerate(pcols):U[col,:]=inverse[i,:]
    mon2=list(combinations_with_replacement(range(3),2));mon3=list(combinations_with_replacement(range(3),3));indices={mon:i for i,mon in enumerate(mon3)};Q2=s.zeros(135,6);numeric=[[(i,j,s.sympify(co).subs(sample))for(i,j),co in ts.items()]for ts in terms]
    for row,ts in enumerate(numeric):
        for i,j,co in ts:
            for k,(u,v)in enumerate(mon2):Q2[row,k]+=co*(T[i,u]*T[j,v]+(T[i,v]*T[j,u]if u!=v else 0))
    W2=-U*Q2[prows,:];R2=Q2+A0*W2;Q3=s.zeros(135,10)
    for row,ts in enumerate(numeric):
        for i,j,co in ts:
            for k in range(3):
                for l,mon in enumerate(mon2):Q3[row,indices[tuple(sorted((k,)+mon))]]+=co*(T[i,k]*W2[j,l]+W2[i,l]*T[j,k])
    R3=Q3-A0*U*Q3[prows,:];rr=list(R2.T.rref()[1]);cc=list(R2.rref()[1]);pure=R3-R2[:,cc]*R2[rr,cc].inv()*R3[rr,:];cr=list(pure.T.rref()[1]);v=s.symbols('v:3');forms=[s.expand(sum(co*s.prod(v[i]for i in mon)for co,mon in zip(row,mons)))for mat,mons in[(R2[rr,:],mon2),(pure[cr,:],mon3)]for row in mat.tolist()];assert len(rr)==2 and len(cr)==7
    saved=list(map(s.sympify,r['initial_bound']['initial_forms']));assert len(forms)==len(saved)and all(s.expand(a-b)==0 for a,b in zip(forms,saved));section=s.groebner([f.subs(v[0],0)for f in forms],*v[1:],order='grevlex');localized=s.groebner([f.subs(v[0],1)for f in forms],*v[1:],order='grevlex');section_length,_=length(section,2);generic_length,counts=length(localized,2);assert generic_length==3 and counts==r['initial_bound']['standard_degrees']
    # The omitted Jordan coefficient is a true proper beta image, including all a.
    zero=r['diagonal_zero_jordan_coefficient'];zb=M(zero['beta']);nv=sorted(set().union(*(x.free_symbols for x in zb)),key=lambda x:int(str(x)[1:]));ZH=zb.jacobian(nv);zq=coords(tensor(0,0,1));ZBM=s.zeros(134,16)
    for row,eq in enumerate(enc[:134]):
        for mm,bm,co in eq:
            if len(mm)==2 and len(bm)==1:ZBM[row,bm[0]]+=co*zq[mm[0]]*zq[mm[1]]
    assert ZBM.rank()==7 and ZH.rank()==9 and ZBM*ZH==s.zeros(134,9);ZD=ZH.row_join(act.xreplace(dict(zip(bu,zb))));left=M(zero['left_kernel']);assert (ZD.T.to_DM()*left.to_DM()).is_zero_matrix;assert left.cols==1 and left.subs(dict(zip(nv,range(1,len(nv)+1)))).rank()==1;assert ZD.subs(dict(zip(nv,range(1,len(nv)+1)))).rank()==15
    C=M(r['normalization_matrix']);E=M(r['exchange_matrix']);assert (C.T*J[1:,1:]*C-J[1:,1:]).applyfunc(s.cancel)==s.zeros(4);assert (C.inv()*tensor(a,bb,1).subs(dict(zip(V,C*V)),simultaneous=True)*C.inv().T-Q).applyfunc(s.cancel)==s.zeros(4);S=s.eye(4);S[1,2]=S[3,0]=-a/2;assert (S.inv()*tensor(a,0,1).subs(dict(zip(V,S*V)),simultaneous=True)*S.inv().T-tensor(0,0,1)).applyfunc(s.expand)==s.zeros(4);assert E.T*J[1:,1:]*E==J[1:,1:];assert (E.inv()*tensor(a,bb,-1).subs(dict(zip(V,E*V)),simultaneous=True)*E.inv().T-tensor(bb,-a,1)).applyfunc(s.expand)==s.zeros(4)
    paths=[Path(__file__),ROOT/'code/s61_minimal_main_jordan_removal_20261101.py',ROOT/'reports/S61_MINIMAL_MAIN_JORDAN_REMOVAL_20261101.json',ROOT/'research/notes/S61_MINIMAL_MAIN_JORDAN_REMOVAL_20261101.md',ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json',ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_MINIMAL_MAIN_JORDAN_REVIEW_20261101',all_passed=True,all135_highest_and_sample_jacobian_reconstructed=True,beta_kernel_dimension=9,orbit_dimension=7,stratum_beta_rank=16,stratum_joint_rank=16,all_H_main_formula_recovery=True,same_beta_main_source_rank=16,same_beta_main_joint_rank=17,main_tangent=list(map(str,tangent)),raw_quadratic_rank=2,pure_cubic_rank=7,initial_section_length=section_length,initial_generic_length=generic_length,full_finite_normalizations_verified=True,zero_jordan_beta_image_dimension=15,point_t_nilpotence_and_generic16D_extension_accepted=True,scope='M0 main c=+/-1 only; no flatness, exact complete local ring, or full t^2=0 assertion.',seconds=time.time()-started,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths})
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds']}))

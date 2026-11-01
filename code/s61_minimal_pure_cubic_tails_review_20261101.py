#!/usr/bin/env python3
"""Independent all-finite coverage and raw135 rank audit of six minimal tails."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)
def run():
    start=time.time();r=json.loads((ROOT/'reports/S61_MINIMAL_PURE_CUBIC_TAILS_REMOVAL_20261101.json').read_text());m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};ui={i:j for j,i in enumerate(one)};eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];encoded=[[(tuple(qi[i] for i in mm if i in qi),tuple(ui[i] for i in mm if i in ui),s.Rational(nu,de)) for mm,nu,de in eq] for eq in eqs]
    x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);a,bb,c,d,h=s.symbols('a b c d h');b=s.symbols('b0:9');field=s.QQ.frac_field(*b)
    def tensor(kind,a,bb,c,d=0,h=1):
        Q=s.zeros(4);Q[0,1]=x*x if kind=='M0' else x*x/2;Q[2,3]=0 if kind=='M0' else x*x/2;Q[1,2]=-bb*x*x-(c+int(kind=='M1'))*x*z;Q[1,3]=a*x*x+c*x*w+2*d*x*z+3*h*z*z;return Q-Q.T
    def coords(Q):
        P=X[0]**2*J;P[1:,1:]+=Q;return s.Matrix(m.normalize(P)['parameters'])[two,:]
    for kind in ['M0','M1']:
        norm=r['normalizations'][kind];C=M(norm['matrix']);assert (C.T*J[1:,1:]*C-J[1:,1:]).applyfunc(s.cancel)==s.zeros(4);Q=tensor(kind,a,bb,c,d,h);pull=C.inv()*Q.subs(dict(zip(V,C*V)),simultaneous=True)*C.inv().T/h**2;assert (pull-tensor(kind,0,0,c)).applyfunc(s.cancel)==s.zeros(4)
        E=M(norm['paired_exchange']);paired=-c if kind=='M0' else -c-1;assert E.T*J[1:,1:]*E==J[1:,1:];Q=tensor(kind,a,bb,c,d,h);G=bb*x*x*z-a*x*x*w+paired*x*z*w+d*x*w*w-h*w**3;target=s.zeros(4);target[0,1]=x*x if kind=='M0' else x*x/2;target[2,3]=0 if kind=='M0' else x*x/2;target[1,2]=-s.diff(G,w)-int(kind=='M1')*x*z;target[1,3]=s.diff(G,z);target-=target.T;assert (E.inv()*Q.subs(dict(zip(V,E*V)),simultaneous=True)*E.inv().T-target).applyfunc(s.expand)==s.zeros(4)
    zz,bu,be=primitive_beta(m);vec=s.Matrix([be[i,j].diff(zz[k]) for (i,j),k in COORDS]);base=vec.jacobian(bu);kr=list(base.T.rref()[1]);left=base[kr,:].inv();acts=s.Matrix.hstack(*[left*infinitesimal_action(vec,A)[kr,:] for A in sp_basis()]);cases=[]
    for rec in r['cases']:
        kind=rec['kind'];c0=s.Rational(rec['c']);src=rec['universal'];H=M(src['beta_kernel']);beta=H*s.Matrix(b);q=coords(tensor(kind,0,0,c0));assert q==M(src['q_point']);fixed=next(i for i,v in enumerate(q) if v);cols=[i for i in range(45) if i!=fixed];assert cols==src['chart_columns'];A=s.zeros(135,44);tc=s.zeros(135,1);BM=s.zeros(64,16);products={}
        for row,eq in enumerate(encoded):
            total=[]
            for mm,bm,coef in eq:
                if bm not in products:products[bm]=s.prod(beta[i] for i in bm)
                v=coef*products[bm]
                if len(mm)==2:
                    i,j=mm;total.append(v*q[i]*q[j])
                    if i!=fixed:A[row,cols.index(i)]+=v*q[j]
                    if j!=fixed:A[row,cols.index(j)]+=v*q[i]
                    if 70<=row<134:
                        assert len(bm)==1;BM[row-70,bm[0]]+=coef*q[i]*q[j]
                elif len(mm)==1:tc[row]+=v*q[mm[0]]
            assert s.expand(s.Add(*total))==0
        A=A.applyfunc(s.expand);tc=tc.applyfunc(s.expand);assert A==M(src['leading_jacobian']) and tc==M(src['t_column']);assert BM.rank()==7 and H.rank()==9 and BM*H==s.zeros(64,9)
        C=M(src['first_kernel']);pr=src['first_pivot_rows'];pc=src['first_pivot_columns'];assert A[:70,:].rank()==src['first_rank'] and C.rank()==44-src['first_rank'] and A[:70,:]*C==s.zeros(70,C.cols);U=s.zeros(44,len(pc));iv=A[pr,pc].inv()
        for i,col in enumerate(pc):U[col,:]=iv[i,:]
        RR=(A*C).applyfunc(s.expand);ct=(tc-A*U*tc[pr,:]).applyfunc(s.expand);assert RR==M(src['residual_matrix']) and ct==M(src['residual_t_column']) and ct[:70,:]!=s.zeros(70,1)
        EE=RR.row_join(ct);NN=M(src['residual_kernel']);assert (EE.to_DM().convert_to(field)*NN.to_DM().convert_to(field)).is_zero_matrix;assert NN[src['free_columns'],:]==s.eye(1) and NN[-1,0]==0
        fullframe=(C*NN[:-1,:]-U*tc[pr,:]*NN[-1,:]).col_join(NN[-1,:]);assert (A.row_join(tc).to_DM()*fullframe.to_DM()).is_zero_matrix and fullframe.rank()==1
        DM=EE.extract(src['full_minor_rows'],src['full_minor_columns']).to_DM();det=s.factor(DM.domain.to_sympy(DM.det()));assert det==s.sympify(src['full_minor']) and len(src['full_minor_rows'])+src['first_rank']==44
        assert {s.factor(f) for f,e in s.factor_list(det)[1]}=={s.sympify(f['factor']) for f in rec['factors']}
        Q=tensor(kind,0,0,c0);orbit=[]
        for B in sp_basis():
            dq=coords(sum(((B*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-B*Q-Q*B.T);dq-=q*dq[fixed]/q[fixed];orbit.append(dq)
        OD=s.Matrix.hstack(*orbit);st=s.Matrix.hstack(*OD.nullspace());assert OD.rank()==8 and st.cols==2;ar=acts.xreplace(dict(zip(bu,beta)));hr=list(H.T.rref()[1]);hi=H[hr,:].inv();ba=(hi*(ar*st)[hr,:]).applyfunc(s.expand);assert (H*ba-ar*st).applyfunc(s.expand)==s.zeros(16,2);sample=dict(zip(b,range(1,10)));BD=H.row_join(ar);joint=BD.col_join(s.zeros(45,9).row_join(OD));assert BD.xreplace(sample).rank()==16 and joint.xreplace(sample).rank()==17
        factors=[]
        for f in rec['factors']:
            factor=s.sympify(f['factor']);var=s.sympify(f['variable']);val=s.sympify(f['solved_value']);assert s.degree(factor,var)==1 and s.gcd(s.diff(factor,var),factor.subs(var,0))==1 and s.cancel(factor.subs(var,val))==0
            for j in range(2):
                derivative=s.expand(sum(s.diff(factor,bb)*ba[i,j] for i,bb in enumerate(b)));quot=s.cancel(derivative/factor);assert not quot.free_symbols.intersection(b)
            subbeta=beta.subs(var,val);pars=[bb for bb in b if bb!=var];D=subbeta.jacobian(pars).row_join(ar.subs(var,val));data=dict(factor=str(factor),primitive_linear_geometrically_irreducible=True,stabilizer_preserved=True,type=f['type'])
            if f['type']=='proper_beta_image':
                N=M(f['certificate']['left_kernel']);assert (D.T.to_DM()*N.to_DM()).is_zero_matrix;assert N.subs(sample).rank()==N.cols;rank=D.subs(sample).rank();assert rank==16-N.cols==f['certificate']['image_dimension_upper'];data['universal_beta_image_dimension']=rank
            else:
                assert D.subs(sample).rank()==16;data['restricted_beta_rank']=16
                if f['type']=='smooth_curve_artifact_divisor':
                    beta0=subbeta.subs(sample);bsample=hi*beta0[hr,:];sub=dict(zip(b,bsample));ld=f['local_data'];point=list(map(s.Rational,ld['point61']));assert s.Matrix([point[i] for i in one])==beta0 and s.Matrix([point[i] for i in two])==q;J0=A.row_join(tc).subs(sub);assert J0==M(ld['jacobian']) and J0.rank()==44 and J0[:,:44].rank()==43;data['raw135_full_rank']=44
                else:assert factor==b[0]
            factors.append(data)
        cases.append(dict(kind=kind,c=str(c0),all135_universal_matrix_coefficients_verified=True,full_beta_kernel_from_raw_weight5=True,first70_rank=src['first_rank'],top70_t_compatibility_retained=True,minor_order=len(src['full_minor_rows']),minor=str(det),universal_leading_rank_upper=43,universal_full_rank_upper=44,orbit_dimension=8,joint_dimension=17,factors=factors));print(kind,c0,'PASS',flush=True)
    old=json.loads((ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json').read_text());covered=[];k=s.Symbol('k')
    for fam in old['families']:
        for tail in fam['rational_tails']:
            if tail['axis'] not in [5,8]:continue
            cc=s.Rational(tail['coordinate_formulas']['g3']);matched=[z for z in r['cases'] if z['kind']==fam['label'] and (cc==s.Rational(z['c']) if tail['axis']==5 else cc==(-s.Rational(z['c']) if fam['label']=='M0' else -s.Rational(z['c'])-1))];assert len(matched)==1
            param=s.Symbol(tail['parameter']);full=dict(tail['coordinate_formulas']);full[tail['parameter']]=str(param);assert len(full)==9
            for name,expression in full.items():assert s.cancel(s.sympify(expression).subs(param,0)-s.sympify(fam['baseline_solution'][name]).subs(k,s.Rational(tail['k'])))==0
            assert cc*(cc-1)*(cc+1)!=0 if fam['label']=='M0' else cc*(cc+1)*(2*cc+1)!=0;covered.append((fam['label'],tail['k'],tail['axis']))
    assert len(covered)==6
    paths=[Path(__file__),ROOT/'code/s61_minimal_pure_cubic_tails_removal_20261101.py',ROOT/'reports/S61_MINIMAL_PURE_CUBIC_TAILS_REMOVAL_20261101.json',ROOT/'research/notes/S61_MINIMAL_PURE_CUBIC_TAILS_REMOVAL_20261101.md',ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json',ROOT/'reports/S61_MINIMAL_JUNCTION_GENERIC_REVIEW_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_MINIMAL_PURE_CUBIC_TAILS_REVIEW_20261101',all_passed=True,full_normalizations_verified=True,cases=cases,covered_six_presentations=covered,all_finite_points_including_zero_covered=True,accepted_general_beta_removal=True,scope='Only the six stated finite g5/g8 tails. Parameter-infinity binary endpoints, mixed tails, and affine exhaustion remain separate.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))

#!/usr/bin/env python3
"""Independent universal coefficient and coverage audit of finite rank-one A1."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)

def run():
    start=time.time();path=ROOT/'reports/S61_MINIMAL_RANKONE_A1_REMOVAL_20261101.json';src=json.loads(path.read_text());assert src['all_passed'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    encoded=[[(tuple(qi[i] for i in mon if i in qi),tuple(one.index(i) for i in mon if i in one),s.Rational(nu,de)) for mon,nu,de in eq] for eq in eqs]
    def polys(beta,top=False):
        fs=[]
        products={}
        for eq in encoded:
            f={}
            for mm,bm,c in eq:
                if top and len(mm)!=2:continue
                if bm not in products:products[bm]=s.S.Zero if any(beta[i]==0 for i in bm) else s.prod(beta[i] for i in bm)
                value=products[bm]
                if value:f.setdefault(mm,[]).append(c*value)
            fs.append({mm:s.expand(s.Add(*values)) for mm,values in f.items()})
        return fs
    x,y,z,w=X[1:];Y=s.Matrix([x,y,z,w]);bs=s.symbols('b0:9');h=s.Symbol('h');sample=dict(zip(bs,range(1,10)))
    def tensor(d,hh,bb=0):
        Q=s.zeros(4);Q[0,1]=x*x;Q[1,2]=-bb*x*x-2*d*x*w-3*hh*w*w;return Q-Q.T
    def coords(Q):
        P=X[0]**2*J;P[1:,1:]+=Q;return s.Matrix(m.normalize(P)['parameters'])[two,:]
    def jac(F,q,cols):
        D=s.zeros(135,45)
        for row,f in enumerate(F):
            assert s.expand(sum(c*s.prod(q[i] for i in mm) for mm,c in f.items() if len(mm)==2))==0
            for mm,c in f.items():
                if len(mm)==2:
                    for col,i in enumerate(cols):D[row,col]+=c*sum(q[mm[1-j]] for j in range(2) if mm[j]==i)
                elif len(mm)==1:D[row,-1]+=c*q[mm[0]]
        return D
    def hess(F,T,cols):
        ET=s.zeros(46,2)
        for j,i in enumerate(cols):ET[i,:]=T[j,:]
        ET[45,:]=T[-1,:];HQ=s.zeros(135,3)
        for row,f in enumerate(F):
            for mm,c in f.items():
                i,j=(*mm,*([45]*(2-len(mm))))
                if ET[i,:]==s.zeros(1,2) or ET[j,:]==s.zeros(1,2):continue
                HQ[row,0]+=c*ET[i,0]*ET[j,0];HQ[row,1]+=c*(ET[i,0]*ET[j,1]+ET[i,1]*ET[j,0]);HQ[row,2]+=c*ET[i,1]*ET[j,1]
        return HQ
    uni=src['uniform'];BK=M(uni['beta_kernel']);q=coords(tensor(1,1));fixed=next(i for i,v in enumerate(q) if v);cols=[i for i in range(45) if i!=fixed];F=polys(BK*s.Matrix(bs));D=jac(F,q,cols);T=M(uni['kernel']);assert (D.to_DM()*T.to_DM()).is_zero_matrix and T[-1,:]==s.zeros(1,2)
    D0=D[:70,:44];pc=list(D0.rref()[1]);pr=list(D0[:,pc].T.rref()[1]);assert len(pc)==30;C=s.Matrix.hstack(*D0.nullspace());U=s.zeros(44,30);iv=D0[pr,pc].inv()
    for j,i in enumerate(pc):U[i,:]=iv[j,:]
    E=(D[:,:44]*C).row_join(D[:,-1]-D[:,:44]*U*D[pr,-1]);assert (E.to_DM()-M(uni['residual']).to_DM()).is_zero_matrix
    rp=uni['pivot_rows'];cp=uni['pivot_columns'];dm=E.extract(rp,cp).to_DM();det=s.factor(dm.domain.to_sympy(dm.det()));assert det==s.sympify(uni['jacobian_minor'])
    denoms=set()
    for value in set(T):
        for f,power in s.factor_list(s.denom(s.cancel(value)))[1]:
            assert s.rem(s.Poly(det,*bs),s.Poly(f,*bs)).is_zero;denoms.add(str(f))
    # Independence of the tangent basis on the whole minor open.
    pair=[23,43];assert T[pair,:]==s.eye(2);tminor=s.S.One
    for part in s.fraction(tminor):
        for fac,power in s.factor_list(part)[1]:assert s.rem(s.Poly(det,*bs),s.Poly(fac,*bs)).is_zero
    print('full universal rank and denominator audit passed',flush=True)
    HQ=hess(F,T,cols).applyfunc(s.cancel);YQ=(HQ-D[:,:44]*U*HQ[pr,:]).applyfunc(s.cancel);CE=E[:,cp]*E.extract(rp,cp).inv(method='DM');res=(YQ.to_DM()-CE.to_DM()*YQ[rp,:].to_DM()).to_Matrix().applyfunc(s.factor)
    difference=res.to_DM()-M(uni['full_residual']).to_DM()
    assert difference.is_zero_matrix,[(i,j,str(s.factor(v))) for i,row in enumerate(difference.to_Matrix().tolist()) for j,v in enumerate(row) if v][:4]
    nz=[i for i in range(135) if any(res[i,:])];assert res[nz,:].rank()==1 and res[nz[0],:]==s.Matrix([[0,0,-bs[7]/12]])
    print('universal quadratic obstruction passed',flush=True)
    z,bu,Beta=primitive_beta(m);vec=s.Matrix([Beta[i,j].diff(z[k]) for (i,j),k in COORDS]);basis=vec.jacobian(bu);rr=list(basis.T.rref()[1]);left=basis[rr,:].inv();actions=s.Matrix.hstack(*[left*infinitesimal_action(vec,H)[rr,:] for H in sp_basis()]);bad=[]
    specs=[({bs[3]:0},3),({bs[6]:0},6),({bs[7]:0},7),({bs[7]:2*bs[6]},7)]
    for rec,(sub,drop) in zip(src['bad_beta_images'],specs):
        beta=(BK*s.Matrix(bs)).xreplace(sub);pars=[b for i,b in enumerate(bs) if i!=drop];BD=beta.jacobian(pars).row_join(actions.xreplace(dict(zip(bu,beta))));N=M(rec['left_kernel']);assert BD==M(rec['beta_differential']) and (BD.T.to_DM()*N.to_DM()).is_zero_matrix
        free=[i for i in range(16) if any(N[i,j]==1 and all(N[i,k]==0 for k in range(N.cols) if k!=j) for j in range(N.cols))];assert N[free,:].rank()==N.cols
        rank=BD.xreplace(sample).rank();assert rank==16-N.cols==rec['image_dimension'];bad.append({'factor':rec['factor'],'universal_image_dimension':rank})
    curve_results=[]
    allbeta=s.symbols('beta0:16');generic_forms=polys(allbeta,top=True)
    for base,block in [(1,{'beta_kernel':uni['beta_kernel'],'curve':src['canonical_actual_curve']}),(0,src['pure_cubic_orbit'])]:
        K=M(block['beta_kernel']);d=base+bs[7]*(h-1)/(2*bs[6]);shape=2*bs[7]**2*(h-1)/(9*bs[6]**2);qc=coords(tensor(d,h,shape));assert qc==M(block['curve']['q']);ff=polys(K*s.Matrix(bs),top=True);scaled=[s.cancel(9*bs[6]**2*v) for v in qc]
        center=coords(tensor(base,1));bc=s.Matrix([s.expand(sum(c*s.prod(center[i] for i in mm) for mm,c in f.items() if len(mm)==2)) for f in generic_forms[:134]]).jacobian(allbeta)
        assert bc.rank()==7 and K.rank()==9 and bc*K==s.zeros(134,9)
        assert all(s.expand(sum(c*s.prod(scaled[i] for i in mm) for mm,c in f.items() if len(mm)==2))==0 for f in ff)
        assert qc.subs(h,1)==coords(tensor(base,1))
        if base==1:
            dq=qc.diff(h).subs(h,1);dq-=q*dq[fixed]/q[fixed];td=dq[cols,:].col_join(s.zeros(1,1));tv=T[pair,:].inv()*td[pair,:];assert (T*tv-td).applyfunc(s.cancel)==s.zeros(45,1) and s.cancel(tv[1])==0 and s.cancel(tv[0])!=0
        else:
            ld=block['local_data'];point=list(map(s.Rational,ld['point61']));qp=s.Matrix([point[i] for i in two]);bp=s.Matrix([point[i] for i in one]);cc=ld['chart_columns'];fi=next(i for i in range(45) if i not in cc);pf=polys(bp);DD=jac(pf,qp,cc);TT=M(ld['tangent_basis']);assert DD==M(ld['jacobian']) and DD.rank()==43 and DD[:,:44].rank()==42
            LL=s.Matrix.hstack(*DD.T.nullspace()).T;ob=LL*hess(pf,TT,cc);assert ob.rank()==1;row=next(r for r in ob.tolist() if any(r));Hess=s.Matrix([[2*row[0],row[1]],[row[1],2*row[2]]]);assert Hess.rank()==1
            dq=qc.diff(h).subs(h,1).xreplace(sample);dq-=qp*dq[fi]/qp[fi];td=dq[cc,:].col_join(s.zeros(1,1));tr=list(TT.T.rref()[1]);tv=TT[tr,:].inv()*td[tr,:];assert TT*tv==td and Hess*tv==s.zeros(2,1) and tv!=s.zeros(2,1)
            BD=K.row_join(actions.xreplace(dict(zip(bu,bp))));assert BD==M(block['beta_differential']) and BD.rank()==16
            Q=tensor(0,1);orbit=[]
            for HH in sp_basis():
                dq=coords(sum(((HH*Y)[i]*Q.diff(Y[i]) for i in range(4)),s.zeros(4))-HH*Q-Q*HH.T);dq-=qp*dq[fi]/qp[fi];orbit.append(dq[cc,:])
            OD=s.Matrix.hstack(*orbit);assert OD==M(block['projective_orbit_differential']) and OD.rank()==7
        curve_results.append({'canonical_d':base,'all135_universal_curve_identities':True,'actual_nonzero_radical_tangent_verified':True})
    paths=[Path(__file__),path,ROOT/'code/s61_minimal_rankone_A1_removal_20261101.py',ROOT/'research/notes/S61_MINIMAL_RANKONE_A1_REMOVAL_20261101.md',ROOT/'reports/S61_MINIMAL_RANKONE_JUNCTION_REVIEW_20261101.json',ROOT/'reports/S61_MINIMAL_VERTEX_REVIEW_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',m.source/'symplectic_slice_reduced_equations.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
    return dict(id='S61_MINIMAL_RANKONE_A1_REMOVAL_REVIEW_20261101',all_passed=True,uniform_full_rank=43,uniform_minor=str(det),tangent_denominator_factors=sorted(denoms),tangent_independence_rows=pair,tangent_independence_minor=str(tminor),universal_square_obstruction='-b7*v1^2/12',bad_beta_images=bad,curves=curve_results,accepted_entire_finite_A1_removed_over_general_beta=True,scope='Finite A1 parameter only; h=infinity binary endpoint and full affine exhaustion remain separate.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))

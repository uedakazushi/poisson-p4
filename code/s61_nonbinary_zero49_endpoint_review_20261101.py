#!/usr/bin/env python3
"""Independent raw135, local-cone, Artin quotient and branch audit."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X,J
from weight1_structure import COORDS,sp_basis,infinitesimal_action
from primitive_jet_elliptic_curve_20261101 import primitive_beta
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)

def standard(gb):
    assert gb.is_zero_dimensional
    lm=[p.LM(order=gb.order).exponents for p in gb.polys];result=[];counts=[]
    for d in range(20):
        now=[]
        for inds in combinations_with_replacement(range(len(gb.gens)),d):
            e=tuple(inds.count(i) for i in range(len(gb.gens)))
            if not any(all(a<=b for a,b in zip(l,e)) for l in lm):now.append(e)
        counts.append(len(now));result+=now
        if not now:return result,counts
    raise AssertionError('unbounded standard monomials')

def run():
    started=time.time();mainpath=ROOT/'reports/S61_NONBINARY_ZERO49_ENDPOINT_20261101.json'
    branchpath=ROOT/'reports/S61_MINIMAL_ZERO49_BRANCHES_20261101.json'
    report=json.loads(mainpath.read_text());branches=json.loads(branchpath.read_text());local=report['local'];ld=local['local_data']
    assert report['all_passed'] and branches['all_passed'] and ld['point61']==branches['point61']
    m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];index={i:j for j,i in enumerate(two)}
    point=list(map(s.Rational,ld['point61']));q=s.Matrix([point[i] for i in two]);cols=ld['chart_columns'];fixed=next(i for i in range(45) if i not in cols)
    source=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    def equations(beta):
        ans=[]
        for eq in source:
            f={}
            for mon,n,d in eq:
                mm=tuple(index[i] for i in mon if i in index);c=s.Rational(n,d)*s.prod(beta[one.index(i)] for i in mon if i in one)
                f[mm]=f.get(mm,s.S.Zero)+c
            ans.append({mon:s.expand(c) for mon,c in f.items() if c})
        return ans
    F=equations([point[i] for i in one]);Jac=s.zeros(135,45)
    for row,f in enumerate(F):
        assert sum(c*s.prod(q[i] for i in mon) for mon,c in f.items() if len(mon)==2)==0
        for mon,c in f.items():
            if len(mon)==2:
                for j,i in enumerate(cols):Jac[row,j]+=c*sum(q[mon[1-a]] for a in range(2) if mon[a]==i)
            elif len(mon)==1:Jac[row,-1]+=c*q[mon[0]]
    assert Jac==M(ld['jacobian']) and Jac.rank()==39 and Jac[:,:44].rank()==39
    T=M(ld['tangent_basis']);assert Jac*T==s.zeros(135,6) and T.rank()==6
    rows=local['jacobian_rows'];piv=local['jacobian_columns'];U=s.zeros(45,39);iv=Jac.extract(rows,piv).inv(method='DM')
    for a,b in enumerate(piv):U[b,:]=iv[a,:]
    assert Jac[rows,:]*U==s.eye(39)
    v=s.symbols('v0:6');mon2=list(combinations_with_replacement(range(6),2));mon3=list(combinations_with_replacement(range(6),3))
    # Rebuild the quadratic map in the full 46 homogenized coordinates.
    E=s.zeros(46,6)
    for j,i in enumerate(cols):E[i,:]=T[j,:]
    E[45,:]=T[-1,:]
    Q=s.zeros(135,21)
    for row,f in enumerate(F):
        for mon,c in f.items():
            ii,jj=(*mon,*([45]*(2-len(mon))))
            for k,(a,b) in enumerate(mon2):Q[row,k]+=c*(E[ii,a]*E[jj,b]+(E[ii,b]*E[jj,a] if a!=b else 0))
    lift2=-U*Q[rows,:];A2=Q+Jac*lift2
    assert A2==M(local['A2'])
    E2=s.zeros(46,21)
    for j,i in enumerate(cols):E2[i,:]=lift2[j,:]
    E2[45,:]=lift2[-1,:]
    Q3=s.zeros(135,56);lookup={mon:i for i,mon in enumerate(mon3)}
    for row,f in enumerate(F):
        for mon,c in f.items():
            ii,jj=(*mon,*([45]*(2-len(mon))))
            for a in range(6):
                for b,mm in enumerate(mon2):Q3[row,lookup[tuple(sorted((a,)+mm))]]+=c*(E[ii,a]*E2[jj,b]+E[jj,a]*E2[ii,b])
    lift3=-U*Q3[rows,:];A3=Q3+Jac*lift3
    assert A3==M(local['A3']) and A2.to_DM().rank()==9
    # A different row-space basis and lex order independently verify the initial cone.
    left=s.Matrix.hstack(*A2.T.nullspace()).T;pure=left*A3
    assert pure.to_DM().rank()==39
    ar=list(A2.T.to_DM().rref()[1]);br=list(pure.T.to_DM().rref()[1])
    quadrics=[sum(c*s.prod(v[i] for i in mm) for c,mm in zip(row,mon2)) for row in A2[ar,:].tolist()]
    cubics=[sum(c*s.prod(v[i] for i in mm) for c,mm in zip(row,mon3)) for row in pure[br,:].tolist()]
    gb0=s.groebner([f.subs(v[0],0) for f in quadrics],*reversed(v[1:]),order='grlex');basis0,h0=standard(gb0)
    assert h0==[1,5,6,2,0]
    gb1=s.groebner([f.subs(v[0],1) for f in quadrics+cubics],*reversed(v[1:]),order='lex');basis1,h1=standard(gb1)
    assert len(basis1)==7
    # Preserve mixed quadratic+cubic relations in the Artin quotient.
    m2=[mm for mm in mon2 if 0 not in mm];m3=[mm for mm in mon3 if 0 not in mm];mons=[(i,) for i in range(1,6)]+m2+m3;ix={mm:i for i,mm in enumerate(mons)}
    ix2=[mon2.index(mm) for mm in m2];ix3=[mon3.index(mm) for mm in m3];relations=[]
    for row in range(135):
        relations.append([s.S.Zero]*5+list(A2[row,ix2])+list(A3[row,ix3]))
        for a in range(1,6):
            rr=[s.S.Zero]*55
            for j,mm in enumerate(m2):rr[ix[tuple(sorted((a,)+mm))]]+=A2[row,ix2[j]]
            relations.append(rr)
    Rel=s.Matrix(relations);rank=Rel.to_DM().rank();assert rank==44
    tc=list(T[-1,1:])+list(lift2[-1,ix2])+list(lift3[-1,ix3]);assert list(map(str,tc))==local['transverse_t_class']
    def mult(a,b):
        out=[s.S.Zero]*55
        for i,c in enumerate(a):
            if not c:continue
            for j,d in enumerate(b):
                if d and len(mons[i])+len(mons[j])<=3:out[ix[tuple(sorted(mons[i]+mons[j]))]]+=c*d
        return out
    assert Rel.col_join(s.Matrix([tc])).to_DM().rank()==45
    assert Rel.col_join(s.Matrix([mult(tc,tc)])).to_DM().rank()==44
    tr=[tc]
    for mm in [(i,) for i in range(1,6)]+m2:
        rr=[s.S.Zero]*55;rr[ix[mm]]=1;tr.append(mult(tc,rr))
    assert Rel.col_join(s.Matrix(tr)).to_DM().rank()==46
    # Verify the actual extended minimal families directly against raw leading equations.
    A,B,K,L,UU,VV=s.symbols('A B K L U V');us=s.symbols('u1:62');bv=M(branches['source_beta']);base=list(map(s.sympify,branches['source_variables']));BF=equations(list(bv))
    sample={A:1,B:0,K:0,L:0,UU:point[25],VV:point[34],**{us[i]:point[i] for i in [24,26,33,35,36,47]}}
    betaD=M(branches['beta_differential']);assert betaD[:,:12]==bv.jacobian(base) and betaD.rank()==16
    zz,bu,beta_form=primitive_beta(m);bvect=s.Matrix([beta_form[i,j].diff(zz[k]) for (i,j),k in COORDS]);bbasis=bvect.jacobian(bu);rr=list(bbasis.T.rref()[1]);binv=bbasis[rr,:].inv()
    independent_action=s.Matrix.hstack(*[binv*infinitesimal_action(bvect,H)[rr,:] for H in sp_basis()]).subs(dict(zip(bu,[point[i] for i in one])))
    assert betaD[:,12:]==independent_action
    tangents=[M(report['geometry']['main_tangent'])]
    x,y,z,w=X[1:];Y=s.Matrix([x,y,z,w]);a,b,c=s.symbols('a b c')
    Q0=s.zeros(4);Q0[1,2]=x*z;Q0[1,3]=-x*w;Q0-=Q0.T;orbit=[]
    for HH in sp_basis():
        dQ=sum(((HH*Y)[i]*Q0.diff(Y[i]) for i in range(4)),s.zeros(4))-HH*Q0-Q0*HH.T
        P=X[0]**2*J;P[1:,1:]+=dQ;pp=m.normalize(P)['parameters'];dq=s.Matrix([pp[i] for i in two]);dq-=q*dq[fixed]/q[fixed];orbit.append(dq[cols,:])
    orbit=s.Matrix.hstack(*orbit);assert orbit.rank()==7
    H0=M(report['geometry']['beta_kernel']);BD=H0.row_join(independent_action);QD=s.zeros(44,10).row_join(orbit)
    assert BD.rank()==16 and H0.rank()==10
    r0_tangent=QD*s.Matrix.hstack(*BD.nullspace());assert r0_tangent.rank()==1
    actual_r0=(T*tangents[0])[:44,:];assert r0_tangent.row_join(actual_r0).rank()==1
    for kind,entry in branches['cases'].items():
        qr=s.Matrix(list(map(s.sympify,entry['q_regular'])));assert qr.subs(sample)==q
        assert all(s.denom(s.cancel(f)).subs(sample)!=0 for f in qr)
        assert all(s.cancel(sum(cc*s.prod(qr[i] for i in mon) for mon,cc in eq.items() if len(mon)==2))==0 for eq in BF)
        direction=M(entry['source_direction']);assert betaD*direction==s.zeros(16,1) and direction[1]==s.Rational(14164,1563)
        DQ=qr.jacobian(base).subs(sample);DQ-=q*DQ[fixed,:]/q[fixed]
        qD=M(entry['q_differential']);assert qD[:,:12]==DQ[cols,:] and qD[:,12:]==orbit
        assert betaD.col_join(qD).rank()==17 and (qD*s.Matrix.hstack(*betaD.nullspace())).rank()==1
        tangent=s.Matrix(list(map(s.Rational,entry['tangent_coordinates'])));assert T*tangent==(qD*direction).col_join(s.zeros(1,1)) and tangent[0]!=0
        tangents.append(tangent)
        H=M(entry['normalizing_shear']);assert (H*J[1:,1:]*H.T-J[1:,1:]).applyfunc(s.cancel)==s.zeros(4)
        G=a*x*x*z+b*x*x*w+c*x*z*w;Q=s.zeros(4);Q[0,1]=x*x if kind=='M0' else x*x/2;Q[2,3]=0 if kind=='M0' else x*x/2;Q[1,2]=-s.diff(G,w)-(x*z if kind=='M1' else 0);Q[1,3]=s.diff(G,z);Q-=Q.T
        assert (H.inv()*Q.subs(dict(zip(Y,H*Y)),simultaneous=True)*H.inv().T-Q.subs({a:0,b:0})).applyfunc(s.cancel)==s.zeros(4)
    assert s.Matrix.hstack(*tangents).rank()==3 and all(vv[0]!=0 for vv in tangents)
    paths=[Path(__file__),mainpath,branchpath,ROOT/'code/s61_nonbinary_zero49_endpoint_20261101.py',ROOT/'code/s61_minimal_zero49_branches_20261101.py',ROOT/'research/notes/S61_NONBINARY_ZERO49_ENDPOINT_20261101.md',ROOT/'research/notes/S61_MINIMAL_ZERO49_BRANCHES_20261101.md',ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json',ROOT/'reports/S61_NONBINARY_FORMAL_RING_20261101.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',m.source/'symplectic_slice_reduced_equations.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
    return dict(id='S61_NONBINARY_ZERO49_ENDPOINT_REVIEW_20261101',all_passed=True,field='QQ',full_jacobian_rank=39,quadratic_rank=9,pure_cubic_rank=39,transverse_quadratic_hilbert=h0,localized_cone_order='lex with reversed variables',localized_cone_standard_monomials=basis1,localized_cone_length=7,artin_relation_rank=44,transverse_length=12,transverse_infinity_length=10,transverse_t_nonzero=True,transverse_t_square_zero=True,three_independent_branch_tangents=True,source_B_derivative='14164/1563',accepted_generic_lengths=[2,3,2],accepted_torsion_minimal_generator_count=5,accepted_pointwise_t_nilpotence=True,accepted_generic16D_stratum_removal=True,limitations=['Universal rank bounds and equivariance are explicit prior dependencies.','No full t exponent, global torsion length, or removal of other endpoints is inferred.'],seconds=time.time()-started,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','localized_cone_length','transverse_length','seconds']}))

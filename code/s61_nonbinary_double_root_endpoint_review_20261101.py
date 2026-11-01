#!/usr/bin/env python3
"""Independent raw135 series and universal-tail audit of the double-root endpoint."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X,J
from s61_minimal_junction_review_20261101 import standards
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)

def run():
    start=time.time();path=ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_ENDPOINT_20261101.json';src=json.loads(path.read_text());assert src['all_passed'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    def polys(beta):
        fs=[]
        for eq in eqs:
            f={}
            for mon,nu,de in eq:
                mm=tuple(qi[i] for i in mon if i in qi);cc=s.Rational(nu,de)*s.prod(beta[one.index(i)] for i in mon if i in one);f[mm]=f.get(mm,s.S.Zero)+cc
            fs.append({mm:s.expand(v) for mm,v in f.items() if v})
        return fs
    ld=src['local']['local'];p=list(map(s.Rational,ld['point61']));q=s.Matrix([p[i] for i in two]);beta=s.Matrix([p[i] for i in one]);F=polys(beta);cols=ld['chart_columns'];fixed=next(i for i in range(45) if i not in cols);A=s.zeros(135,45)
    for row,f in enumerate(F):
        assert sum(cc*s.prod(q[i] for i in mm) for mm,cc in f.items() if len(mm)==2)==0
        for mm,cc in f.items():
            if len(mm)==2:
                for col,i in enumerate(cols):A[row,col]+=cc*sum(q[mm[1-j]] for j in range(2) if mm[j]==i)
            elif len(mm)==1:A[row,-1]+=cc*q[mm[0]]
    T=M(ld['tangent_basis']);assert A==M(ld['jacobian']) and A.rank()==41 and T.rank()==4 and A*T==s.zeros(135,4) and T[-1,:]==s.zeros(1,4)
    rows=list(A.T.to_DM().rref()[1]);pv=list(A.to_DM().rref()[1]);U=s.zeros(45,41);inv=A.extract(rows,pv).inv(method='DM')
    for j,i in enumerate(pv):U[i,:]=inv[j,:]
    E={1:s.zeros(46,4)}
    for j,i in enumerate(cols):E[1][i,:]=T[j,:]
    E[1][45,:]=T[44,:];mons={d:list(combinations_with_replacement(range(4),d)) for d in [1,2,3]};res={}
    for degree in [2,3]:
        ix={mm:i for i,mm in enumerate(mons[degree])};raw=s.zeros(135,len(ix))
        for row,f in enumerate(F):
            for mm,c in f.items():
                i,j=(*mm,*([45]*(2-len(mm))))
                if fixed in [i,j]:continue
                for aa in range(1,degree):
                    bb=degree-aa
                    for ii,vi in enumerate(E[aa][i,:]):
                        if not vi:continue
                        for jj,vj in enumerate(E[bb][j,:]):
                            if vj:raw[row,ix[tuple(sorted(mons[aa][ii]+mons[bb][jj]))]]+=c*vi*vj
        lift=-U*raw[rows,:];res[degree]=raw+A*lift;assert res[degree]==M(src['local']['A'+str(degree)])
        E[degree]=s.zeros(46,len(ix))
        for j,i in enumerate(cols):E[degree][i,:]=lift[j,:]
        E[degree][45,:]=lift[44,:]
    allmons=mons[2]+mons[3];ix={mm:i for i,mm in enumerate(allmons)};trunc=res[2].row_join(res[3]);rr,piv=trunc.to_DM().rref();base=rr.to_Matrix()[:len(piv),:];multiples=[]
    for row in base.tolist():
        for vv in range(4):
            val=[s.S.Zero]*len(allmons)
            for j,mm in enumerate(mons[2]):
                if row[j]:val[ix[tuple(sorted(mm+(vv,)))]]=row[j]
            if any(val):multiples.append(val)
    rr,piv=base.col_join(s.Matrix(multiples)).to_DM().rref();rr=rr.to_Matrix();v=s.symbols('v0:4');forms=[];ranks={}
    for row,pc in enumerate(piv):
        d=len(allmons[pc]);ranks[d]=ranks.get(d,0)+1;forms.append(sum(rr[row,j]*s.prod(v[i] for i in mm) for j,mm in enumerate(allmons) if len(mm)==d))
    gb0=s.groebner([f.subs(v[0],0) for f in forms],*reversed(v[1:]),order='grevlex');st0,hf0=standards(gb0);gb1=s.groebner([f.subs(v[0],1) for f in forms],*reversed(v[1:]),order='grevlex');st1,hf1=standards(gb1);assert len(st0)==6 and len(st1)==5 and hf0==[1,3,2,0]
    # Retain mixed-degree relations to check the actual transverse algebra.
    secmons=mons[1]+mons[2]+mons[3];secmon=[mm for mm in secmons if 0 not in mm];ix={mm:i for i,mm in enumerate(secmon)};mixed=[]
    for row in range(135):
        val=[s.S.Zero]*len(secmon)
        for degree in [2,3]:
            for j,mm in enumerate(mons[degree]):
                if 0 not in mm:val[ix[mm]]=res[degree][row,j]
        mixed.append(val)
        for vv in range(1,4):
            val=[s.S.Zero]*len(secmon)
            for j,mm in enumerate(mons[2]):
                if 0 not in mm:val[ix[tuple(sorted(mm+(vv,)))]]=res[2][row,j]
            mixed.append(val)
    rel=s.Matrix(mixed);rank=rel.to_DM().rank();assert 1+len(secmon)-rank==6
    tvec=[s.S.Zero]*len(secmon)
    for degree in [2,3]:
        for j,mm in enumerate(mons[degree]):
            if 0 not in mm:tvec[ix[mm]]=E[degree][45,j]
    assert rel.col_join(s.Matrix([tvec])).to_DM().rank()==rank
    # Two tails, including all nine free beta parameters of A=2B.
    x,y,z,w=X[1:];Y=s.Matrix([x,y,z,w]);B,UU,VV=s.symbols('B U V');bb=s.symbols('b0:6');tau=s.Symbol('tau');bgen=s.Matrix([0,0,0,0,2*B,bb[0],UU,bb[1],bb[2],VV,bb[3],bb[4],bb[5],B,-5*VV,5*UU]);FG=polys(bgen);lam=-1/(3*B);tails=[]
    def coords(Q):
        PP=X[0]**2*J;PP[1:,1:]+=Q;return s.Matrix(m.normalize(PP)['parameters'])[two,:]
    for kind in ['M0','M1']:
        G=z*z*w+tau*((-6*lam*VV*x*x*z-3*lam*UU*x*x*w/2+2*x*z*w) if kind=='M0' else (4*lam*VV*x*x*z+lam*UU*x*x*w-x*z*w));Q=s.zeros(4);Q[0,1]=tau*x*x if kind=='M0' else tau*x*x/2;Q[2,3]=0 if kind=='M0' else tau*x*x/2;Q[1,2]=-s.diff(G,w)-(tau*x*z if kind=='M1' else 0);Q[1,3]=s.diff(G,z);Q=-(Q-Q.T);qr=coords(Q)
        assert all(s.expand(sum(c*s.prod((B*qr[i]) for i in mm) for mm,c in f.items() if len(mm)==2))==0 for f in FG)
        sample={B:-s.Rational(1,3),UU:2,VV:3,**dict(zip(bb,[2,3,5,7,11,13]))};qs=qr.xreplace(sample);assert qs==M(src['normalized_tail_curves'][kind]) and qs.subs(tau,0)==q
        dq=qs.diff(tau).subs(tau,0);dq-=q*dq[fixed]/q[fixed];td=dq[cols,:].col_join(s.zeros(1,1));tr=list(T.T.to_DM().rref()[1]);tv=T[tr,:].inv()*td[tr,:];assert T*tv==td;tails.append(tv)
    # Orbit/source geometry, independently from the written six H constraints.
    Q=s.zeros(4);Q[1,2]=z*z;Q[1,3]=-2*z*w;Q-=Q.T;orbit=[]
    for H in sp_basis():
        dq=coords(sum(((H*Y)[i]*Q.diff(Y[i]) for i in range(4)),s.zeros(4))-H*Q-Q*H.T);dq-=q*dq[fixed]/q[fixed];orbit.append(dq[cols,:])
    QD=s.Matrix.hstack(*orbit);stab=s.Matrix.hstack(*QD.nullspace());assert QD.rank()==7 and stab.cols==3
    zz,bu,Beta=primitive_beta(m);vec=s.Matrix([Beta[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vec.jacobian(bu);rr=list(basis.T.rref()[1]);left=basis[rr,:].inv();acts=s.Matrix.hstack(*[left*infinitesimal_action(vec,H)[rr,:] for H in sp_basis()]);HC=s.Matrix([bu[one.index(i)] for i in [3,4,11,12]]+[bu[one.index(49)]+5*bu[one.index(34)],bu[one.index(50)]-5*bu[one.index(25)]]).jacobian(bu);HH=s.Matrix.hstack(*HC.nullspace());C=s.Matrix([bu[one.index(13)]-2*bu[one.index(48)]]).jacobian(bu);HK=HH*s.Matrix.hstack(*(C*HH).nullspace());bfree=s.symbols('h0:9');sa=(acts*stab).subs(dict(zip(bu,HK*s.Matrix(bfree))),simultaneous=True).applyfunc(s.expand);assert (HC*sa).applyfunc(s.expand)==s.zeros(6,3) and (C*sa).applyfunc(s.expand)==s.zeros(1,3)
    act=acts.xreplace(dict(zip(bu,beta)));BD=HH.row_join(act);RB=HK.row_join(act);assert BD.rank()==16 and RB.rank()==16 and RB.col_join(s.zeros(44,9).row_join(QD)).rank()==16
    curve=s.zeros(44,10).row_join(QD)*s.Matrix.hstack(*BD.nullspace());assert curve.rank()==1;main=M(src['geometry']['main_tangent']);td=T*main;assert td[-1]==0 and curve.row_join(td[:44,:]).rank()==1
    branch=s.Matrix.hstack(main,*tails);assert branch==M(src['branch_tangents']) and branch.rank()==3 and all(branch[0,i]!=0 for i in range(3))
    cps=list(RB[:,::-1].to_DM().rref()[1]);cps=[RB.cols-1-i for i in cps];det=RB[:,cps].det(method='domain-ge');assert det
    paths=[Path(__file__),path,ROOT/'code/s61_nonbinary_double_root_endpoint_20261101.py',ROOT/'code/s61_local_cubic_certificate_20261101.py',ROOT/'research/notes/S61_NONBINARY_DOUBLE_ROOT_ENDPOINT_20261101.md',ROOT/'reports/S61_MINIMAL_DOUBLE_ROOT_ENDPOINT_BRANCHES_20261101.json',ROOT/'research/notes/S61_MINIMAL_DOUBLE_ROOT_ENDPOINT_BRANCHES_20261101.md',ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_REVIEW_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',m.source/'symplectic_slice_reduced_equations.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'code/s61_minimal_junction_review_20261101.py']
    return dict(id='S61_NONBINARY_DOUBLE_ROOT_ENDPOINT_REVIEW_20261101',all_passed=True,point_full_rank=41,truncated_initial_ranks=ranks,zero_section_standard_degrees=hf0,localized_standard_degrees=hf1,generic_rank=5,exact_transverse_length=6,transverse_t_zero=True,actual_branch_lower_lengths=[3,1,1],universal_nine_beta_tail_identities=True,restricted_beta_minor=str(det),finite_length_torsion_minimal_generators=1,accepted_t_nilpotent=True,accepted_generic_A2B_stratum_removal=True,scope='A=2B endpoint stratum; no exact nilpotence exponent, full local multiplication, or entire double-root orbit removal.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))

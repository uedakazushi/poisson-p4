#!/usr/bin/env python3
"""Independent raw135 quartic, mixed-branch, and pure-double-root coverage audit."""
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
    start=time.time();path=ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_MIXED_JUNCTIONS_20261101.json';src=json.loads(path.read_text());assert src['all_passed'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];encoded=[[(tuple(qi[i] for i in mon if i in qi),tuple(one.index(i) for i in mon if i in one),s.Rational(nu,de)) for mon,nu,de in eq] for eq in eqs]
    def polys(beta,top=False):
        out=[];products={}
        for eq in encoded:
            f={}
            for mm,bm,c in eq:
                if top and len(mm)!=2:continue
                if bm not in products:products[bm]=s.S.Zero if any(beta[i]==0 for i in bm) else s.prod(beta[i] for i in bm)
                if products[bm]:f.setdefault(mm,[]).append(c*products[bm])
            out.append({mm:s.expand(s.Add(*values)) for mm,values in f.items()})
        return out
    x,y,z,w=X[1:];Y=s.Matrix([x,y,z,w]);tau=s.Symbol('tau');alpha=s.Symbol('alpha',nonzero=True)
    def tensor(G):
        Q=s.zeros(4);Q[1,2]=s.diff(G,w);Q[1,3]=-s.diff(G,z);return Q-Q.T
    def coords(Q):
        PP=X[0]**2*J;PP[1:,1:]+=Q;return s.Matrix(m.normalize(PP)['parameters'])[two,:]
    Q=tensor(z*z*w);q0=coords(Q);fixed=next(i for i,v in enumerate(q0) if v);cols=[i for i in range(45) if i!=fixed];orbit=[]
    for H in sp_basis():
        dq=coords(sum(((H*Y)[i]*Q.diff(Y[i]) for i in range(4)),s.zeros(4))-H*Q-Q*H.T);dq-=q0*dq[fixed]/q0[fixed];orbit.append(dq[cols,:])
    QD=s.Matrix.hstack(*orbit);stab=s.Matrix.hstack(*QD.nullspace());assert QD.rank()==7 and stab.cols==3
    zz,bu,Beta=primitive_beta(m);vec=s.Matrix([Beta[i,j].diff(zz[k]) for (i,j),k in COORDS]);bas=vec.jacobian(bu);rr=list(bas.T.rref()[1]);left=bas[rr,:].inv();acts=s.Matrix.hstack(*[left*infinitesimal_action(vec,H)[rr,:] for H in sp_basis()]);HC=s.Matrix([bu[one.index(i)] for i in [3,4,11,12]]+[bu[one.index(49)]+5*bu[one.index(34)],bu[one.index(50)]-5*bu[one.index(25)]]).jacobian(bu);HH=s.Matrix.hstack(*HC.nullspace());torus_action=left*infinitesimal_action(vec,s.diag(1,-1,0,0))[rr,:]
    u=s.symbols('u1:62');out=[]
    for label,rec in src['cases'].items():
        timer=time.time();ld=rec['cubic']['local'];p=list(map(s.Rational,ld['point61']));beta=s.Matrix([p[i] for i in one]);q=s.Matrix([p[i] for i in two]);assert q==q0 and ld['chart_columns']==cols;F=polys(beta);A=s.zeros(135,45)
        for row,f in enumerate(F):
            assert sum(c*s.prod(q[i] for i in mm) for mm,c in f.items() if len(mm)==2)==0
            for mm,c in f.items():
                if len(mm)==2:
                    for col,i in enumerate(cols):A[row,col]+=c*sum(q[mm[1-j]] for j in range(2) if mm[j]==i)
                elif len(mm)==1:A[row,-1]+=c*q[mm[0]]
        T=M(ld['tangent_basis']);n=T.cols;assert A==M(ld['jacobian']) and A*T==s.zeros(135,n) and T.rank()==n
        rows=list(A.T.to_DM().rref()[1]);piv=list(A.to_DM().rref()[1]);U=s.zeros(45,len(rows));iv=A.extract(rows,piv).inv(method='DM')
        for j,i in enumerate(piv):U[i,:]=iv[j,:]
        E={1:s.zeros(46,n)}
        for j,i in enumerate(cols):E[1][i,:]=T[j,:]
        E[1][45,:]=T[44,:];mons={d:list(combinations_with_replacement(range(n),d)) for d in range(1,5)};res={}
        for degree in [2,3,4]:
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
            lift=-U*raw[rows,:];res[degree]=raw+A*lift;assert res[degree]==M((rec['cubic'] if degree<4 else rec)['A'+str(degree)]);E[degree]=s.zeros(46,len(ix))
            for j,i in enumerate(cols):E[degree][i,:]=lift[j,:]
            E[degree][45,:]=lift[44,:]
        # Full truncated ideal, with all degree<=2 multipliers, then one RREF.
        allmons=mons[2]+mons[3]+mons[4];ix={mm:i for i,mm in enumerate(allmons)};rr,pv=s.Matrix.hstack(res[2],res[3],res[4]).to_DM().rref();base=rr.to_Matrix()[:len(pv),:];multiples=[]
        for row in base.tolist():
            for fac in mons[1]+mons[2]:
                val=[s.S.Zero]*len(allmons)
                for j,mm in enumerate(allmons):
                    if row[j] and len(mm)+len(fac)<=4:val[ix[tuple(sorted(mm+fac))]]=row[j]
                if any(val):multiples.append(val)
        rr,pv=base.col_join(s.Matrix(multiples)).to_DM().rref();rr=rr.to_Matrix();v=s.symbols('v0:'+str(n));forms=[];initial_ranks={}
        for row,pc in enumerate(pv):
            d=len(allmons[pc]);initial_ranks[d]=initial_ranks.get(d,0)+1;forms.append(sum(rr[row,j]*s.prod(v[i] for i in mm) for j,mm in enumerate(allmons) if len(mm)==d))
        counts=[]
        for value in [0,1]:
            gb=s.groebner([f.subs(v[0],value-sum(v[1:])) for f in forms],*reversed(v[1:]),order='grevlex');st,hf=standards(gb);assert len(st)==5;counts.append(hf)
        # Actual mixed curve and its complete nine-variable leading-rank bound.
        extra=x*x*z if label=='B0' else x*w*w;Qmix=tensor(z*z*w+tau*extra);qc=coords(Qmix);assert qc==M(src['branches'][label]['q_curve']) and qc.subs(tau,0)==q
        bv={i:u[i] for i in one};bv.update({3:0,4:0,11:0,12:0,49:-5*u[34],50:5*u[25]});bv.update({48:0} if label=='B0' else {13:3*u[48]/5});FG=polys(s.Matrix([bv[i] for i in one]),top=True)
        assert all(s.expand(sum(c*s.prod(qc[i] for i in mm) for mm,c in f.items()))==0 for f in FG)
        qm=qc.subs(tau,1);D=s.zeros(135,45)
        for row,f in enumerate(FG):
            for (i,j),c in f.items():D[row,i]+=c*qm[j];D[row,j]+=c*qm[i]
        known=src['universal_mixed_rank'][label];K=M(known['kernel_first70']);assert K.rank()==17 and D[:70,:].rank()==28 and D[:70,:]*K==s.zeros(70,17);R=D[70:,:]*K;assert R==M(known['residual']);assert len(R.to_DM().rref()[1])==14
        S=s.diag(alpha,1/alpha,1,1);power=2 if label=='B0' else 1;pulled=S.inv()*Qmix.subs(dict(zip(Y,S*Y)),simultaneous=True)*S.inv().T/alpha;assert (pulled-tensor(z*z*w+tau*alpha**power*extra)).applyfunc(s.expand)==s.zeros(4)
        cond=bu[one.index(48)] if label=='B0' else 5*bu[one.index(13)]-3*bu[one.index(48)];C=s.Matrix([cond]).jacobian(bu);HK=HH*s.Matrix.hstack(*(C*HH).nullspace());b=s.symbols('b0:9');arg=dict(zip(bu,HK*s.Matrix(b)));st=(acts*stab).subs(arg,simultaneous=True).applyfunc(s.expand);assert (HC*st).applyfunc(s.expand)==s.zeros(6,3) and (C*st).applyfunc(s.expand)==s.zeros(1,3)
        da=torus_action.subs(arg,simultaneous=True).applyfunc(s.expand);assert (HC*da).applyfunc(s.expand)==s.zeros(6,1) and (C*da).applyfunc(s.expand)==s.zeros(1,1)
        act=acts.xreplace(dict(zip(bu,beta)));BD=HH.row_join(act);RB=HK.row_join(act);assert BD.rank()==16 and RB.rank()==16 and RB.col_join(s.zeros(44,9).row_join(QD)).rank()==16
        curve=s.zeros(44,10).row_join(QD)*s.Matrix.hstack(*BD.nullspace());assert curve.rank()==1;main=M(src['geometry'][label]['main_tangent']);mainq=T*main;assert mainq[-1]==0 and curve.row_join(mainq[:44,:]).rank()==1
        dq=qc.diff(tau);dq-=q*dq[fixed]/q[fixed];td=dq[cols,:].col_join(s.zeros(1,1));tr=list(T.T.rref()[1]);tail=T[tr,:].inv()*td[tr,:];assert T*tail==td and tail==M(src['branches'][label]['tangent']);assert s.Matrix.hstack(main,tail).rank()==2 and sum(main)!=0 and sum(tail)!=0
        for tangent in [main,tail]:assert all(s.expand(f.subs(dict(zip(v,tangent))))==0 for f in forms)
        out.append(dict(label=label,raw135_quartic_series_verified=True,initial_ranks=initial_ranks,zero_section_standard_degrees=counts[0],localized_standard_degrees=counts[1],universal_mixed_leading_rank=42,universal_nine_beta_curve=True,finite_torus_normalization=True,source_beta_rank=16,restricted_joint_rank=16,actual_branch_tangents=[list(map(str,main)),list(map(str,tail))],accepted_free_rank=5,accepted_t_zero=True,accepted_square_zero_nilradical=True,seconds=time.time()-timer));print(json.dumps(out[-1]),flush=True)
    # Independent squarefree coverage of every factor of the universal minor.
    op=json.loads((ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101.json').read_text());det=s.sympify(op['jacobian_minor']);AA,BB,UU,VV=u[13],u[48],u[25],u[34];factors=[UU,VV,BB,AA-2*BB,AA-BB,AA+BB,2*AA-BB,5*AA-7*BB,5*AA-4*BB,5*AA-3*BB,5*AA-BB,20*AA**2-59*AA*BB+65*BB**2,11*AA**3-39*AA**2*BB+51*AA*BB**2-7*BB**3]
    actual=[f for f,e in s.factor_list(det)[1]];assert len(actual)==13
    for f in actual:assert any(s.cancel(f/g).is_number and s.cancel(f/g)!=0 for g in factors)
    deps=['S61_NONBINARY_DOUBLE_ROOT_OPEN','S61_NONBINARY_DOUBLE_ROOT_ARTIFACTS','S61_NONBINARY_DOUBLE_ROOT_ENDPOINT','S61_NONBINARY_DOUBLE_ROOT_TORSION_DIVISORS','S61_NONBINARY_RESONANCE_EXCLUSIONS','S61_NONBINARY_DOUBLE_ROOT_EXHAUSTION'];paths=[Path(__file__),path,ROOT/'code/s61_nonbinary_double_root_mixed_junctions_20261101.py',ROOT/'code/s61_local_quartic_certificate_20261101.py',ROOT/'research/notes/S61_NONBINARY_DOUBLE_ROOT_MIXED_JUNCTIONS_20261101.md',ROOT/'research/notes/S61_NONBINARY_DOUBLE_ROOT_EXHAUSTION_20261101.md',ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_REVIEW_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',m.source/'symplectic_slice_reduced_equations.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'code/s61_minimal_junction_review_20261101.py']
    for stem in deps:
        pth=ROOT/f'reports/{stem}_20261101.json';assert json.loads(pth.read_text())['all_passed'];paths.append(pth)
    return dict(id='S61_NONBINARY_DOUBLE_ROOT_MIXED_JUNCTIONS_REVIEW_20261101',all_passed=True,cases=out,all13_QQ_minor_factors_covered=True,accepted_all_pure_double_root_points_removed_over_general_beta=True,scope='Two mixed junction strata and assembly for the pure z^2*w / z*w^2 orbits; not all mixed tails or full infinity/affine exhaustion.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))

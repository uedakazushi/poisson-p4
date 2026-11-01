#!/usr/bin/env python3
"""Identify all remaining minimal flag endpoints with four dominant16D strata."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_minimal_rankone_A1_removal_20261101 import beta_actions
from s61_minimal_baseline_thickenings_20261101 import coords,enc
from s61_minimal_infinity_local_data_20261101 import calculate
from weight1_structure import COORDS,sp_basis,minor_witness
ROOT=Path(__file__).resolve().parents[1]
def run():
    start=time.time();m=S61();source=json.loads((ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json').read_text());zz,u,B=primitive_beta(m);vec=s.Matrix([B[i,j].diff(zz[h])for(i,j),h in COORDS]);basis=vec.jacobian(u);rows=list(basis.T.to_DM().rref()[1]);iv=basis[rows,:].inv();u,action=beta_actions(m);two=[i for i,wg in enumerate(m.reduced['parameter_weights'])if wg==2];one=[i for i,wg in enumerate(m.reduced['parameter_weights'])if wg==1];x,y,z,w=X[1:];xx=s.Matrix([x,y,z,w]);Q=s.zeros(4);Q[1,3]=-x*x;Q-=Q.T;q0=s.Matrix(coords(m,Q))[two,:];fixed=next(i for i,v in enumerate(q0)if v);orbit=[]
    for A in sp_basis():
        D=sum(((A*xx)[i]*Q.diff(xx[i])for i in range(4)),s.zeros(4))-A*Q-Q*A.T;dq=s.Matrix(coords(m,D))[two,:];dq-=q0*dq[fixed]/q0[fixed];orbit.append(dq)
    O=s.Matrix.hstack(*orbit);stab=s.Matrix.hstack(*O.nullspace());assert O.rank()==4 and stab.cols==6
    e,K,L,U,V,k=s.symbols('e K L U V k');n=s.symbols('n0:6');params=[e,K,L,U,V]+list(n);base=s.Matrix([0,0,0,0,k,n[0],U,n[1],n[2],V,n[3],n[4],n[5],2*k+1,K-5*V,5*U-L]);strata={}
    for label,kk in [('free_M0',-s.Rational(1,3)),('free_M1',-s.Rational(4,9)),('main_pole',-s.Rational(7,9))]:
        beta=e*base.subs(k,kk);# use a constant linear frame for its identical image H11
        beta=beta.subs({a:a/e for a in [K,L,U,V]+list(n)},simultaneous=True).applyfunc(s.cancel);H=beta.jacobian(params);strata[label]=(beta,H,params)
    c=s.symbols('c0:10');beta=s.Matrix([0,0,c[0],c[1],4*c[2],0,c[3],0,c[4],c[5],c[6],c[7],c[8],5*c[2],c[9],3*c[3]]);strata['pure_cubic']=(beta,beta.jacobian(c),list(c));out={}
    for label,(beta,H,pars)in strata.items():
        AB=action.subs(dict(zip(u,beta)));ann=s.Matrix.hstack(*H.T.nullspace()).T;N=(ann*AB*stab).applyfunc(s.factor);nr=N.to_DM().rank();assert nr==(2 if label=='pure_cubic'else 1);D=H.row_join(AB);joint=D.col_join(s.zeros(45,H.cols).row_join(O));sample=dict(zip(pars,range(1,len(pars)+1)));bw=minor_witness(D.subs(sample));jw=minor_witness(joint.subs(sample));assert bw['rank']==jw['rank']==16 and 4+H.cols+nr==16
        point=s.Matrix(coords(m,Q))
        for i,vv in zip(one,beta.subs(sample)):point[i]=vv
        local=calculate(m,True,point);out[label]=dict(beta=list(map(str,beta)),parameter_order=list(map(str,pars)),beta_kernel=enc(H),annihilator=enc(ann),stabilizer_normal=enc(N),stabilizer_normal_rank=nr,beta_minor=bw,joint_minor=jw,closure_dimension=16,sample={str(a):str(v)for a,v in sample.items()},local_data=local)
    def pullbeta(T,kv):
        BB=B.subs(dict(zip(u,base.subs(k,kv))));C=T.T*BB.subs(dict(zip(zz,T*s.Matrix(zz))),simultaneous=True)*T;return(iv*s.Matrix([C[i,j].diff(zz[h])for(i,j),h in COORDS])[rows,:]).applyfunc(s.factor)
    def pullq(T,P):return(T.inv()*P.subs(dict(zip(xx,T*xx)),simultaneous=True)*T.inv().T).applyfunc(s.factor)
    E=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,-1],[0,0,1,0]]);R=s.Matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])
    def pureT(a):
        S=s.eye(4);S[1,3]=-a;S[2,0]=a;return S*R
    g=s.symbols('g0:9');mon=[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3];records=[]
    for family in source['families']:
        label=family['label']
        for row in family['baseline_free_lines']:
            axis=row['axis'];kv=s.Rational(row['k']);G=mon[axis];P=s.zeros(4);P[1,2]=-s.diff(G,w);P[1,3]=s.diff(G,z);P-=P.T;T=s.eye(4)if axis==0 else E;PB=pullbeta(T,kv);tr=pullq(T,P);assert tr==-Q;target='free_'+label;ann=s.Matrix(out[target]['annihilator']);assert(ann*PB).applyfunc(s.factor)==s.zeros(ann.rows,1);records.append(dict(kind=label,source='free',axis=axis,k=str(kv),G_infinity=str(G),matrix=enc(T),projective_factor='-1',normalized_beta=list(map(str,PB)),target=target))
        for kv,axis in [(-s.Rational(7,9),0),(-s.Rational(2,9),1)]:
            G=sum(s.sympify(family['baseline_solution'].get(str(a),str(a)))*mon[j]for j,a in enumerate(g));res=s.factor(s.cancel((k-kv)*G).subs(k,kv));assert s.factor(res/mon[axis]).free_symbols.isdisjoint({x,y,z,w});P=s.zeros(4);P[1,2]=-s.diff(res,w);P[1,3]=s.diff(res,z);P-=P.T;T=s.eye(4)if axis==0 else E;tr=pullq(T,P);factor=s.factor(tr[1,3]/Q[1,3]);assert(tr-factor*Q).applyfunc(s.factor)==s.zeros(4);PB=pullbeta(T,kv);ann=s.Matrix(out['main_pole']['annihilator']);assert(ann*PB).applyfunc(s.factor)==s.zeros(ann.rows,1);records.append(dict(kind=label,source='main_pole',axis=axis,k=str(kv),residue_G=str(res),matrix=enc(T),projective_factor=str(factor),normalized_beta=list(map(str,PB)),target='main_pole'))
        for row in family['rational_tails']:
            if not((row['axis']==5 and row['k']=='-1')or(row['axis']==8 and row['k']=='0')):continue
            axis=row['axis'];kv=s.Rational(row['k']);par=g[axis];G=sum(s.sympify(row['coordinate_formulas'].get(str(a),str(a)))*mon[j]for j,a in enumerate(g));Gi=s.diff(G,par);assert not Gi.has(par);P=s.zeros(4);P[1,2]=-s.diff(Gi,w);P[1,3]=s.diff(Gi,z);P-=P.T;T=pureT(L)if axis==5 else E*pureT(-K);tr=pullq(T,P);factor=s.factor(tr[1,3]/Q[1,3]);assert(tr-factor*Q).applyfunc(s.factor)==s.zeros(4);PB=pullbeta(T,kv);ann=s.Matrix(out['pure_cubic']['annihilator']);assert(ann*PB).applyfunc(s.factor)==s.zeros(ann.rows,1);records.append(dict(kind=label,source='pure_cubic',axis=axis,k=str(kv),G_infinity=str(Gi),matrix=enc(T),projective_factor=str(factor),normalized_beta=list(map(str,PB)),target='pure_cubic'))
    assert len(records)==12
    # The source pure-cubic image is dense in the fullH10 used above.
    PB=e*pullbeta(pureT(L),-1);H=s.Matrix(out['pure_cubic']['beta_kernel']);DD=PB.jacobian(params);num=DD.subs(dict(zip(params,range(1,12))));pw=minor_witness(num);assert pw['rank']==10
    for r in records:
        T=s.Matrix(r['matrix']);assert(T*J[1:,1:]*T.T-J[1:,1:]).applyfunc(s.factor)==s.zeros(4)
    deps=['code/s61_model.py','code/primitive_jet_elliptic_curve_20261101.py','code/s61_minimal_rankone_A1_removal_20261101.py','code/s61_minimal_baseline_thickenings_20261101.py','code/s61_minimal_infinity_local_data_20261101.py','code/weight1_structure.py','reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json','reports/S61_UNIVERSAL_WEIGHT6_20261031.json','reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json'];paths=[Path(__file__)]+[ROOT/f for f in deps]
    return dict(id='S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101',all_passed=True,flag_q=list(map(str,q0)),projective_orbit_dimension=4,stabilizer=enc(stab),strata=out,presentations=records,pure_source_dense_minor=pw,conclusion='The remaining four free-line, four pure-cubic, and four main-pole presentations lie in four dominant16D flag strata.',scope=['A finite endpoint identification and dimension certificate only; no t-nilpotence or saturation removal claimed.','All three free/main strata have u5=u12=u13=0; the four pure flags have u5=0 and generic u12!=0.','The b0!=0,b1=0 invariant-r chart does not cover these strata.'],proof_note='research/notes/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.md',input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in paths},seconds=time.time()-start)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run()
    with a.output.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))

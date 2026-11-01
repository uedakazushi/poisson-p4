#!/usr/bin/env python3
"""Two genuine minimal-modular branches through a zero-modular infinity point."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_nil22_infinity_20261101 import beta_matrix
from s61_minimal_infinity_local_data_20261101 import calculate
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action,minor_witness
ROOT=Path(__file__).resolve().parents[1]
def enc(M):return [[str(v) for v in row] for row in M.tolist()]
def run():
    start=time.time();m=S61();x,y,z,w=X[1:];V4=s.Matrix([x,y,z,w]);one=[i for i,a in enumerate(m.reduced['parameter_weights']) if a==1];two=[i for i,a in enumerate(m.reduced['parameter_weights']) if a==2];u=s.symbols('u1:62');A,B,K,L,U,V=s.symbols('A B K L U V');a,b,c=s.symbols('a b c');G=a*x*x*z+b*x*x*w+c*x*z*w
    bsub={u[3]:0,u[4]:0,u[11]:0,u[12]:0,u[13]:A,u[48]:B,u[25]:U,u[34]:V,u[49]:K-5*V,u[50]:5*U-L};base=[A,B,K,L,U,V]+[u[i] for i in [24,26,33,35,36,47]];bv=s.Matrix([u[i].subs(bsub) for i in one]);p=s.zeros(61,1)
    for j,i in enumerate(one):p[i]=j+2
    for i in [3,4,11,12]:p[i]=0
    p[13]=1;p[48]=0;p[49]=-5*p[34];p[50]=5*p[25];p[40]=-1
    sample={A:1,B:0,K:0,L:0,U:p[25],V:p[34],**{u[i]:p[i] for i in [24,26,33,35,36,47]}};local=calculate(m,True,p);assert local['homogenized_jacobian_rank']==39
    zz,bu,Beta=primitive_beta(m);vec=s.Matrix([Beta[i,j].diff(zz[k]) for (i,j),k in COORDS]);bas=vec.jacobian(bu);rows=list(bas.T.rref()[1]);left=bas[rows,:].inv();BA=s.Matrix.hstack(*[left*infinitesimal_action(vec,H)[rows,:] for H in sp_basis()]).subs(dict(zip(bu,[p[i] for i in one])));BD=bv.jacobian(base).row_join(BA);bm=minor_witness(BD);assert bm['rank']==16
    extra=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json';G6=json.loads(extra.read_text())['universal_new_equation'];G6top=[term for term in G6 if sum(i in two for i in term[0])==2];old_rank=json.loads((ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json').read_text());out={}
    for kind in ['M0','M1']:
        Q=s.zeros(4);Q[0,1]=x*x if kind=='M0' else x*x/2;Q[2,3]=0 if kind=='M0' else x*x/2;Q[1,2]=-s.diff(G,w)-(x*z if kind=='M1' else 0);Q[1,3]=s.diff(G,z);Q-=Q.T
        sa=-a/(c+1) if kind=='M0' else -a/c;sb=-b/(c-1) if kind=='M0' else -b/(c+1);shear=s.Matrix([[1,0,0,0],[0,1,sa,-sb],[sb,0,1,0],[sa,0,0,1]])
        assert (shear*J[1:,1:]*shear.T-J[1:,1:]).applyfunc(s.cancel)==s.zeros(4)
        assert (shear.inv()*Q.subs(dict(zip(V4,shear*V4)),simultaneous=True)*shear.inv().T-Q.subs({a:0,b:0})).applyfunc(s.cancel)==s.zeros(4)
        canonical_c=c if kind=='M0' else 2*c+1;oldQ=s.zeros(4);oldQ[0,1]=x*x;oldQ[2,3]=0 if kind=='M0' else x*x;oldQ[1,2]=-canonical_c*x*z-(x*z if kind=='M1' else 0);oldQ[1,3]=canonical_c*x*w-(x*w if kind=='M1' else 0);oldQ-=oldQ.T
        assert (oldQ-(1 if kind=='M0' else 2)*Q.subs({a:0,b:0})).applyfunc(s.expand)==s.zeros(4);assert old_rank['cases'][kind]['rank_certificates']['135']['full_rank']==(42 if kind=='M0' else 43)
        F=[s.factor(f) for f in beta_matrix(m,Q)*bv];assert [i for i,f in enumerate(F) if f!=0]==[70,84,87];cs=s.solve(F[87],c)[0];aa=s.solve(F[84].subs(c,cs),a)[0];bb=s.solve(F[70].subs(c,cs),b)[0];sol={a:s.factor(aa),b:s.factor(bb),c:s.factor(cs)};assert all(s.cancel(f.subs(sol))==0 for f in F)
        P=X[0]**2*J;P[1:,1:]+=Q;raw=m.normalize(P)['parameters'];values=[bv[one.index(i)] if i in one else raw[i] for i in range(61)];assert s.expand(m.evaluate(G6top,values))==0
        Qreg=(B*Q.subs(sol)).applyfunc(s.cancel);assert all(s.denom(v).subs(sample)!=0 for v in Qreg);Qend=Qreg.subs(sample);P=X[0]**2*J;P[1:,1:]+=Qreg;pars=m.normalize(P)['parameters'];qr=s.Matrix([s.cancel(pars[i]) for i in two]);qe=qr.subs(sample);fixed=two.index(local['chart_fixed_u']-1);target=s.Matrix([p[i] for i in two]);scale=target[fixed]/qe[fixed];qr*=scale;qe*=scale;assert qe==target
        cols=local['chart_columns'];DQ=qr.jacobian(base).subs(sample);DQ-=qe*(DQ[fixed,:]/qe[fixed]);qd=[]
        for H in sp_basis():
            dQ=sum(((H*V4)[i]*Qend.diff(V4[i]) for i in range(4)),s.zeros(4))-H*Qend-Qend*H.T;PP=X[0]**2*J;PP[1:,1:]+=scale*dQ;qq=m.normalize(PP)['parameters'];dq=s.Matrix([qq[i] for i in two]);dq-=qe*dq[fixed]/qe[fixed];qd.append(dq[cols,:])
        QD=DQ[cols,:].row_join(s.Matrix.hstack(*qd));N=s.Matrix.hstack(*BD.nullspace());C=QD*N;joint=minor_witness(BD.col_join(QD));assert C.rank()==1 and joint['rank']==17
        j=C.rref()[1][0];direction=N[:,j];assert direction[1]!=0;tv=C[:,j].col_join(s.zeros(1,1));T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);rr=list(T.T.rref()[1]);v=T[rr,:].inv()*tv[rr,:];assert T*v==tv and s.Matrix(local['jacobian']).applyfunc(s.Rational)*tv==s.zeros(135,1);assert v[0]!=0
        out[kind]=dict(solution={str(k):str(v) for k,v in sol.items()},regular_Q=enc(Qreg),projective_scale=str(scale),q_regular=list(map(str,qr)),q_differential=enc(QD),joint_minor=joint,fixed_beta_q_rank=1,source_direction=list(map(str,direction)),B_derivative=str(direction[1]),tangent_coordinates=list(map(str,v)),normalizing_shear=enc(shear),old_canonical_parameter=str(canonical_c),old_canonical_scale=1 if kind=='M0' else 2,generic_full_rank_upper_bound=42 if kind=='M0' else 43,generic_full_length_lower_bound=3 if kind=='M0' else 2)
    assert s.Matrix.hstack(*[s.Matrix(v['tangent_coordinates']).applyfunc(s.Rational) for v in out.values()]).rank()==2
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',extra,ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_MINIMAL_ZERO49_BRANCHES_20261101',all_passed=True,field='QQ',point61=list(map(str,p)),local_data=local,source_beta=list(map(str,bv)),source_variables=list(map(str,base)),beta_differential=enc(BD),beta_minor=bm,cases=out,scope='Two actual smooth fixed-beta infinity branches through the specified zero-modular point, with distinct tangent lines and generic full lengths at least3 and2. No branch exhaustion or completed-ring determination.',proof_note='research/notes/S61_MINIMAL_ZERO49_BRANCHES_20261101.md',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
if __name__=='__main__':main()

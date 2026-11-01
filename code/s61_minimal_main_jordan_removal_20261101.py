#!/usr/bin/env python3
"""Remove the M0 c=+/-1 main-curve Jordan stratum over general beta."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_minimal_baseline_thickenings_20261101 import tensor,coords,enc
from s61_nil22_infinity_20261101 import beta_matrix
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_rankone_A1_removal_20261101 import beta_actions
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_minimal_semisimple_tails_removal_20261101 import minor_witness
from s61_minimal_mixed_tails_removal_20261101 import proper_image
from s61_minimal_local_initial_bound_20261101 import cubic_bound
from weight1_structure import sp_basis
ROOT=Path(__file__).resolve().parents[1]
def orbit(m,Q,two,fi):
 q=s.Matrix(coords(m,Q))[two,:];V=s.Matrix(X[1:]);ds=[]
 for H in sp_basis():
  dQ=sum(((H*V)[i]*Q.diff(V[i])for i in range(4)),s.zeros(4))-H*Q-Q*H.T;dq=s.Matrix(coords(m,dQ))[two,:];dq-=q*dq[fi]/q[fi];ds.append(dq)
 return s.Matrix.hstack(*ds)
def run():
 st=time.time();m=S61();Q=tensor('M0',0,1,1);K=s.Matrix.hstack(*beta_matrix(m,Q).nullspace());assert K.shape==(16,9);b=s.symbols('b0:9');u,act=beta_actions(m);bv=K*s.Matrix(b);ar=act.subs(dict(zip(u,bv)));D=K.row_join(ar);two=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==2];one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1];q0=s.Matrix(coords(m,Q))[two,:];fi=next(i for i,v in enumerate(q0)if v);O=orbit(m,Q,two,fi);stab=O.nullspace();assert len(stab)==3 and O.to_DM().rank()==7;kr=list(K.T.to_DM().rref()[1]);actions=[]
 for ss in stab:
  db=ar*ss;v=K[kr,:].inv()*db[kr,:];assert (K*v-db).applyfunc(s.expand)==s.zeros(16,1);actions.append(list(map(str,v)))
 _,F=specialize(m,dict(zip(one,bv)),True);assert all(s.expand(sum(v*s.prod(q0[i]for i in mon)for mon,v in f.items()if len(mon)==2))==0 for f in F)
 sample=dict(zip(b,range(1,10)));point=s.Matrix(coords(m,Q))
 for i,v in zip(one,bv.subs(sample)):point[i]=v
 loc=calculate(m,True,point);assert (loc['leading_jacobian_rank'],loc['homogenized_jacobian_rank'],loc['quadratic_span_rank'],loc['common_hessian_rank'])==(42,42,2,2)
 source_path=ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json';source=json.loads(source_path.read_text());pars=list(map(s.sympify,source['source_variables']));sourcebeta=s.Matrix(source['source_beta']).applyfunc(s.sympify);q=s.Matrix(source['cases']['M0']['q_coordinates']).applyfunc(s.sympify);beta=point[one,:];universal_sourcepoint=s.solve(list(sourcebeta-bv),pars,dict=True)[0];assert (q.subs(universal_sourcepoint)-q0).applyfunc(s.cancel)==s.zeros(45,1);sourcepoint=s.solve(list(sourcebeta-beta),pars,dict=True)[0];assert q.subs(sourcepoint)==q0
 BD=sourcebeta.jacobian(pars).subs(sourcepoint).row_join(act.subs(dict(zip(u,beta))));co=loc['chart_columns'];assert fi==next(i for i in range(45)if i not in co);Dq=q.jacobian(pars).subs(sourcepoint);Dq-=q0*Dq[fi,:]/q0[fi];QD=Dq[co,:].row_join(O[co,:]);N=s.Matrix.hstack(*BD.nullspace());C=QD*N;assert C.to_DM().rank()==1;j=next(j for j in range(C.cols)if C[:,j]!=s.zeros(44,1));tv=C[:,j].col_join(s.zeros(1,1));T=s.Matrix(loc['tangent_basis']).applyfunc(s.Rational);rr=list(T.T.to_DM().rref()[1]);v=T[rr,:].inv()*tv[rr,:];assert T*v==tv;branches=[dict(kind='M0_main',lower=3,tangent=list(map(str,v)))];bound=cubic_bound(m,loc,branches,3)
 strata_beta=minor_witness(D.subs(sample));strata_joint=minor_witness(D.col_join(s.zeros(45,K.cols).row_join(O)).subs(sample));branch_beta=minor_witness(BD);branch_joint=minor_witness(BD.col_join(QD));assert strata_beta['rank']==strata_joint['rank']==16 and branch_beta['rank']==16 and branch_joint['rank']==17
 # The zero Jordan coefficient is not silently omitted.
 Qd=tensor('M0',0,0,1);Kd=s.Matrix.hstack(*beta_matrix(m,Qd).nullspace());zero=proper_image(Kd,act,u);assert zero['image_dimension_upper']==15
 a,bb=s.symbols('a b');V=s.Matrix(X[1:]);S=s.eye(4);S[1,2]=-a/2;S[3,0]=-a/2;R=s.diag(1,1,bb,1/bb);Cmat=S*R;assert (Cmat*J[1:,1:]*Cmat.T-J[1:,1:]).applyfunc(s.cancel)==s.zeros(4);pull=Cmat.inv()*tensor('M0',a,bb,1).subs(dict(zip(V,Cmat*V)),simultaneous=True)*Cmat.inv().T;assert (pull-Q).applyfunc(s.cancel)==s.zeros(4)
 E=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,-1],[0,0,1,0]]);assert E*J[1:,1:]*E.T==J[1:,1:];pull=E.inv()*tensor('M0',a,bb,-1).subs(dict(zip(V,E*V)),simultaneous=True)*E.inv().T;assert (pull-tensor('M0',bb,-a,1)).applyfunc(s.expand)==s.zeros(4)
 paths=[Path(__file__),source_path,ROOT/'code/s61_model.py',ROOT/'code/s61_minimal_baseline_thickenings_20261101.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_minimal_rankone_A1_removal_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_minimal_semisimple_tails_removal_20261101.py',ROOT/'code/s61_minimal_mixed_tails_removal_20261101.py',ROOT/'code/s61_minimal_local_initial_bound_20261101.py',ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
 return dict(id='S61_MINIMAL_MAIN_JORDAN_REMOVAL_20261101',all_passed=True,field='QQ',all135_highest_identities=True,universal_main_source_point={str(k):str(v)for k,v in universal_sourcepoint.items()},beta=list(map(str,bv)),beta_kernel=enc(K),orbit_dimension=7,stabilizer=enc(s.Matrix.hstack(*stab)),stabilizer_actions=actions,stratum_beta_minor=strata_beta,stratum_joint_minor=strata_joint,local_data=loc,main_source_point={str(k):str(x)for k,x in sourcepoint.items()},main_beta_differential=enc(BD),main_q_differential=enc(QD),main_beta_minor=branch_beta,main_joint_minor=branch_joint,main_source_direction=list(map(str,N[:,j])),main_branch_tangent=list(map(str,v)),initial_bound=bound,diagonal_zero_jordan_coefficient=zero,normalization_matrix=enc(Cmat),normalization_denominator='b',exchange_matrix=enc(E),generic_rank=3,t_nilpotent=True,conclusion='The whole M0 main c=+/-1 locus, including zero Jordan coefficient, is removed over sufficiently general beta.',limitations=['No flatness or exact full completed ring at the Jordan sample is claimed.','The rest of the main c-parameter family is not removed here.','No full generic affine-fiber exhaustion.'],seconds=time.time()-st,proof_note='research/notes/S61_MINIMAL_MAIN_JORDAN_REMOVAL_20261101.md',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths})
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run()
 with a.output.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
if __name__=='__main__':main()

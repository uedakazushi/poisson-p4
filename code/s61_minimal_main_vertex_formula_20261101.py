#!/usr/bin/env python3
"""Nondiagonal minimal main curves and the rank-one M1 point stratum."""
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
 start=time.time();m=S61();x,y,z,w=X[1:];vv=s.Matrix([x,y,z,w]);one=[i for i,d in enumerate(m.reduced['parameter_weights']) if d==1];two=[i for i,d in enumerate(m.reduced['parameter_weights']) if d==2];u=s.symbols('u1:62');A,B,P,R,U,V,K,L=s.symbols('A B P R U V K L');a,b=s.symbols('a b');nuisance=[u[i] for i in [24,26,33,35,36,47]];params=[A,B,P,R,U,V,K,L]+nuisance;bs={u[3]:0,u[11]:0,u[4]:P,u[12]:R,u[13]:A,u[48]:B,u[25]:U,u[34]:V,u[49]:K-5*V,u[50]:5*U-L};bv=s.Matrix([u[i].subs(bs) for i in one]);qv=R*z*z+(2*A-B)*z*w-P*w*w;disc=(2*A-B)**2+4*P*R;den=25*A*A-25*A*B-14*B*B+25*P*R;assert s.expand(4*den-25*disc+81*B*B)==0
 z4,bu,Beta=primitive_beta(m);vec=s.Matrix([Beta[i,j].diff(z4[h]) for (i,j),h in COORDS]);basis=vec.jacobian(bu);rr=list(basis.T.rref()[1]);left=basis[rr,:].inv();actions=s.Matrix.hstack(*[left*infinitesimal_action(vec,H)[rr,:] for H in sp_basis()]);source_path=ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json';source=json.loads(source_path.read_text());extra_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json';extra=json.loads(extra_path.read_text())['universal_new_equation'];extra=[term for term in extra if sum(i in two for i in term[0])==2];sample={A:1,B:2,P:-1,R:0,U:2,V:3,K:0,L:0,**dict(zip(nuisance,[2,3,5,7,11,13]))};out={}
 for kind in ['M0','M1']:
  H=2*qv/(3*B) if kind=='M0' else -qv/(6*B)-z*w/2;G=a*x*x*z+b*x*x*w+x*H;Q=s.zeros(4);Q[0,1]=x*x if kind=='M0' else x*x/2;Q[2,3]=0 if kind=='M0' else x*x/2;Q[1,2]=-s.diff(G,w)-(x*z if kind=='M1' else 0);Q[1,3]=s.diff(G,z);Q-=Q.T;F=[s.factor(f) for f in beta_matrix(m,Q)*bv];M,ct=s.linear_eq_to_matrix(F,a,b);rows=[70,84];sol=(M[rows,:].inv()*ct[rows,:]).applyfunc(s.factor);assert (M*sol-ct).applyfunc(s.cancel)==s.zeros(134,1);pivot=s.factor(M[rows,:].det());Qsol=Q.subs({a:sol[0],b:sol[1]},simultaneous=True);Qsol=Qsol.applyfunc(s.cancel)
  PP=X[0]**2*J;PP[1:,1:]+=Q;raw=m.normalize(PP)['parameters'];values=[bv[one.index(i)] if i in one else raw[i] for i in range(61)];g6=s.factor(m.evaluate(extra,values));assert s.cancel(g6.subs({a:sol[0],b:sol[1]},simultaneous=True))==0
  k=s.Symbol('k');block=next(c for c in source['families'] if c['label']==kind);baseline={s.sympify(j):s.sympify(v) for j,v in block['baseline_solution'].items()};gs=s.symbols('g0:9');diagonal={P:0,R:0,A:k,B:2*k+1};assert s.cancel(sol[0].subs(diagonal)-baseline[gs[0]])==0;assert s.cancel(sol[1].subs(diagonal)-baseline[gs[1]])==0;assert s.cancel(s.expand(H).coeff(z,1).coeff(w,1).subs(diagonal)-baseline[gs[3]])==0
  PP=X[0]**2*J;PP[1:,1:]+=Qsol;raw=m.normalize(PP)['parameters'];q=s.Matrix([s.cancel(raw[i]) for i in two]);rankone={P:-1,R:0,A:1,B:2,K:0,L:0};Gsol=s.cancel(G.subs({a:sol[0],b:sol[1]},simultaneous=True));rkG=s.expand(Gsol.subs(rankone));rankone_saved={s.sympify(j):s.sympify(v).subs(gs[8],0) for j,v in block['rankone_solution'].items()};rankone_saved[gs[8]]=0;mons=[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3];assert s.expand(rkG-sum(rankone_saved[j]*h for j,h in zip(gs,mons)))==0
  point=s.zeros(61,1);beta=bv.subs(sample);q0=q.subs(sample)
  for i,v in zip(one,beta):point[i]=v
  for i,v in zip(two,q0):point[i]=v
  local=calculate(m,True,point);fixed=two.index(local['chart_fixed_u']-1);cols=local['chart_columns'];BD=bv.jacobian(params).subs(sample).row_join(actions.subs(dict(zip(bu,beta))));Dq=q.jacobian(params).subs(sample);Dq-=q0*Dq[fixed,:]/q0[fixed];Q0=Qsol.subs(sample);qd=[]
  for Hsp in sp_basis():
   dQ=sum(((Hsp*vv)[i]*Q0.diff(vv[i]) for i in range(4)),s.zeros(4))-Hsp*Q0-Q0*Hsp.T;PP=X[0]**2*J;PP[1:,1:]+=dQ;raw=m.normalize(PP)['parameters'];dq=s.Matrix([raw[i] for i in two]);dq-=q0*dq[fixed]/q0[fixed];qd.append(dq[cols,:])
  QD=Dq[cols,:].row_join(s.Matrix.hstack(*qd));N=s.Matrix.hstack(*BD.nullspace());C=QD*N;assert len(C.to_DM().rref()[1])==1;dd=s.Matrix([[s.diff(disc,v).subs(sample) for v in params]+[0]*10]);dn=dd*N;j=next(j for j in range(N.cols) if dn[j]!=0 and C[:,j]!=s.zeros(44,1));tv=C[:,j].col_join(s.zeros(1,1));T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);tr=list(T.T.to_DM().rref()[1]);tangent=T[tr,:].inv()*tv[tr,:];assert T*tangent==tv;bd=minor_witness(BD);restricted=[i for i in range(BD.cols) if i!=3];rbd=minor_witness(BD[:,restricted]);joint=minor_witness(BD.col_join(QD));assert bd['rank']==16 and rbd['rank']==16 and joint['rank']==17
  result=dict(H=str(H),qv=str(qv),solution_a=str(sol[0]),solution_b=str(sol[1]),linear_pivot_rows=rows,linear_pivot_matrix=enc(M[rows,:]),uncanceled_pivot_determinant=str(pivot),all134_identities=True,G6_before_solving=str(g6),G6_identity=True,q_coordinates=list(map(str,q)),rankone_G=str(rkG),local_data=local,beta_differential=enc(BD),q_differential=enc(QD),beta_minor=bd,rankone_beta_minor=rbd,rankone_source_column_indices=restricted,joint_minor=joint,source_direction=list(map(str,N[:,j])),discriminant_derivative=str(dn[j]),main_branch_tangent=list(map(str,tangent)))
  if kind=='M1':
   assert local['leading_jacobian_rank']==42 and local['homogenized_jacobian_rank']==43 and local['quadratic_span_rank']==1 and local['common_hessian_rank']==1;rad=s.Matrix(local['quadratic_common_radical']).applyfunc(s.Rational);assert s.Matrix.hstack(rad,tangent).rank()==1;result.update(completed_local_ring='C[[u,e]]/(e^2)',t_zero=True,whole_generic_beta_rankone_point_stratum_removed=True)
  out[kind]=result;print(kind,'all135 identities; rankone ranks',local['leading_jacobian_rank'],local['homogenized_jacobian_rank'],'dDisc',dn[j],flush=True)
 paths=[Path(__file__),source_path,extra_path,ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json']
 return dict(id='S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101',all_passed=True,field='QQ',source_variables=list(map(str,params)),source_beta=list(map(str,bv)),vertex_discriminant=str(disc),denominator_factor=str(den),sample={str(k):str(v) for k,v in sample.items()},cases=out,scope='Universal rational main formulas regular across rank-one vertex quadrics; M1 rank-one point stratum removed over general beta. M0 rank-one A1/junction not removed here.',proof_note='research/notes/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.md',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
if __name__=='__main__':main()

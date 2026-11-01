#!/usr/bin/env python3
"""Two double-root/mixed intersections, free rank5, square-zero nilradical, t=0."""
import argparse,json,time,hashlib
from pathlib import Path
import sympy as s
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_nil22_infinity_20261101 import beta_matrix
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
from s61_local_quartic_certificate_20261101 import quartic_profile
J=J0
ROOT=Path(__file__).resolve().parents[1]
def geometry(locals):
 m=S61();x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);Q=s.zeros(4);Q[1,2]=z*z;Q[1,3]=-2*z*w;Q=Q-Q.T
 one=[i for i,a in enumerate(m.reduced['parameter_weights']) if a==1];two=[i for i,a in enumerate(m.reduced['parameter_weights']) if a==2];M=beta_matrix(m,Q);H=s.Matrix.hstack(*M.nullspace());zz,bu,B=primitive_beta(m);vec=s.Matrix([B[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vec.jacobian(bu);rows=list(basis.T.rref()[1]);left=basis[rows,:].inv();acts=s.Matrix.hstack(*[(left*infinitesimal_action(vec,A)[rows,:]).applyfunc(s.expand) for A in sp_basis()])
 def coords(Q):
  P=X[0]**2*J0;P[1:,1:]+=Q;return s.Matrix([m.normalize(P)['parameters'][i] for i in two])
 q0=coords(Q);fixed=next(i for i,v in enumerate(q0) if v);chart=[i for i in range(45) if i!=fixed];orbit=[]
 for A in sp_basis():
  dQ=sum(((A*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T;dq=coords(dQ);dq-=q0*dq[fixed]/q0[fixed];orbit.append(dq[chart,:])
 QD=s.Matrix.hstack(*orbit);stab=s.Matrix.hstack(*QD.nullspace());enc=lambda M:[[str(x) for x in row] for row in M.tolist()];out={}
 for label,cons in [('B0',bu[one.index(48)]),('5A3B',5*bu[one.index(13)]-3*bu[one.index(48)])]:
  r=locals[label];point=s.Matrix(r['local']['point61']).applyfunc(s.Rational);sub=dict(zip(bu,[point[i] for i in one]));BD=H.row_join(acts.subs(sub));sourceQ=s.zeros(44,10).row_join(QD);C=sourceQ*s.Matrix.hstack(*BD.nullspace());cp=list(C.to_DM().rref()[1]);T=s.Matrix(r['local']['tangent_basis']).applyfunc(s.Rational);tr=list(T.T.to_DM().rref()[1]);tangent=T[tr,:].inv()*C[:,cp[0]].col_join(s.zeros(1,1))[tr,:]
  constraint=s.Matrix([cons]).jacobian(bu);HH=H*s.Matrix.hstack(*(constraint*H).nullspace());b=s.symbols('b:9');sa=(acts*stab).subs(dict(zip(bu,HH*s.Matrix(b))),simultaneous=True).applyfunc(s.expand);restricted=HH.row_join(acts.subs(sub));joint=restricted.col_join(s.zeros(44,9).row_join(QD));preserved=(constraint*sa).applyfunc(s.expand);cpp=list(restricted.to_DM().rref()[1]);dm=restricted[:,cpp].to_DM();minor=dm.domain.to_sympy(dm.det())
  assert BD.rank()==16 and len(cp)==1 and sum(tangent)!=0
  assert restricted.rank()==joint.rank()==16 and preserved==s.zeros(1,3) and minor!=0
  print(label,'curve',len(cp),'tan',tangent.T,'restricted',restricted.rank(),joint.rank(),'preserved',preserved,'minor',minor,flush=True)
  out[label]=dict(main_tangent=enc(tangent),main_beta_rank=BD.rank(),main_curve_rank=len(cp),restricted_beta_rank=restricted.rank(),restricted_joint_rank=joint.rank(),stabilizer_constraint=enc(preserved),minor=str(minor),restricted_beta_differential=enc(restricted),minor_columns=cpp)
 return out
def branches(locals):
 m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];x,y,z,w=X[1:];tau=s.Symbol('tau');out={}
 for label in ['B0','5A3B']:
  st=time.time();loc=locals[label]['local'];point=s.Matrix(loc['point61']).applyfunc(s.Rational);_,F=specialize(m,{i:point[i] for i in one},True);G=z*z*w+tau*(x*x*z if label=='B0' else x*w*w);Q=s.zeros(4);Q[1,2]=s.diff(G,w);Q[1,3]=-s.diff(G,z);Q-=Q.T;P=X[0]**2*J;P[1:,1:]+=Q;qs=m.normalize(P)['parameters'];q=s.Matrix([qs[i] for i in two]);q0=s.Matrix([point[i] for i in two]);assert q.subs(tau,0)==q0
  assert all(s.expand(sum(c*s.prod(q[i] for i in mon) for mon,c in f.items() if len(mon)==2))==0 for f in F)
  fixed=next(i for i,v in enumerate(q0) if v);cols=loc['chart_columns'];dq=q.diff(tau);dq-=q0*dq[fixed]/q0[fixed];T=s.Matrix(loc['tangent_basis']).applyfunc(s.Rational);tr=list(T.T.to_DM().rref()[1]);tangent=T[tr,:].inv()*dq[cols,:].col_join(s.zeros(1,1))[tr,:];assert T*tangent==dq[cols,:].col_join(s.zeros(1,1))
  p=point.copy()
  for i,v in zip(two,q.subs(tau,1)):p[i]=v
  sample=calculate(m,True,p);print(label,'sample',sample['leading_jacobian_rank'],sample['homogenized_jacobian_rank'],sample['quadratic_span_rank'],'tangent',tangent.T,flush=True)
  D=s.zeros(135,45)
  for row,f in enumerate(F):
   for mon,c in f.items():
    if len(mon)==2:
     i,j=mon;D[row,i]+=c*q[j];D[row,j]+=c*q[i]
  rk=len(D.to_DM().rref()[1]);print(label,'universal leading rank',rk,'seconds',time.time()-st,flush=True)
  out[label]=dict(G=str(G),q_curve=list(map(str,q)),tangent=list(map(str,tangent)),generic_leading_rank=rk,sample=sample,seconds=time.time()-st)
 return out
def universal():
 m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];u=s.symbols('u1:62');x,y,z,w=X[1:];out={}
 for label in ['B0','5A3B']:
  st=time.time();beta={i:u[i] for i in one}
  for i in [3,4,11,12]:beta[i]=0
  beta[49]=-5*u[34];beta[50]=5*u[25]
  if label=='B0':beta[48]=0
  else:beta[13]=3*u[48]/5
  _,F=specialize(m,beta,True);G=z*z*w+(x*x*z if label=='B0' else x*w*w);Q=s.zeros(4);Q[1,2]=s.diff(G,w);Q[1,3]=-s.diff(G,z);Q-=Q.T;P=X[0]**2*J;P[1:,1:]+=Q;qs=m.normalize(P)['parameters'];q=s.Matrix([qs[i] for i in two]);assert all(s.expand(sum(c*s.prod(q[i] for i in mon) for mon,c in f.items() if len(mon)==2))==0 for f in F)
  D=s.zeros(135,45)
  for row,f in enumerate(F):
   for mon,c in f.items():
    if len(mon)==2:
     i,j=mon;D[row,i]+=c*q[j];D[row,j]+=c*q[i]
  K=s.Matrix.hstack(*D[:70,:].nullspace());E=(D[70:,:]*K).applyfunc(s.expand);print(label,'first rank',45-K.cols,'E',E.shape,flush=True)
  dm=E.to_DM();red,piv=dm.rref();rank=45-K.cols+len(piv);print(label,'universal',rank,'seconds',time.time()-st,flush=True);out[label]={'first_rank':45-K.cols,'residual_rank':len(piv),'rank':rank,'kernel_first70':[[str(v) for v in row] for row in K.tolist()],'residual':[[str(v) for v in row] for row in E.tolist()],'seconds':time.time()-st};assert rank==42
 return out

def verify():
 start=time.time();m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];cases={}
 for label,sub in [('B0',{48:0}),('5A3B',{13:6})]:
  point=s.zeros(61,1)
  for j,i in enumerate(one):point[i]=j+2
  for i in [3,4,11,12]:point[i]=0
  point[13]=3;point[48]=10
  for i,v in sub.items():point[i]=v
  point[49]=-5*point[34];point[50]=5*point[25];point[45]=-2;r=quartic_profile(point)
  assert all(c['zero_dimensional'] and c['length']==5 for c in r['sections']['sum'].values());cases[label]=r
 local={k:v['cubic'] for k,v in cases.items()};geo=geometry(local);tails=branches(local);uni=universal()
 for label in cases:
  main=s.Matrix(geo[label]['main_tangent']).applyfunc(s.Rational);tail=s.Matrix(tails[label]['tangent']).applyfunc(s.Rational);assert s.Matrix.hstack(main,tail).rank()==2 and sum(main)!=0 and sum(tail)!=0
  sample=tails[label]['sample'];assert sample['leading_jacobian_rank']==42 and sample['homogenized_jacobian_rank']==43 and sample['quadratic_span_rank']==1 and sample['common_hessian_rank']==1 and sample['homogenizer_zero_on_tangent']
 paths=[Path(__file__),ROOT/'code/s61_local_quartic_certificate_20261101.py',ROOT/'code/s61_local_cubic_certificate_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101.json']
 return dict(id='S61_NONBINARY_DOUBLE_ROOT_MIXED_JUNCTIONS_20261101',all_passed=True,field='QQ',cases=cases,geometry=geo,branches=tails,universal_mixed_rank=uni,branch_generic_lengths=[3,2],completed_free_rank=5,homogenizer_zero=True,nilradical_square_zero=True,reduced_ordinary_node=True,whole_generic_beta_strata_removed=True,limitations=['Only the two double-root endpoint16D strata; not all mixed tails.','No explicit full multiplication table.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({k:r[k] for k in ['all_passed','completed_free_rank','seconds']}))

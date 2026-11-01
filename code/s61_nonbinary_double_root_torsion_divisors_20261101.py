#!/usr/bin/env python3
"""Four double-root divisors: generic length3, transverse length4, nilpotent t."""
import argparse,sys,json,time,hashlib
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X,J as J0
from s61_nil22_infinity_20261101 import beta_matrix
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
from s61_local_cubic_certificate_20261101 import profile
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
 for label,cons in [('AeqB',bu[one.index(13)]-bu[one.index(48)]),('5A7B',5*bu[one.index(13)]-7*bu[one.index(48)]),('U0',bu[one.index(25)]),('V0',bu[one.index(34)])]:
  r=locals[label];point=s.Matrix(r['local']['point61']).applyfunc(s.Rational);sub=dict(zip(bu,[point[i] for i in one]));BD=H.row_join(acts.subs(sub));sourceQ=s.zeros(44,10).row_join(QD);C=sourceQ*s.Matrix.hstack(*BD.nullspace());cp=list(C.to_DM().rref()[1]);T=s.Matrix(r['local']['tangent_basis']).applyfunc(s.Rational);tr=list(T.T.to_DM().rref()[1]);tangent=T[tr,:].inv()*C[:,cp[0]].col_join(s.zeros(1,1))[tr,:]
  constraint=s.Matrix([cons]).jacobian(bu);HH=H*s.Matrix.hstack(*(constraint*H).nullspace());b=s.symbols('b:9');sa=(acts*stab).subs(dict(zip(bu,HH*s.Matrix(b))),simultaneous=True).applyfunc(s.expand);restricted=HH.row_join(acts.subs(sub));joint=restricted.col_join(s.zeros(44,9).row_join(QD));preserved=(constraint*sa).applyfunc(s.expand);cpp=list(restricted.to_DM().rref()[1]);dm=restricted[:,cpp].to_DM();minor=dm.domain.to_sympy(dm.det())
  assert BD.rank()==16 and len(cp)==1 and tangent[0]!=0
  assert restricted.rank()==joint.rank()==16 and preserved==s.zeros(1,3) and minor!=0
  print(label,'curve',len(cp),'tan',tangent.T,'restricted',restricted.rank(),joint.rank(),'preserved',preserved,'minor',minor,flush=True)
  out[label]=dict(main_tangent=enc(tangent),main_beta_rank=BD.rank(),main_curve_rank=len(cp),restricted_beta_rank=restricted.rank(),restricted_joint_rank=joint.rank(),stabilizer_constraint=enc(preserved),minor=str(minor),restricted_beta_differential=enc(restricted),minor_columns=cpp)
 return out

def transverse(r):
 n=len(r['local']['tangent_basis'][0]);ms2=list(map(tuple,r['quadratic_monomials']));ms3=list(map(tuple,r['cubic_monomials']));ii2=[i for i,v in enumerate(ms2) if 0 not in v];ii3=[i for i,v in enumerate(ms3) if 0 not in v]
 mon=[(i,) for i in range(1,n)]+[ms2[i] for i in ii2]+[ms3[i] for i in ii3];ix={v:i for i,v in enumerate(mon)};A2=s.Matrix(r['A2']).applyfunc(s.Rational);A3=s.Matrix(r['A3']).applyfunc(s.Rational);rows=[]
 for k in range(135):
  rows.append([0]*(n-1)+list(A2[k,ii2])+list(A3[k,ii3]))
  for a in range(1,n):
   row=[s.S.Zero]*len(mon)
   for i in ii2:row[ix[tuple(sorted((a,)+ms2[i]))]]+=A2[k,i]
   if any(row):rows.append(row)
 R=s.Matrix(rows);rank=len(R.to_DM().rref()[1]);length=1+len(mon)-rank;assert length==4
 tq=list(map(s.Rational,r['homogenizer_quadratic']));tc=list(map(s.Rational,r['homogenizer_cubic']));t=[0]*(n-1)+[tq[i] for i in ii2]+[tc[i] for i in ii3]
 return dict(length=length,relation_rank=rank,t_nonzero=len(R.col_join(s.Matrix([t])).to_DM().rref()[1])>rank)
def verify():
 start=time.time();m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];locals={}
 for label,sub in [('AeqB',{13:10}),('5A7B',{13:14}),('U0',{25:0}),('V0',{34:0})]:
  point=s.zeros(61,1)
  for j,i in enumerate(one):point[i]=j+2
  for i in [3,4,11,12]:point[i]=0
  point[13]=3;point[48]=10
  for i,v in sub.items():point[i]=v
  point[49]=-5*point[34];point[50]=5*point[25];point[45]=-2
  r=profile(point)
  assert r['sections']['0']['zero_dimensional'] and r['sections']['0']['length']==4
  assert r['sections']['1']['zero_dimensional'] and r['sections']['1']['length']==3
  r['actual_transverse']=transverse(r);locals[label]=r
 geo=geometry(locals)
 paths=[Path(__file__),ROOT/'code/s61_local_cubic_certificate_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101.json']
 return dict(id='S61_NONBINARY_DOUBLE_ROOT_TORSION_DIVISORS_20261101',all_passed=True,field='QQ',cases=locals,geometry=geo,generic_rank=3,transverse_length=4,torsion_minimal_generators=1,homogenizer_nilpotent=True,whole_generic_beta_strata_removed=True,limitations=['No full t=0 claim or specified nilpotence exponent.','Only A=B,5A=7B,U=0,V=0 strata; B=0 and5A=3B remain separate.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({k:r[k] for k in ['all_passed','generic_rank','transverse_length','seconds']}))

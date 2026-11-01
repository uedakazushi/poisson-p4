#!/usr/bin/env python3
"""Uniform flag open and reduction of the residual locus to two invariant divisors."""
import argparse,json,time,hashlib
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import S61,X,J as J0
from s61_nil22_infinity_20261101 import beta_matrix
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
def rank(M):return len(M.to_DM().rref()[1])
def enc(M):return [[str(v) for v in row] for row in M.tolist()]
def flag_data():
 m=S61();x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);one=[i for i,v in enumerate(m.reduced['parameter_weights']) if v==1];two=[i for i,v in enumerate(m.reduced['parameter_weights']) if v==2];zz,bu,beta=primitive_beta(m);v=s.Matrix([beta[i,j].diff(zz[k]) for (i,j),k in COORDS]);bs=v.jacobian(bu);rows=list(bs.T.to_DM().rref()[1]);left=bs[rows,:].inv();acts=s.Matrix.hstack(*[(left*infinitesimal_action(v,A)[rows,:]).applyfunc(s.expand) for A in sp_basis()]);out={}
 def coords(Q):
  P=X[0]**2*J0;P[1:,1:]+=Q;return s.Matrix([m.normalize(P)['parameters'][i] for i in two])
 for label,G in [('FLAG',x*x*z)]:
  start=time.time();Q=s.zeros(4);Q[1,2]=s.diff(G,w);Q[1,3]=-s.diff(G,z);Q-=Q.T;M=beta_matrix(m,Q);H=s.Matrix.hstack(*M.nullspace());b=H*s.Matrix(range(2,H.cols+2));sub=dict(zip(bu,b));BD=H.row_join(acts.subs(sub));q0=coords(Q);fixed=next(i for i,v in enumerate(q0) if v);chart=[i for i in range(45) if i!=fixed];orbit=[]
  for A in sp_basis():
   dQ=sum(((A*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T;dq=coords(dQ);dq-=q0*dq[fixed]/q0[fixed];orbit.append(dq[chart,:])
  QD=s.Matrix.hstack(*orbit);point=s.zeros(61,1)
  for i,v in zip(one,b):point[i]=v
  for i,v in zip(two,q0):point[i]=v
  joint=BD.col_join(s.zeros(44,H.cols).row_join(QD));stats=dict(Hdim=H.cols,orbitdim=rank(QD),beta_rank=rank(BD),joint_rank=rank(joint))
  print(label,str(G),stats,time.time()-start,flush=True);out[label]=dict(G=str(G),stats=stats,beta_kernel=enc(H),q_orbit_differential=enc(QD),source_beta_differential=enc(BD),point61=list(map(str,point)))
  bb=s.symbols('b:'+str(H.cols));beta_all=H*s.Matrix(bb);_,F=specialize(m,dict(zip(one,beta_all)),True);assert all(s.expand(sum(c*s.prod(q0[i] for i in mon) for mon,c in f.items() if len(mon)==2))==0 for f in F)
  if label=='T3':
   for constraint in [bu[0],bu[1],bu[2],bu[3],5*bu[4]-bu[13],5*bu[9]+bu[14],5*bu[6]-bu[15]]:assert s.expand(constraint.subs(dict(zip(bu,beta_all)),simultaneous=True))==0
   assert H.cols==9
 return out

def linear(source):
 r=source;m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
 for label in ['FLAG']:
  start=time.time();case=r[label];H=s.Matrix(case['beta_kernel']).applyfunc(s.Rational);b=s.symbols('b:'+str(H.cols));q0=s.Matrix([s.Rational(case['point61'][i]) for i in two]);_,F=specialize(m,dict(zip(one,H*s.Matrix(b))),True);fixed=next(i for i,v in enumerate(q0) if v);cols=[i for i in range(45) if i!=fixed];J=s.zeros(135,45)
  for row,f in enumerate(F):
   for mon,c in f.items():
    if len(mon)==2:
     i,j=mon
     if i!=fixed:J[row,cols.index(i)]+=c*q0[j]
     if j!=fixed:J[row,cols.index(j)]+=c*q0[i]
    elif len(mon)==1:J[row,44]+=c*q0[mon[0]]
  rr,piv=J.to_DM().rref();rr=rr.to_Matrix();assert len(piv)==38;free=[i for i in range(45) if i not in piv];T=s.zeros(45,len(free))
  for j,f in enumerate(free):
   T[f,j]=1
   for k,p in enumerate(piv):T[p,j]=-rr[k,f]
  assert T[-1,:]==s.zeros(1,len(free));print(label,'rank38',time.time()-start,flush=True)
  lr,lp=J.T.to_DM().rref();lr=lr.to_Matrix();lf=[i for i in range(135) if i not in lp];L=s.zeros(len(lf),135)
  for j,f in enumerate(lf):
   L[j,f]=1
   for k,p in enumerate(lp):L[j,p]=-lr[k,f]
  dm=J[list(lp),list(piv)].to_DM();det=s.factor(dm.domain.to_sympy(dm.det()));print(label,'minor',det,time.time()-start,flush=True)
  data=dict(label=label,beta_kernel=enc(H),jacobian=enc(J),tangent=enc(T),left_kernel=enc(L),pivot_rows=list(lp),pivot_columns=list(piv),jacobian_minor=str(det),seconds=time.time()-start)

 return data
def quadratic(r,case):
 start=time.time();m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];b=s.symbols('b:15');H=s.Matrix(r['beta_kernel']).applyfunc(s.sympify);T=s.Matrix(r['tangent']).applyfunc(s.sympify);L=s.Matrix(r['left_kernel']).applyfunc(s.sympify);sample=dict(zip(b,range(2,17)));Te=T.subs(sample);Le=L.subs(sample);BD=s.Matrix(case['source_beta_differential']).applyfunc(s.Rational);QD=s.Matrix(case['q_orbit_differential']).applyfunc(s.Rational);C=s.zeros(44,15).row_join(QD)*s.Matrix.hstack(*BD.nullspace());cp=list(C.to_DM().rref()[1]);C=C[:,cp].col_join(s.zeros(1,len(cp)));tr=list(Te.T.to_DM().rref()[1]);CT=Te[tr,:].inv()*C[tr,:];assert Te*CT==C;params=list(CT.T.to_DM().rref()[1]);normal=[i for i in range(T.cols) if i not in params];assert len(params)==3 and len(normal)==4;T=T[:,normal]
 q=s.Matrix([s.Rational(case['point61'][i]) for i in two]);fixed=next(i for i,v in enumerate(q) if v);chart=[i for i in range(45) if i!=fixed];_,F=specialize(m,dict(zip(one,H*s.Matrix(b))),True);mons=list(combinations_with_replacement(range(4),2));Q=s.zeros(135,10);Te=T.subs(sample)
 for row,f in enumerate(F):
  for mon,c in f.items():
   if fixed in mon:continue
   i,j=[chart.index(k) for k in mon]+[44]*(2-len(mon));c=c.subs(sample) if hasattr(c,'subs') else c
   for k,(a,z) in enumerate(mons):Q[row,k]+=c*(Te[i,a]*Te[j,z]+(Te[i,z]*Te[j,a] if a!=z else 0))
 Ae=Le*Q;rows=list(Ae.T.to_DM().rref()[1]);assert len(rows)==10
 print('params',params,'normal',normal,'rows',rows,time.time()-start,flush=True)
 # Clear all rational denominators before multiplication, and use a polynomial domain.
 Ls=L[rows,:];td=s.lcm([s.denom(v) for v in T]);ld=s.lcm([s.denom(v) for v in Ls]);Tp=T.applyfunc(lambda v:s.cancel(v*td));Lp=Ls.applyfunc(lambda v:s.cancel(v*ld));Q=s.zeros(135,10)
 for row,f in enumerate(F):
  for mon,c in f.items():
   if fixed in mon:continue
   i,j=[chart.index(k) for k in mon]+[44]*(2-len(mon))
   for k,(a,z) in enumerate(mons):Q[row,k]+=c*(Tp[i,a]*Tp[j,z]+(Tp[i,z]*Tp[j,a] if a!=z else 0))
 Q=Q.applyfunc(s.expand);A=(Lp.to_DM()*Q.to_DM()).to_Matrix();A=A.applyfunc(lambda v:s.factor(s.cancel(v/(ld*td**2))));print('A',time.time()-start,flush=True)
 # Row-denominator cancellation keeps this determinant small.
 rowden=[s.lcm([s.denom(v) for v in A[i,:]]) for i in range(10)];AP=s.Matrix([[s.cancel(A[i,j]*rowden[i]) for j in range(10)] for i in range(10)]);dm=AP.to_DM();det=s.factor(s.cancel(dm.domain.to_sympy(dm.det())/s.prod(rowden)));print('det',det,time.time()-start,flush=True)
 enc=lambda M:[[str(v) for v in row] for row in M.tolist()];return dict(parameter_rows=params,normal_columns=normal,selected_left_rows=rows,quadratic_matrix=enc(A),quadratic_minor=str(det),tangent_clearing_factor=str(td),left_clearing_factor=str(ld),seconds=time.time()-start)
def orbit_frames(case):
 m=S61();z,bu,beta=primitive_beta(m);v=s.Matrix([beta[i,j].diff(z[k]) for (i,j),k in COORDS]);B=v.jacobian(bu);rr=list(B.T.to_DM().rref()[1]);left=B[rr,:].inv();acts=s.Matrix.hstack(*[(left*infinitesimal_action(v,A)[rr,:]).applyfunc(s.expand) for A in sp_basis()]);H=s.Matrix(case['beta_kernel']).applyfunc(s.Rational);QD=s.Matrix(case['q_orbit_differential']).applyfunc(s.Rational);stab=s.Matrix.hstack(*QD.nullspace());assert stab.cols==6
 b=s.symbols('b:15');hr=list(H.T.to_DM().rref()[1]);hl=H[hr,:].inv();act=acts.subs(dict(zip(bu,H*s.Matrix(b))),simultaneous=True);db=(hl*(act*stab)[hr,:]).applyfunc(s.expand)
 assert all(s.cancel(v/b[0]).is_number for v in db[0,:]);assert all(s.cancel(v/b[1]).is_number for v in db[1,:])
 sp=sp_basis();A=sum((stab[j,4]*sp[j] for j in range(10)),s.zeros(4));assert A*A==s.zeros(4)
 assert db[0,4]==db[1,4]==db[12,4]==0 and db[3,4]==b[0]
 ss=s.Symbol('s');g=s.eye(4)+ss*A;J4=J0[1:,1:];assert g.T*J4*g==J4
 x,y,z,w=X[1:];Q=s.zeros(4);Q[1,3]=-x*x;Q-=Q.T;assert A*Q+Q*A.T==s.zeros(4)
 boundary=[]
 for zero in [0,1]:
  HC=H[:,[j for j in range(15) if j!=zero]];vals=dict(zip(b,range(2,17)));vals[b[zero]]=0;BD=HC.row_join(act.subs(vals,simultaneous=True));joint=BD.col_join(s.zeros(44,14).row_join(QD));assert BD.rank()==16 and joint.rank()==18
  boundary.append(dict(zero_coordinate=zero,source_beta_differential=enc(BD),beta_rank=16,joint_rank=18))
 return dict(boundary=boundary,stabilizer=enc(stab),beta_derivative=enc(db),unipotent_generator=enc(A),translation='b3 -> b3+s*b0; b0,b1,b12 fixed',remaining_divisors=['b0=0','b1=0'])
def verify():
 start=time.time();source=flag_data();case=source['FLAG'];assert [case['stats'][k] for k in ['Hdim','orbitdim','beta_rank','joint_rank']]==[15,4,16,19]
 lin=linear(source);q=quadratic(lin,case);frames=orbit_frames(case);b=s.symbols('b:15');factors=lambda e:{str(s.Poly(f,*b).monic().as_expr()) for f,k in s.factor_list(e,*b)[1]}
 D=s.sympify(lin['jacobian_minor']);allowed=factors(D)
 for key in ['tangent','left_kernel']:
  for row in lin[key]:
   for value in set(row):assert factors(s.denom(s.cancel(s.sympify(value))))<=allowed
 det=s.sympify(q['quadratic_minor']);assert factors(s.numer(det))==factors(b[0]*b[1]*(5*b[3]-13*b[12]));assert factors(s.denom(det))<=allowed
 old=ROOT/'reports/S61_FLAG_FORMAL_RING_20261101.json';old_data=json.loads(old.read_text());assert old_data['all_passed'] and old_data['generic_ranks'][1]['generic_leading_jacobian_rank']==37
 paths=[Path(__file__),old,ROOT/'reports/S61_FLAG_FORMAL_RING_REVIEW_20261101.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
 return dict(id='S61_FLAG_UNIFORM_OPEN_20261101',all_passed=True,flag=case,linear=lin,quadratic=q,frame_change=frames,conclusion='Every flag point with b0*b1 !=0 has completed full135 ring k[[u0,u1,u2,e0,e1,e2,e3]]/(e0,e1,e2,e3)^2 and t=0. All remaining flags lie in the two invariant divisors b0=0 or b1=0.',proof_note='research/notes/S61_FLAG_UNIFORM_OPEN_20261101.md',limitations=['The two divisors have dominant beta images and remain unresolved.','No entire flag-orbit or full infinity removal.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({k:r[k] for k in ['all_passed','seconds']}))

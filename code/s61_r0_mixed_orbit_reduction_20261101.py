#!/usr/bin/env python3
"""Finite mixed-orbit cover of the remaining R0 normal-form faces."""
import argparse,json,time,hashlib
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
def orbit_data():
 m=S61();x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);one=[i for i,v in enumerate(m.reduced['parameter_weights']) if v==1];two=[i for i,v in enumerate(m.reduced['parameter_weights']) if v==2];zz,bu,beta=primitive_beta(m);v=s.Matrix([beta[i,j].diff(zz[k]) for (i,j),k in COORDS]);bs=v.jacobian(bu);rows=list(bs.T.to_DM().rref()[1]);left=bs[rows,:].inv();acts=s.Matrix.hstack(*[(left*infinitesimal_action(v,A)[rows,:]).applyfunc(s.expand) for A in sp_basis()]);out={}
 def coords(Q):
  P=X[0]**2*J0;P[1:,1:]+=Q;return s.Matrix([m.normalize(P)['parameters'][i] for i in two])
 for label,G in [('T1',z*z*w+x*x*z),('T2',z*z*w+x*w*w),('T3',z*(w*w+x*x)),('T4',x*z*w+w**3),('T5',x*x*z+x*w*w),('T6',x*z*z+w**3),('T7',x*x*z+w**3)]:
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

def verify():
 start=time.time();a,b,c,x,z,w=s.symbols('a b c x z w');source=ROOT/'reports/S61_R0_UNIFORM_VERTEX_20261101.json';r=json.loads(source.read_text());assert r['all_passed']
 swap={0:1,1:0,2:4,4:2,3:3,5:8,8:5,6:7,7:6}
 groups={(0,6):'T1',(4,6):'T2',(0,3,7):'T3-or-double-or-canonical-or-flag',(3,8):'T4',(0,1,4):'T5-or-binary/flag',(2,8):'T6',(0,1,8):'T7-or-binary/flag',(2,6):'double-or-binary'}
 covered=[]
 for k,faces in r['extra_vertical_prime_faces'].items():
  for face in faces:
   f=tuple(sorted(face));sf=tuple(sorted(swap[i] for i in f));key=f if f in groups else sf;assert key in groups,(k,face);covered.append(dict(k=k,face=face,canonical_face=list(key),reduction=groups[key]))
 assert len(covered)==16
 # All transformations below are Darboux shears or the SL2(z,w) shear.
 identities={
  'double_pair':s.expand((z*z*(a*x+b*w)).subs(w,w-a*x/b)-b*z*z*w),
  'three_line_plane':s.expand((z*(a*x*x+b*x*w+c*w*w)).subs(w,w-b*x/(2*c))-z*((a-b*b/(4*c))*x*x+c*w*w)),
  'canonical_boundary':s.expand((x*z*(a*x+b*w)).subs(w,w-a*x/b)-b*x*z*w),
  'T5_plane':s.expand((a*x*x*z+b*x*x*w+c*x*w*w).subs(z,z-b*w/a)-(a*x*x*z+c*x*w*w)),
  'T7_plane':s.expand((a*x*x*z+b*x*x*w+c*w**3).subs(z,z-b*w/a)-(a*x*x*z+c*w**3))}
 assert all(v==0 for v in identities.values())
 # Distinct torus weights permit nonzero coefficients of each two-term normal form to be normalized over C.
 representatives=[z*z*w+x*x*z,z*z*w+x*w*w,z*(w*w+x*x),x*z*w+w**3,x*x*z+x*w*w,x*z*z+w**3,x*x*z+w**3]
 weights=[]
 for G in representatives:
  mon=s.Poly(G,x,z,w).monoms();ww=[(m[0],m[1]-m[2]) for m in mon];assert len(ww)==2 and ww[0]!=ww[1];weights.append(ww)
 rho,lam=s.symbols('rho lam',nonzero=True);V=s.Matrix(s.symbols('x y z w'));J4=J0[1:,1:];D=s.diag(rho,1/rho,lam,1/lam);assert D.T*J4*D==J4
 shear1=s.Matrix([[1,0,0,0],[0,1,a,-b],[b,0,1,0],[a,0,0,1]]);shear2=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,a],[0,0,0,1]]);assert shear1.T*J4*shear1==J4 and shear2.T*J4*shear2==J4
 data=orbit_data();expected={'T1':[9,8,16,17],'T2':[9,8,16,17],'T3':[9,8,15,17],'T4':[9,8,16,17],'T5':[12,6,16,18],'T6':[9,8,16,17],'T7':[11,7,16,18]}
 for label,d in data.items():assert [d['stats'][k] for k in ['Hdim','orbitdim','beta_rank','joint_rank']]==expected[label]
 exclusion=ROOT/'reports/S61_NONBINARY_RESONANCE_EXCLUSIONS_20261101.json';assert json.loads(exclusion.read_text())['all_passed']
 paths=[Path(__file__),source,exclusion,ROOT/'reports/S61_R0_GENERIC_EXHAUSTION_20261101.json',ROOT/'reports/S61_NONBINARY_CANONICAL_EXHAUSTION_20261101.json',ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_EXHAUSTION_20261101.json',ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
 return dict(id='S61_R0_MIXED_ORBIT_REDUCTION_20261101',all_passed=True,field='QQ for identities; C for nonzero-coefficient torus normalization',covered_vertical_faces=covered,normalization_identities={k:str(v) for k,v in identities.items()},torus_weights=weights,orbits=data,excluded_general_beta_orbits=['T3'],remaining_mixed_orbits=['T1','T2','T4','T5','T6','T7'],conclusion='For sufficiently general beta, all R0 points remaining after the established canonical/double-root removals lie in the binary/flag loci or the six displayed mixed Sp4 orbits.',limitations=['This is a reduced set-theoretic orbit cover, not saturation of the six mixed orbits.','No pairwise distinct component count is inferred.','T5/T7 have generic fixed-beta surfaces; a uniform curve argument does not cover them.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({k:r[k] for k in ['all_passed','remaining_mixed_orbits','seconds']}))

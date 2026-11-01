#!/usr/bin/env python3
"""Universal highest-Jacobian bounds on four mixed R0 curve orbits."""
import argparse,json,time,hashlib
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]
def enc(M):return [[str(v) for v in row] for row in M.tolist()]
def verify():
 start=time.time();path=ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json';r=json.loads(path.read_text());assert r['all_passed'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];cases={}
 for label in ['T6','T4','T1','T2']:
  st=time.time();c=r['orbits'][label];H=s.Matrix(c['beta_kernel']).applyfunc(s.Rational);b=s.symbols('b:'+str(H.cols));point=s.Matrix(c['point61']).applyfunc(s.Rational);q=s.Matrix([point[i] for i in two]);_,F=specialize(m,dict(zip(one,H*s.Matrix(b))),True);D=s.zeros(135,45)
  for row,f in enumerate(F):
   assert s.expand(sum(coef*s.prod(q[i] for i in mon) for mon,coef in f.items() if len(mon)==2))==0
   for mon,coef in f.items():
    if len(mon)==2:
     i,j=mon;D[row,i]+=coef*q[j];D[row,j]+=coef*q[i]
  K=s.Matrix.hstack(*D[:70,:].nullspace());E=(D[70:,:]*K).applyfunc(s.expand);rr,piv=E.to_DM().rref();rr=rr.to_Matrix();free=[i for i in range(E.cols) if i not in piv];N=s.zeros(E.cols,len(free))
  for j,col in enumerate(free):
   N[col,j]=1
   for k,p in enumerate(piv):N[p,j]=-rr[k,col]
  assert len(free)==3 and (E.to_DM()*N.to_DM()).is_zero_matrix;KN=(K.to_DM()*N.to_DM()).to_Matrix();assert (D.to_DM()*KN.to_DM()).is_zero_matrix
  # Three affine homogeneous kernel vectors include the radial vector; hence projective tangent dimension>=2.
  first=45-K.cols;assert first+len(piv)==42
  cases[label]=dict(G=c['G'],beta_kernel=enc(H),first70_kernel=enc(K),residual_matrix=enc(E),residual_kernel=enc(N),full_highest_kernel=enc(KN),first70_rank=first,residual_rank=len(piv),universal_highest_rank_upper=42,generic_beta_fiber_dimension=1,seconds=time.time()-st)
  print(label,first,len(piv),time.time()-st,flush=True)
 paths=[Path(__file__),path,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
 return dict(id='S61_R0_MIXED_CURVE_UNIFORM_RANK_20261101',all_passed=True,field='QQ(b0,...,b8)',cases=cases,conclusion='On all four full beta kernels, every43-minor of the highest135x45 q-Jacobian vanishes identically. Along any actual one-dimensional fixed-beta branch of these orbits, S/(t) has generic local length at least2 whenever that factor is Artinian.',limitations=['No full local length upper bound or saturation removal.','The two surface orbits T5/T7 are separate.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({k:r[k] for k in ['all_passed','seconds']}))

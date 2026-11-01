#!/usr/bin/env python3
"""Exhaust four complete mixed R0 curve orbits over sufficiently general beta."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
def prime_factors(expr,b):
 return {str(s.Poly(f,*b).monic().as_expr()) for f,k in s.factor_list(expr,*b)[1]}
def compute(source):
 r=source['orbits'];cases={};m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
 for label in ['T6','T1','T2','T4']:
  start=time.time();case=r[label];H=s.Matrix(case['beta_kernel']).applyfunc(s.Rational);b=s.symbols('b:'+str(H.cols));q0=s.Matrix([s.Rational(case['point61'][i]) for i in two]);_,F=specialize(m,dict(zip(one,H*s.Matrix(b))),True);fixed=next(i for i,v in enumerate(q0) if v);cols=[i for i in range(45) if i!=fixed];J=s.zeros(135,45)
  for row,f in enumerate(F):
   for mon,c in f.items():
    if len(mon)==2:
     i,j=mon
     if i!=fixed:J[row,cols.index(i)]+=c*q0[j]
     if j!=fixed:J[row,cols.index(j)]+=c*q0[i]
    elif len(mon)==1:J[row,44]+=c*q0[mon[0]]
  rr,piv=J.to_DM().rref();rr=rr.to_Matrix();assert len(piv)==43;free=[i for i in range(45) if i not in piv];T=s.zeros(45,2)
  for j,f in enumerate(free):
   T[f,j]=1
   for k,p in enumerate(piv):T[p,j]=-rr[k,f]
  assert T[-1,:]==s.zeros(1,2);print(label,'rank43',time.time()-start,flush=True)
  lr,lp=J.T.to_DM().rref();lr=lr.to_Matrix();lf=[i for i in range(135) if i not in lp];L=s.zeros(len(lf),135)
  for j,f in enumerate(lf):
   L[j,f]=1
   for k,p in enumerate(lp):L[j,p]=-lr[k,f]
  Q=s.zeros(135,3);mons=[(0,0),(0,1),(1,1)]
  for row,f in enumerate(F):
   for mon,c in f.items():
    if fixed in mon:continue
    i,j=[cols.index(k) for k in mon]+[44]*(2-len(mon))
    for k,(a,z) in enumerate(mons):Q[row,k]+=c*(T[i,a]*T[j,z]+(T[i,z]*T[j,a] if a!=z else 0))
  A=(L.to_DM()*Q.to_DM()).to_Matrix().applyfunc(s.factor);assert len(A.to_DM().rref()[1])==1
  qrows=[i for i in range(A.rows) if any(A[i,:])];chosen=qrows[0];assert s.cancel(A[chosen,1]**2-4*A[chosen,0]*A[chosen,2])==0
  print(label,'quadratic',qrows,[str(x) for x in A[chosen,:]],time.time()-start,flush=True)
  dm=J[list(lp),list(piv)].to_DM();det=s.factor(dm.domain.to_sympy(dm.det()));print(label,'minor',det,time.time()-start,flush=True)
  data=dict(label=label,beta_kernel=enc(H),jacobian=enc(J),tangent=enc(T),left_kernel=enc(L),quadratic=enc(A),quadratic_nonzero_rows=qrows,pivot_rows=list(lp),pivot_columns=list(piv),jacobian_minor=str(det),seconds=time.time()-start);cases[label]=data
  # Every rational certificate denominator is a unit on the stated minor open.
  allowed=prime_factors(det,b)
  for mat in (T,L,A):
   for value in set(mat):
    assert prime_factors(s.cancel(value).as_numer_denom()[1],b)<=allowed
  assert A[:,0:2]==s.zeros(A.rows,2)
  extra=prime_factors(A[chosen,2].as_numer_denom()[0],b)
  data['bad_factors']=sorted(allowed|extra)

 m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];z,bu,beta=primitive_beta(m);v=s.Matrix([beta[i,j].diff(z[k]) for (i,j),k in COORDS]);base=v.jacobian(bu);rr=list(base.T.to_DM().rref()[1]);left=base[rr,:].inv();acts=s.Matrix.hstack(*[(left*infinitesimal_action(v,A)[rr,:]).applyfunc(s.expand) for A in sp_basis()]);source_cases=source['orbits'];out={}
 for label in ['T6','T2','T1','T4']:
  b=s.symbols('b:9');r=source_cases[label];H=s.Matrix(r['beta_kernel']).applyfunc(s.Rational);hr=list(H.T.to_DM().rref()[1]);hl=H[hr,:].inv();C=s.Matrix.hstack(*H.T.nullspace()).T;act=acts.subs(dict(zip(bu,H*s.Matrix(b))),simultaneous=True).applyfunc(s.expand);N0=C*act;db=hl*act[hr,:]
  if label=='T6':factors=[b[6],325*b[4]*b[6]+143*b[5]*b[6]+70*b[7]*b[8]]
  if label=='T2':factors=[b[6],105*b[4]*b[6]+63*b[5]*b[6]+23*b[7]*b[8]]
  if label=='T1':factors=[b[0],b[7],b[8],15*b[0]*b[6]-16*b[7]*b[8]]
  if label=='T4':factors=[b[6],b[7],105*b[4]*b[6]+42*b[5]*b[6]+23*b[7]*b[8]]
  output=[]
  assert {str(s.Poly(f,*b).monic().as_expr()) for f in factors}==set(cases[label]['bad_factors'])
  for f in factors:
   start=time.time();coord=f if f in b else (b[6] if label=='T1' else b[4]);sol=s.solve(f,coord)[0];sub={coord:sol};N=N0.col_join(s.Matrix([f]).jacobian(b)*db).subs(sub,simultaneous=True).applyfunc(s.cancel);rank=len(N.to_DM().rref()[1]);item=dict(factor=str(f),solve_variable=str(coord),solve_value=str(sol),normal_rank=rank,beta_image_dimension=8+rank,normal_action=[[str(x) for x in row] for row in N.tolist()]);print(label,str(f),'image',8+rank,'seconds',time.time()-start,flush=True)
   assert rank<8
   if s.Poly(f,*b).total_degree()==2:
    assert s.hessian(f,b).rank()==4
   item['geometrically_irreducible']=True
   output.append(item)
  out[label]=output
 return cases,out
def verify():
 start=time.time();path=ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json';source=json.loads(path.read_text());assert source['all_passed']
 cases,exceptions=compute(source)
 rank_report=ROOT/'reports/S61_R0_MIXED_CURVE_UNIFORM_RANK_20261101.json'
 assert json.loads(rank_report.read_text())['all_passed']
 paths=[Path(__file__),path,rank_report,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
 return dict(id='S61_R0_MIXED_CURVE_EXHAUSTION_20261101',all_passed=True,cases=cases,exceptional_factors=exceptions,conclusion='For sufficiently general beta every point of each of the four entire mixed curve orbits T1,T2,T4,T6 has full135 completed local ring k[[u,e]]/(e^2), with t=0. All are removed by t-saturation.',proof_note='research/notes/S61_R0_MIXED_CURVE_EXHAUSTION_20261101.md',limitations=['The two mixed surfaces T5,T7 and the binary/flag endpoints are separate.','This does not prove affine finiteness or five-point exhaustion.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({k:r[k] for k in ['all_passed','seconds']}))

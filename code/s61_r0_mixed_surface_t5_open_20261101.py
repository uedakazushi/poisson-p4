#!/usr/bin/env python3
"""Exact T5 transverse length23, dimension18 component, and generic ambient removal."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate
ROOT=Path(__file__).resolve().parents[1]
def enc(M):return [[str(v) for v in row] for row in M.tolist()]
def point_data():
 r=json.loads((ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json').read_text())['orbits']['T5'];m=S61();p=s.Matrix(r['point61']).applyfunc(s.Rational);loc=calculate(m,True,p)
 assert loc['homogenized_jacobian_rank']==40 and loc['quadratic_span_rank']==1
 J=s.Matrix(loc['jacobian']).applyfunc(s.Rational);T=s.Matrix(loc['tangent_basis']).applyfunc(s.Rational);H=s.Matrix(r['beta_kernel']).applyfunc(s.Rational);BD=s.Matrix(r['source_beta_differential']).applyfunc(s.Rational);QD=s.Matrix(r['q_orbit_differential']).applyfunc(s.Rational)
 C=s.zeros(44,H.cols).row_join(QD)*s.Matrix.hstack(*BD.nullspace());cp=list(C.to_DM().rref()[1]);assert len(cp)==2;C=C[:,cp].col_join(s.zeros(1,2));rr=list(T.T.to_DM().rref()[1]);CT=T[rr,:].inv()*C[rr,:];assert T*CT==C;params=list(CT.T.to_DM().rref()[1]);assert params==[0,1]
 joint=BD.col_join(s.zeros(44,H.cols).row_join(QD));assert BD.rank()==16 and joint.rank()==18
 return dict(local=dict(local=loc),parameter_rows=params,actual_surface_tangent=enc(CT),beta_rank=16,joint_rank=18)
def compute():
 r={'T5':point_data()};m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];out={}
 for label in ['T5']:
  start=time.time();case=r[label];loc=case['local']['local'];point=s.Matrix(loc['point61']).applyfunc(s.Rational);_,F=specialize(m,{i:point[i] for i in one},True);cols=loc['chart_columns'];fixed=next(i for i in range(45) if i not in cols);J=s.Matrix(loc['jacobian']).applyfunc(s.Rational);T=s.Matrix(loc['tangent_basis']).applyfunc(s.Rational);normal=[i for i in range(T.cols) if i not in case['parameter_rows']];T=T[:,normal];n=T.cols;R=list(J.T.to_DM().rref()[1]);C=list(J.to_DM().rref()[1]);U=s.zeros(45,len(R));inv=J[R,C].inv(method='DM')
  for k,c in enumerate(C):U[c,:]=inv[k,:]
  JU=J*U;ms={d:list(combinations_with_replacement(range(n),d)) for d in range(1,9)};index={d:{mon:i for i,mon in enumerate(ms[d])} for d in ms};W={1:T};A={};terms=[]
  for f in F:
   terms.append([([cols.index(i) for i in mon]+[44]*(2-len(mon)),coef) for mon,coef in f.items() if fixed not in mon])
  v=s.symbols('v:'+str(n));gens=[];stats=[]
  for degree in range(2,7):
   Q=s.zeros(135,len(ms[degree]))
   for row,ts in enumerate(terms):
    for (i,j),coef in ts:
     for a in range(1,degree):
      b=degree-a
      for k,ma in enumerate(ms[a]):
       if not W[a][i,k]:continue
       for l,mb in enumerate(ms[b]):
        if W[b][j,l]:Q[row,index[degree][tuple(sorted(ma+mb))]]+=coef*W[a][i,k]*W[b][j,l]
   A[degree]=Q-JU*Q[R,:];W[degree]=-U*Q[R,:]
   # Form the truncated actual ideal and its initial degree-d consequence space.
   mons=sum([ms[d] for d in range(2,degree+1)],[]);ix={mon:i for i,mon in enumerate(mons)};base=s.Matrix.hstack(*[A[d] for d in range(2,degree+1)]);rows=[]
   for row in base.tolist():
    if not any(row):continue
    for md in range(degree-1):
     mults=[()] if md==0 else ms[md]
     for mul in mults:
      rr=[s.S.Zero]*len(mons)
      for c,mon in zip(row,mons):
       if c and len(mon+mul)<=degree:rr[ix[tuple(sorted(mon+mul))]]+=c
      if any(rr):rows.append(rr)
   if rows:
    MM=s.Matrix(rows);lo=MM[:,:len(mons)-len(ms[degree])];hi=MM[:,len(mons)-len(ms[degree]):];pr=list(lo.T.to_DM().rref()[1]);pc=list(lo.to_DM().rref()[1]);rel=hi-lo[:,pc]*lo[pr,pc].inv(method='DM')*hi[pr,:];dr=list(rel.T.to_DM().rref()[1]);new=[s.expand(sum(c*s.prod(v[i] for i in mon) for c,mon in zip(row,ms[degree]))) for row in rel[dr,:].tolist()];gens+=new
   else:dr=[]
   gb=s.groebner(gens,*v,order='grevlex') if gens else None;counts=[]
   if gb is not None and gb.is_zero_dimensional:
    lead=[f.LM(order=gb.order).exponents for f in gb.polys]
    for d in range(20):
     exps=[tuple(mon.count(i) for i in range(n)) for mon in combinations_with_replacement(range(n),d)];count=sum(not any(all(a<=b for a,b in zip(l,e)) for l in lead) for e in exps);counts.append(count)
     if not count:break
   info=dict(degree=degree,initial_rank=len(dr),transverse_counts=counts,length=sum(counts),homogenizer=list(map(str,W[degree][-1,:])),seconds=time.time()-start);stats.append(info);print(label,'degree',degree,'initial_rank',len(dr),'length',sum(counts),'seconds',time.time()-start,flush=True)

  MM=s.Matrix(rows);rk=len(MM.to_DM().rref()[1]);length=1+n+len(mons)-rk;trow=s.Matrix.hstack(*[W[d][-1,:] for d in range(2,7)]);rkt=len(MM.col_join(trow).to_DM().rref()[1]);print('ACTUAL',length,'rank',rk,'t_nonzero',rkt>rk,flush=True)
  assert length==23 and rkt==rk and stats[-1]['transverse_counts']==[1,3,5,6,5,3,0]
  out[label]=dict(point_data=case,stats=stats,initial_generators=list(map(str,gens)),actual_transverse_length=length,actual_relation_rank=rk,homogenizer_nonzero=rkt>rk,residual_coefficients={str(d):[[str(v) for v in row] for row in A[d].tolist()] for d in A})
 return out
def verify():
 start=time.time();r=compute()['T5'];paths=[Path(__file__),ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py']
 return dict(id='S61_R0_MIXED_SURFACE_T5_OPEN_20261101',all_passed=True,field='QQ',point=r,fixed_beta_local_dimension=2,total_local_dimension=18,incidence_closure_is_irreducible_component=True,ambient_nonempty_open_removed=True,whole_generic_beta_surface_removed=False,proof_note='research/notes/S61_R0_MIXED_SURFACE_T5_OPEN_20261101.md',limitations=['Only the transverse quotient has proved length23 and t=0.','No free rank23 description or t=0 at the full specified local ring.','Dimension comparison removes a nonempty ambient open of the18D incidence; its17D complement can still dominate beta.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({k:r[k] for k in ['all_passed','fixed_beta_local_dimension','total_local_dimension','seconds']}))

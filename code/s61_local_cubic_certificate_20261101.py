#!/usr/bin/env python3
"""Exact cubic initial-ideal certificate on a specified local chart."""
import sys,time,json
import sympy as s
from itertools import combinations_with_replacement
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate

def profile(point):
 start=time.time();m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];local=calculate(m,True,point);_,F=specialize(m,{i:point[i] for i in one},True);cols=local['chart_columns'];fixed=next(i for i in range(45) if i not in cols)
 J=s.Matrix(local['jacobian']).applyfunc(s.Rational);T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);n=T.cols;rows=list(J.T.to_DM().rref()[1]);piv=list(J.to_DM().rref()[1]);U=s.zeros(45,len(rows));iv=J[rows,piv].inv(method='DM')
 for k,c in enumerate(piv):U[c,:]=iv[k,:]
 JU=J*U;ms={d:list(combinations_with_replacement(range(n),d)) for d in [2,3]};ix={mon:i for i,mon in enumerate(ms[3])};Q=s.zeros(135,len(ms[2]));terms=[]
 for r,f in enumerate(F):
  ts=[]
  for mon,c in f.items():
   if fixed in mon:continue
   i,j=[cols.index(k) for k in mon]+[44]*(2-len(mon));ts.append((i,j,c))
   for k,(a,b) in enumerate(ms[2]):Q[r,k]+=c*(T[i,a]*T[j,b]+(T[i,b]*T[j,a] if a!=b else 0))
  terms.append(ts)
 w2=-U*Q[rows,:];A2=Q-JU*Q[rows,:];C3=s.zeros(135,len(ms[3]))
 for r,ts in enumerate(terms):
  for i,j,c in ts:
   for a in range(n):
    for k,mon in enumerate(ms[2]):C3[r,ix[tuple(sorted((a,)+mon))]]+=c*(T[i,a]*w2[j,k]+T[j,a]*w2[i,k])
 A3=C3-JU*C3[rows,:];qr=list(A2.T.to_DM().rref()[1]);qc=list(A2.to_DM().rref()[1]);C=A2[:,qc]*A2[qr,qc].inv(method='DM');assert C*A2[qr,:]==A2;B=A3-C*A3[qr,:];cr=list(B.T.to_DM().rref()[1]);v=s.symbols('v:'+str(n));f2=[sum(c*s.prod(v[i] for i in mon) for c,mon in zip(row,ms[2])) for row in A2[qr,:].tolist()];f3=[sum(c*s.prod(v[i] for i in mon) for c,mon in zip(row,ms[3])) for row in B[cr,:].tolist()];print('initial',n,len(qr),len(cr),time.time()-start,flush=True)
 result=dict(local=local,A2=[[str(x) for x in row] for row in A2.tolist()],A3=[[str(x) for x in row] for row in A3.tolist()],initial_quadrics=list(map(str,f2)),initial_cubics=list(map(str,f3)),parameter='v0',sections={})
 for value in [0,1]:
  sub={v[0]:value};gb=s.groebner([s.expand(f.subs(sub)) for f in f2+f3],*v[1:],order='grevlex');counts=[]
  if gb.is_zero_dimensional:
   leads=[f.LM(order=gb.order).exponents for f in gb.polys]
   for d in range(30):
    exps=[tuple(mon.count(i) for i in range(n-1)) for mon in combinations_with_replacement(range(n-1),d)];count=sum(not any(all(a<=b for a,b in zip(l,e)) for l in leads) for e in exps);counts.append(count)
    if count==0:break
   assert counts[-1]==0
  result['sections'][str(value)]=dict(zero_dimensional=gb.is_zero_dimensional,counts=counts,length=sum(counts),groebner=list(map(str,gb)))
  print('section',value,result['sections'][str(value)]['length'],counts,time.time()-start,flush=True)
 result['homogenizer_quadratic']=list(map(str,w2[-1,:]));result['homogenizer_cubic']=list(map(str,(-U*C3[rows,:])[-1,:]));result['quadratic_monomials']=ms[2];result['cubic_monomials']=ms[3];result['seconds']=time.time()-start;return result

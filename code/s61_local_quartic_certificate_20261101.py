#!/usr/bin/env python3
"""Exact quartic initial relations, including lower-degree syzygies."""
import sys,json,time
import sympy as s
from itertools import combinations_with_replacement
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_local_cubic_certificate_20261101 import profile

def quartic_profile(point):
 start=time.time();old=profile(point);loc=old['local'];label='quartic';m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];_,F=specialize(m,{i:point[i] for i in one},True);cols=loc['chart_columns'];fixed=next(i for i in range(45) if i not in cols)
 J=s.Matrix(loc['jacobian']).applyfunc(s.Rational);T=s.Matrix(loc['tangent_basis']).applyfunc(s.Rational);n=T.cols;rows=list(J.T.to_DM().rref()[1]);piv=list(J.to_DM().rref()[1]);U=s.zeros(45,len(rows));iv=J[rows,piv].inv(method='DM')
 for k,c in enumerate(piv):U[c,:]=iv[k,:]
 JU=J*U;ms={d:list(combinations_with_replacement(range(n),d)) for d in [2,3,4]};ix={d:{mon:i for i,mon in enumerate(ms[d])} for d in [3,4]};Q=s.zeros(135,len(ms[2]));terms=[]
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
    for k,mon in enumerate(ms[2]):C3[r,ix[3][tuple(sorted((a,)+mon))]]+=c*(T[i,a]*w2[j,k]+T[j,a]*w2[i,k])
 w3=-U*C3[rows,:];A3=C3-JU*C3[rows,:];C4=s.zeros(135,len(ms[4]))
 for r,ts in enumerate(terms):
  for i,j,c in ts:
   for a in range(n):
    for k,mon in enumerate(ms[3]):C4[r,ix[4][tuple(sorted((a,)+mon))]]+=c*(T[i,a]*w3[j,k]+T[j,a]*w3[i,k])
   for k,m1 in enumerate(ms[2]):
    for l,m2 in enumerate(ms[2]):C4[r,ix[4][tuple(sorted(m1+m2))]]+=c*w2[i,k]*w2[j,l]
 A4=C4-JU*C4[rows,:];print(label,'quartic built',time.time()-start,flush=True)
 mon=ms[2]+ms[3]+ms[4];inds={a:i for i,a in enumerate(mon)};Rs=[]
 for row in A2.row_join(A3).row_join(A4).tolist():
  if not any(row):continue
  for mult in [()]+[(i,) for i in range(n)]+ms[2]:
   r=[s.S.Zero]*len(mon)
   for c,mm in zip(row,mon):
    if c and len(mult+mm)<=4:r[inds[tuple(sorted(mult+mm))]]+=c
   if any(r):Rs.append(r)
 R=s.Matrix(Rs);low=R[:,:len(ms[2])+len(ms[3])];high=R[:,len(ms[2])+len(ms[3]):];rr=list(low.T.to_DM().rref()[1]);cc=list(low.to_DM().rref()[1]);C=low[:,cc]*low[rr,cc].inv(method='DM');assert C*low[rr,:]==low;B4=high-C*high[rr,:];pr=list(B4.T.to_DM().rref()[1]);v=s.symbols('v:'+str(n));f4=[sum(c*s.prod(v[i] for i in mon) for c,mon in zip(row,ms[4])) for row in B4[pr,:].tolist()];print(label,'initial4',len(pr),time.time()-start,flush=True)
 forms=list(map(s.sympify,old['initial_quadrics']+old['initial_cubics']))+f4;out={'initial_quartics':list(map(str,f4)),'sections':{}}
 for param in ['sum']:
  results={}
  for value in [0,1]:
   idx=param if isinstance(param,int) else 0;others=[v[i] for i in range(n) if i!=idx];sub={v[idx]:value if isinstance(param,int) else value-sum(others)};gb=s.groebner([s.expand(f.subs(sub)) for f in forms],*others,order='grevlex');counts=[]
   if gb.is_zero_dimensional:
    leads=[f.LM(order=gb.order).exponents for f in gb.polys]
    for d in range(30):
     exps=[tuple(mon.count(i) for i in range(n-1)) for mon in combinations_with_replacement(range(n-1),d)];count=sum(not any(all(a<=b for a,b in zip(l,e)) for l in leads) for e in exps);counts.append(count)
     if not count:break
   results[str(value)]={'zero_dimensional':gb.is_zero_dimensional,'counts':counts,'length':sum(counts),'groebner':list(map(str,gb))}
  out['sections'][str(param)]=results;print(label,param,[(k,r['length'],r['zero_dimensional']) for k,r in results.items()],flush=True)
 out.update(cubic=old,A4=[[str(x) for x in row] for row in A4.tolist()],quartic_monomials=ms[4],homogenizer_quartic=list(map(str,(-U*C4[rows,:])[-1,:])),seconds=time.time()-start)
 return out

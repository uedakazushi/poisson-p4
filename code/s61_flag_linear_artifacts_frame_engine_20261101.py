#!/usr/bin/env python3
"""Polynomial frames at the three linear pivot-artifact cores."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
import json,pickle,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
st=time.time();WD=Path(sys.argv[1]);base=json.load(open(WD/'flag-b1-normalized-frame.json'));spec=json.load(open(ROOT/'reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json'))['core_linear_ranks'][:3];v=list(map(s.Symbol,base['variables']));r=v[0];J=s.Matrix(base['jacobian']).applyfunc(s.sympify);enc=lambda M:[[str(x)for x in row]for row in M.tolist()]
for k,item in enumerate(spec):
 value=s.Rational(item['core']);jj=J.subs(r,value);rows=item['pivot_rows'];cols=item['pivot_columns'];n=len(cols);B=jj[rows,cols];B0=B.subs({x:0 for x in v[1:]});inv0=B0.inv();N=(inv0*(B-B0)).applyfunc(s.expand);edges={i:{j for j in range(n)if N[i,j]}for i in range(n)};rem=set(range(n));order=[]
 while rem:
  leaves=[i for i in rem if not(edges[i]&rem)]
  if not leaves:break
  order+=leaves;rem-=set(leaves)
 print('core',value,'rank',n,'acyclic',not rem,'edges',sum(map(len,edges.values())),time.time()-st,flush=True)
 assert not rem
 inverse=s.zeros(n)
 for i in order:
  for j in range(n):inverse[i,j]=s.expand(inv0[i,j]-sum(N[i,l]*inverse[l,j]for l in edges[i]))
 assert (B*inverse-s.eye(n)).applyfunc(s.expand)==s.zeros(n)
 free=[i for i in range(45)if i not in cols];T=s.zeros(45,len(free));block=-(inverse*jj[rows,free]).applyfunc(s.expand)
 for i,col in enumerate(free):T[col,i]=1
 for i,col in enumerate(cols):T[col,:]=block[i,:]
 assert (jj*T).applyfunc(s.expand)==s.zeros(135,len(free));wd=WD/('core-'+str(k));wd.mkdir(exist_ok=True)
 subs={'b0':'1','b1':'0','b3':'0','b5':'0','b6':'0','b12':'1','b14':'0','b2':str(value)}
 (wd/'flag-b1-normalized-frame.json').write_text(json.dumps(dict(substitutions=subs,variables=list(map(str,v[1:])),jacobian=enc(jj),tangent=enc(T),pivot_rows=rows,pivot_columns=cols,free=free,minor=item['minor'],tangent_denominator='1'))+'\n')
 (wd/'flag-b1-triangular-inverse.json').write_text(json.dumps(dict(inverse=enc(inverse),core_denominator='1',dependency_order=order))+'\n')
 print('CORE FRAME PASS',value,'Tdegree',max(s.Poly(x,*v[1:]).total_degree()for x in T),'seconds',time.time()-st,flush=True)

"""Exact cubic initial-ideal bound for a specified local tangent chart."""
from itertools import combinations_with_replacement
import sympy as s
from s61_flag_transverse_obstruction_20261101 import specialize
def cubic_bound(m,local,branches,expected):
 p=list(map(s.Rational,local['point61']));beta={i:p[i] for i,w in enumerate(m.reduced['parameter_weights']) if w==1};_,F=specialize(m,beta,True);cols=local['chart_columns'];fixed=next(i for i in range(45) if i not in cols);Jc=s.Matrix(local['jacobian']).applyfunc(s.Rational);T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);n=T.cols;rows=list(Jc.T.to_DM().rref()[1]);piv=list(Jc.to_DM().rref()[1]);U=s.zeros(45,len(piv));iv=Jc.extract(rows,piv).inv(method='DM')
 for i,c in enumerate(piv):U[c,:]=iv[i,:]
 JU=Jc*U;mon2=list(combinations_with_replacement(range(n),2));mon3=list(combinations_with_replacement(range(n),3));ix={a:i for i,a in enumerate(mon3)};Q=s.zeros(135,len(mon2));terms=[]
 for row,f in enumerate(F):
  ts=[]
  for mon,c in f.items():
   if fixed in mon:continue
   i,j=[cols.index(k) for k in mon]+[44]*(2-len(mon));ts.append((i,j,c))
   for h,(a,b) in enumerate(mon2):Q[row,h]+=c*(T[i,a]*T[j,b]+(T[i,b]*T[j,a] if a!=b else 0))
  terms.append(ts)
 w2=-U*Q[rows,:];A2=Q-JU*Q[rows,:];Q3=s.zeros(135,len(mon3))
 for row,ts in enumerate(terms):
  for i,j,c in ts:
   for a in range(n):
    for h,mon in enumerate(mon2):Q3[row,ix[tuple(sorted((a,)+mon))]]+=c*(T[i,a]*w2[j,h]+T[j,a]*w2[i,h])
 A3=Q3-JU*Q3[rows,:];ar=list(A2.T.to_DM().rref()[1]);ac=list(A2.to_DM().rref()[1]);C=A2[:,ac]*A2.extract(ar,ac).inv(method='DM');B=A3-C*A3[ar,:];br=list(B.T.to_DM().rref()[1]);v=s.symbols(f'v0:{n}');forms=[s.expand(sum(c*s.prod(v[i] for i in mon) for c,mon in zip(row,mons))) for M,mons in [(A2[ar,:],mon2),(B[br,:],mon3)] for row in M.tolist()];sec=s.groebner([f.subs(v[0],0) for f in forms],*v[1:],order='grevlex');assert sec.is_zero_dimensional;gb=s.groebner([f.subs(v[0],1) for f in forms],*v[1:],order='grevlex');assert gb.is_zero_dimensional;leads=[p.LM(order=gb.order).exponents for p in gb.polys];counts=[]
 for deg in range(20):
  count=0
  for inds in combinations_with_replacement(range(n-1),deg):
   e=tuple(inds.count(i) for i in range(n-1))
   if not any(all(a<=b for a,b in zip(l,e)) for l in leads):count+=1
  counts.append(count)
  if count==0:break
 assert sum(counts)==expected;assert all(s.Rational(b['tangent'][0])!=0 for b in branches)
 return dict(jacobian_pivot_rows=rows,jacobian_pivot_columns=piv,quadratic_rank=len(ar),pure_cubic_rank=len(br),initial_forms=list(map(str,forms)),section_groebner=list(map(str,sec)),localized_groebner=list(map(str,gb)),standard_degrees=counts,parameter='v0',multiplicity_upper_bound=expected)

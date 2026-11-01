#!/usr/bin/env python3
"""Exact h=0 ordinary-flag fourth-order core cover: cubic."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WD=Path(sys.argv[1])
import json,pickle,time
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix as DM
st=time.time();cache=pickle.load(open(WD/'flag-eight-flint-projected4.pkl','rb'));data=json.load(open(WD/'flag-global-low-quad.json'));variables=list(map(s.Symbol,cache['variables']));R,*_=ring(variables,s.QQ);K=s.QQ.frac_field(*variables);ms={d:list(cwr(range(9),d))for d in [2,3]};active=[4,5,7,8];rad=[0,1,2,3,6];actcols=[j for j,m in enumerate(ms[2])if all(a in active for a in m)];rr=json.load(open(WD/'flag-b1-universal-quadratic.json'))['pivot_rows']
def mat(key):
 out={};rows=cache[key]
 for i,row in enumerate(rows):
  for j,terms in enumerate(row):
   if terms:out[i,j]=K.field.new(R.from_dict({tuple(e):s.QQ(c)for e,c in terms}))
 return DM.from_dok(out,(len(rows),len(rows[0])),K)
def ext(M,rr,cc):return M.extract(list(rr),list(cc))
Qinv=s.Matrix(data['right_inverse']).applyfunc(s.sympify).to_DM().convert_to(K);L=s.diag(*[1/s.sympify(x)for x in data['row_contents']]).to_DM().convert_to(K);A2,A3=map(mat,['A2','A3']);C2=ext(A2,range(135),actcols)*Qinv*L
assert C2*ext(A2,rr,range(45))==A2
print('global quadratic row relations',time.time()-st,flush=True)
mu=s.Matrix(data['mu']).applyfunc(s.sympify).to_DM().convert_to(K).to_dok();p5={};pairs=list(cwr(active,2));othercols=[j for j,m in enumerate(ms[3])if sum(a in active for a in m)<=1]
for(i,j),c in A3.to_dok().items():
 m=ms[3][j];aa=tuple(a for a in m if a in active);bb=tuple(a for a in m if a not in active)
 if len(aa)==2:
  key=(i,rad.index(bb[0]));p5[key]=p5.get(key,K.zero)+c*mu.get((pairs.index(aa),0),K.zero)
P5=DM.from_dok(p5,(135,5),K);H5=P5-C2*ext(P5,rr,range(5));O=ext(A3,range(135),othercols);assert O-C2*ext(O,rr,range(len(othercols)))==DM.zeros((135,len(othercols)),K)
print('all cubic classes socle*radical',time.time()-st,flush=True)
M=H5.to_Matrix();enc=[[str(x)for x in row]for row in M.tolist()];core={s.Symbol('b2'),s.Symbol('b8')};pure=[i for i in range(135)if any(M[i,:])and set().union(*(x.free_symbols for x in M[i,:]))<=core]
print('pure core rows',pure,time.time()-st,flush=True)
Path(WD/'flag-global-cubic-rows.json').write_text(json.dumps(dict(variables=cache['variables'],socle_radical_rows=enc,pure_core_rows=pure,selected_raw_rows=rr,seconds=time.time()-st),indent=2)+'\n')

import json,time,math
from pathlib import Path
from itertools import combinations
import sympy as s
from sympy.polys.rings import ring
st=time.time();data=json.load(open(WD/'flag-global-cubic-rows.json'));r,h=s.symbols('b2 b8');P,x=ring('b2',s.QQ);rs=data['pure_core_rows'];rows={};hs={};contents={};denoms={}
for i in rs:
 row={};hh={};D=P.one
 for j,v in enumerate(data['socle_radical_rows'][i][2:5]):
  if v=='0':continue
  n,dd=s.fraction(s.sympify(v));np=s.Poly(n,r,h);hp={m[1]for m,c in np.terms()};assert len(hp)==1;hh[j]=hp.pop();assert not dd.has(h);row[j]=(P.from_expr(s.expand(n/h**hh[j])),P.from_expr(dd));D=D.lcm(row[j][1])
 row={j:n*D.exquo(dd)for j,(n,dd)in row.items()};g=next(iter(row.values()))
 for p in row.values():g=g.gcd(p)
 g=g.monic();contents[i]=str(s.factor(g.as_expr()));denoms[i]=str(D.as_expr());rows[i]={j:p.exquo(g)for j,p in row.items()};hs[i]=hh
rw={rs[0]:0};cw={};changed=True
while changed:
 changed=False
 for i,hh in hs.items():
  for j,e in hh.items():
   if i in rw and j not in cw:cw[j]=e-rw[i];changed=True
   if j in cw and i not in rw:rw[i]=e-cw[j];changed=True
assert all(e==rw[i]+cw[j]for i,hh in hs.items()for j,e in hh.items());print('weights',rw,cw,'contents',set(contents.values()),flush=True)
choices=[]
enc=lambda p:[[e[0],str(c)]for e,c in p.items()]
Path(WD/'flag-global-cubic-core-matrix.json').write_text(json.dumps(dict(rows={i:{j:enc(p)for j,p in row.items()}for i,row in rows.items()},row_weights=rw,column_weights=cw,row_contents=contents,row_denominators=denoms,minor_rows=choices,source_rows=rs,seconds=time.time()-st),indent=2)+'\n');print('MINORS READY',len(choices),time.time()-st,flush=True)

import json,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
st=time.time();d=json.load(open(WD/'flag-global-cubic-core-matrix.json'));r,h=s.symbols('b2 b8');rs=d['source_rows'];B=[[[[e,s.Rational(c)]for e,c in d['rows'][str(i)].get(str(j),[])]for j in range(3)]for i in rs]
g=s.expand((r-2)*(r+6)**2*(r+s.Rational(4,25))*(r+s.Rational(91,400)))
print('target',s.factor(g),flush=True);deg=max(e for row in B for p in row for e,c in p)
for bound in range(9):
 unknown=len(rs)*(bound+1);nr=3*(deg+bound+1);D={}
 for i,row in enumerate(B):
  for j,p in enumerate(row):
   for e,c in p:
    for k in range(bound+1):D[j*(deg+bound+1)+e+k,i*(bound+1)+k]=s.QQ.convert(c)
 for j in range(3):
  for(e,),c in s.Poly(g,r).terms():D[j*(deg+bound+1)+e,unknown+j]=s.QQ.convert(c)
 A=DM.from_dok(D,(nr,unknown+3),s.QQ);rr,piv=A.rref();print('bound',bound,'rank',len(piv),'consistent',all(p<unknown for p in piv),time.time()-st,flush=True)
 if all(p<unknown for p in piv):break
else:raise RuntimeError('no inverse at bound4')
R=rr.to_dok();X=s.zeros(3,len(rs))
for row,p in enumerate(piv):
 i,k=divmod(p,bound+1)
 for j in range(3):X[j,i]+=s.QQ.to_sympy(R.get((row,unknown+j),s.QQ.zero))*r**k
BM=s.Matrix([[sum(c*r**e for e,c in p)for p in row]for row in B]);assert (X*BM-g*s.eye(3)).applyfunc(s.expand)==s.zeros(3)
L=s.zeros(3,135)
for j in range(3):
 for k,i in enumerate(rs):
  L[j,i]=s.cancel(h**(-d['column_weights'][str(j)])*X[j,k]/g*s.sympify(d['row_denominators'][str(i)])/s.sympify(d['row_contents'][str(i)])*h**(-d['row_weights'][str(i)]))
enc=lambda M:[[str(v)for v in row]for row in M.tolist()]
Path(WD/'flag-global-cubic-low-inverse.json').write_text(json.dumps(dict(left_inverse=enc(L),primitive_numerator=enc(X),primitive_denominator=str(g),polynomial_bound=bound,seconds=time.time()-st),indent=2)+'\n');print('LOW DEGREE GLOBAL TRUE CUBIC INVERSE PASS',time.time()-st,flush=True)

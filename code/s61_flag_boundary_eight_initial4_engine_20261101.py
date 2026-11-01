#!/usr/bin/env python3
"""Exact fourth-order projection; no lower tail is discarded."""
import sys
from pathlib import Path
WD=Path(sys.argv[1])
import json,pickle,time
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix as DM
st=time.time();cache=pickle.load(open(WD/'flag-eight-flint-projected4.pkl','rb'));v=list(map(s.Symbol,cache['variables']));R,*_=ring(v,s.QQ);K=s.QQ.frac_field(*v);ms={d:list(cwr(range(9),d))for d in [2,3,4]}
def mat(key):
 rows=cache[key];out={}
 for i,row in enumerate(rows):
  for j,terms in enumerate(row):
   if terms:out[i,j]=K.field.new(R.from_dict({tuple(e):s.QQ.convert(s.Rational(c))for e,c in terms}))
 return DM.from_dok(out,(len(rows),len(rows[0])),K)
a2,a3,a4=map(mat,['A2','A3','A4_projected']);print('matrices',a2.shape,a3.shape,a4.shape,time.time()-st,flush=True)
def ext(A,rows,cols):return A.extract(list(rows),list(cols))
def rows(A):return list(A.transpose().rref()[1])
def express(A,B):
 piv=list(B.rref()[1]);C=ext(A,range(A.shape[0]),piv)*ext(B,range(B.shape[0]),piv).inv();assert C*B==A;return C,piv
# Multiplication into degree3, or multiplication followed by G times the
# already certified degree4 quotient projection.
projection=[None if p is None else (p[0],K.from_sympy(s.sympify(p[1])))for p in cache['projection']]
def multiples(A,degree,project=False):
 ix={mon:j for j,mon in enumerate(ms[degree+1])};out={}
 for (i,j),c in A.to_dok().items():
  mon=ms[degree][j]
  for k in range(9):
   dst=ix[tuple(sorted(mon+(k,)))];co=c
   if project:
    data=projection[dst]
    if data is None:continue
    dst,f=data;co*=f
   key=(i*9+k,dst);out[key]=out.get(key,K.zero)+co
 return DM.from_dok(out,(A.shape[0]*9,213 if project else len(ms[degree+1])),K)
core=s.QQ.frac_field(s.Symbol('b2'));r2=[64,65,75,101,102,103,115,121,134];rawb2=ext(a2,r2,range(45));quad=json.load(open(WD/'flag-b1-universal-quadratic.json'));p2=quad['pivot_columns'];qiv=ext(rawb2,range(9),p2).inv();print('9x9 inverse ready',time.time()-st,flush=True);b2=qiv*rawb2;print('quad normalized',time.time()-st,flush=True);f3=qiv*ext(a3,r2,range(165));print('quad tails normalized',time.time()-st,flush=True);c2=ext(a2,range(135),p2);assert c2*b2==a2;h3=a3-c2*f3
print('canonical quadratic generator tails ready',time.time()-st,flush=True)
g3=multiples(b2,2);G=g3.to_dok();active={4,5,7,8};Zcols=[ms[3].index(tuple(sorted((a,8,8))))for a in sorted(active)]
prefix=[j for j,m in enumerate(ms[3])if sum(a in active for a in m)>=2 and j not in Zcols and m not in [(a,8,8)for a in [0,1,2,3,6]]];assert len(prefix)==61
pr=[]
for col in prefix:
 candidates=[i for i in range(81)if G.get((i,col),K.zero)==K.one and all(not G.get((i,k),K.zero)for k in prefix if k!=col)];assert candidates;pr.append(candidates[0])
assert len(set(pr))==61
pb=ext(g3,pr,range(165));rem=g3-ext(g3,range(81),prefix)*pb;small=ext(rem,range(81),Zcols)
# Only a4^2*a8,a5^2*a8,a7^2*a8,a8^3 remain. A four-by-four solve suffices.
srows=list(small.transpose().rref()[1]);assert len(srows)==4;siv=ext(small,srows,range(4)).inv();rb=siv*ext(rem,srows,range(165));qb=(pb-ext(pb,range(61),Zcols)*rb).vstack(rb);piv3=prefix+Zcols
assert ext(qb,range(65),piv3)==DM.eye(65,K)
assert ext(g3,range(81),piv3)*qb==g3
print('canonical cubic products rank65',time.time()-st,flush=True)
ch=ext(h3,range(135),piv3);hr=h3-ch*qb;hrrows=[8,9,27];hrcols=[ms[3].index((i,8,8))for i in [1,2,3]];hb=ext(hr,hrrows,range(165));hiv=ext(hr,hrrows,hrcols).inv();hc=ext(hr,range(135),hrcols)*hiv;assert hc*hb==hr
print('all true cubic relations span3',time.time()-st,flush=True)
assert multiples(qb,3,True).is_zero_matrix
assert multiples(hb,3,True).is_zero_matrix
f4=qiv*ext(a4,r2,range(213));h4=a4-c2*f4;g4=multiples(f3,3,True);pt=ext(g4,pr,range(213));rt=siv*ext(g4-ext(g3,range(81),prefix)*pt,srows,range(213));qt=(pt-ext(pb,range(61),Zcols)*rt).vstack(rt);hc4=h4-ch*qt
extra=(g4-ext(g3,range(81),piv3)*qt).vstack(hc4-hc*ext(hc4,hrrows,range(213)))
print('true projected quartics ready',time.time()-st,flush=True)
t3=mat('T3');ct=ext(t3,[0],piv3);tr=t3-ct*qb;tc=ext(tr,[0],hrcols)*hiv;assert tc*hb==tr
tt=mat('T4_projected')-ct*qt-tc*ext(hc4,hrrows,range(213))
encpoly=lambda p:[[list(m),str(c)]for m,c in p.items()]
encoded=lambda A:{str(i)+','+str(j):[encpoly(c.numer),encpoly(c.denom)]for(i,j),c in A.to_dok().items()if c}
# Every division used here is recorded; the all-point domain also has to
# exclude these factors, not just the final Gaussian pivots.
dens={}
for name,A in [('quad_inverse',qiv),('quad',b2),('quad_tail3',f3),('quad_tail4',f4),('cubic_inverse',siv),('cubic_products',qb),('cubic_tail',qt),('true_cubic_inverse',hiv),('true_cubic_relations',hc),('t_cubic_coefficients',ct),('t_true_cubic_coefficients',tc)]:
 ds={str(K.to_sympy(K.field.new(c.denom)))for c in A.to_dok().values()};dens[name]=sorted(ds)
payload=dict(variables=cache['variables'],shape=extra.shape,rows=encoded(extra),corrected_t4=encoded(tt),denominators=dens,projection_multiplier=cache['projection_multiplier'],scaling_denominators=dict(ud=cache['ud'],td=cache['td']),quadratic_rank=9,cubic_product_rank=65,cubic_new_rank=3,quartic_lower_rank=282,quartic_quotient_dimension=213,t3_in_true_cubic_ideal=True,seconds=time.time()-st)
with open(WD/'flag-eight-initial4-field.json','w')as f:json.dump(payload,f);f.write('\n')
print('TRUE QUARTIC AND CORRECTED T4 READY',time.time()-st,flush=True)

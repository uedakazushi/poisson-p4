"""Independent QQ(theta) and Laurent-polynomial arithmetic for T38 review.

theta^4=50theta^2-305. No SymPy algebraic-field or producer arithmetic.
"""
from fractions import Fraction as F
from functools import lru_cache
import itertools

class Q4:
 __slots__=('c',)
 def __init__(self,v=0):
  if isinstance(v,Q4):self.c=v.c
  elif isinstance(v,(list,tuple)):
   assert len(v)<=4;self.c=tuple(F(x)for x in v)+(F(0),)*(4-len(v))
  else:self.c=(F(v),F(0),F(0),F(0))
 def __hash__(self):return hash(self.c)
 def __bool__(self):return any(self.c)
 def __eq__(self,o):return self.c==Q4(o).c
 def __add__(self,o):o=Q4(o);return Q4(tuple(a+b for a,b in zip(self.c,o.c)))
 __radd__=__add__
 def __neg__(self):return Q4(tuple(-x for x in self.c))
 def __sub__(self,o):return self+-Q4(o)
 def __rsub__(self,o):return Q4(o)+-self
 def __mul__(self,o):
  o=Q4(o)
  if not self or not o:return ZERO
  c=[F(0)]*7
  for i,a in enumerate(self.c):
   for j,b in enumerate(o.c):c[i+j]+=a*b
  for k in range(6,3,-1):c[k-2]+=50*c[k];c[k-4]-=305*c[k]
  return Q4(c[:4])
 __rmul__=__mul__
 def __truediv__(self,o):
  o=Q4(o);assert o
  if not any(o.c[1:]):return Q4(tuple(x/o.c[0]for x in self.c))
  return self*inverse(o)
 def __rtruediv__(self,o):return Q4(o)/self
 def __pow__(self,n):
  if n<0:return inverse(self)**(-n)
  p=ONE;q=self
  while n:
   if n&1:p=p*q
   q=q*q;n//=2
  return p
 def encode(self):return list(map(str,self.c))
ZERO=Q4(0);ONE=Q4(1);THETA=Q4([0,1])

def rational_det(A):
 B=[[F(x)for x in row]for row in A];n=len(B);d=F(1)
 for j in range(n):
  i=next((i for i in range(j,n)if B[i][j]),None)
  if i is None:return F(0)
  if i!=j:B[i],B[j]=B[j],B[i];d=-d
  p=B[j][j];d*=p
  for i in range(j+1,n):
   c=B[i][j]/p
   for k in range(j+1,n):B[i][k]-=c*B[j][k]
 return d
def multiplication_matrix(q):return [[(q*THETA**j).c[i]for j in range(4)]for i in range(4)]
def norm(q):return rational_det(multiplication_matrix(q))
@lru_cache(maxsize=None)
def inverse(q):
 assert q;M=[row+[F(i==0)]for i,row in enumerate(multiplication_matrix(q))]
 for j in range(4):
  i=next(i for i in range(j,4)if M[i][j]);M[i],M[j]=M[j],M[i];p=M[j][j];M[j]=[x/p for x in M[j]]
  for i in range(4):
   if i!=j:
    c=M[i][j];M[i]=[a-c*b for a,b in zip(M[i],M[j])]
 answer=Q4([row[4]for row in M]);assert q*answer==ONE;return answer
def determinant(A):
 B=[[Q4(x)for x in row]for row in A];d=ONE;n=len(B)
 for j in range(n):
  i=next((i for i in range(j,n)if B[i][j]),None)
  if i is None:return ZERO
  if i!=j:B[i],B[j]=B[j],B[i];d=-d
  p=B[j][j];d=d*p
  for i in range(j+1,n):
   c=B[i][j]/p
   for k in range(j+1,n):B[i][k]=B[i][k]-c*B[j][k]
 return d

# Sparse Laurent polynomials in t, coefficients in the field above.
def pc(c=1,power=0):c=Q4(c);return {power:c}if c else{}
def pa(p,q):
 d=p.copy()
 for k,c in q.items():
  v=d.get(k,ZERO)+c
  if v:d[k]=v
  else:d.pop(k,None)
 return d
def pn(p):return {k:-c for k,c in p.items()}
def ps(p,c):c=Q4(c);return {k:v*c for k,v in p.items()if v*c}
def pm(p,q):
 d={}
 for k,a in p.items():
  for l,b in q.items():
   e=k+l;c=d.get(e,ZERO)+a*b
   if c:d[e]=c
   else:d.pop(e,None)
 return d
def mzero(n=5):return [[{}for j in range(n)]for i in range(n)]
def mid(n=5):return [[pc(i==j)for j in range(n)]for i in range(n)]
def ma(A,B):return [[pa(a,b)for a,b in zip(r,s)]for r,s in zip(A,B)]
def ms(A,c):return [[ps(a,c)for a in row]for row in A]
def mm(A,B):
 return [[sum_poly(pm(A[i][k],B[k][j])for k in range(len(B)))for j in range(len(B[0]))]for i in range(len(A))]
def sum_poly(polys):
 d={}
 for p in polys:d=pa(d,p)
 return d
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2));BASIS=list(itertools.product(PAIRS,MONS));POS={x:i for i,x in enumerate(BASIS)}
def wedge(A,B):
 result=[{}for _ in BASIS]
 for i,j in PAIRS:
  for u in range(5):
   for v in range(5):
    c=pa(pm(A[i][u],B[j][v]),pn(pm(B[i][u],A[j][v])))
    at=POS[(i,j),tuple(sorted((u,v)))];result[at]=pa(result[at],c)
 return result
def physical_checks(values,strict=True):
 tensor={pair:{}for pair in PAIRS}
 for(pair,mon),c in zip(BASIS,values):
  if c:tensor[pair][mon]=c
 def component(i,j):
  if i==j:return{}
  return tensor[i,j]if i<j else{m:pn(c)for m,c in tensor[j,i].items()}
 def derivative(p,k):
  d={}
  for mon,c in p.items():
   for i,v in enumerate(mon):
    if v==k:
     m=mon[:i]+mon[i+1:];d[m]=pa(d.get(m,{}),c)
  return {m:c for m,c in d.items()if c}
 divergence=[]
 for j in range(5):
  d={}
  for i in range(5):
   for mon,c in derivative(component(i,j),i).items():d[mon]=pa(d.get(mon,{}),c)
  divergence.append({m:c for m,c in d.items()if c})
 if not strict:return divergence,[]
 jacobi=[]
 for i,j,k in itertools.combinations(range(5),3):
  d={}
  for a,b,c in [(i,j,k),(j,k,i),(k,i,j)]:
   for l in range(5):
    for u,cu in component(a,l).items():
     for v,cv in derivative(component(b,c),l).items():
      mon=tuple(sorted(u+v));d[mon]=pa(d.get(mon,{}),pm(cu,cv))
  jacobi.append({m:c for m,c in d.items()if c})
 return divergence,jacobi
def push(values,S,Sinv):
 answer=[{}for _ in BASIS]
 for((a,b),(p,q)),c in zip(BASIS,values):
  if not c:continue
  for i,j in PAIRS:
   vector=pa(pm(S[i][a],S[j][b]),pn(pm(S[i][b],S[j][a])))
   if not vector:continue
   for u in range(5):
    for v in range(5):
     d=pm(pm(pm(c,vector),Sinv[p][u]),Sinv[q][v]);at=POS[(i,j),tuple(sorted((u,v)))];answer[at]=pa(answer[at],d)
 return answer

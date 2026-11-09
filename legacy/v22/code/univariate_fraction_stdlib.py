"""Exact sparse polynomial arithmetic over Q(lambda), standard library only."""
import ast
from fractions import Fraction as Q
import audit_certificates_stdlib as A

def trim(p):
 p=list(map(Q,p))
 while p and not p[-1]:p.pop()
 return tuple(p)
def ua(p,q):
 r=[Q(0)]*max(len(p),len(q))
 for i,v in enumerate(p):r[i]+=v
 for i,v in enumerate(q):r[i]+=v
 return trim(r)
def um(p,q):
 if not p or not q:return ()
 r=[Q(0)]*(len(p)+len(q)-1)
 for i,a in enumerate(p):
  for j,b in enumerate(q):r[i+j]+=a*b
 return trim(r)
def ud(p,q):
 if not q:raise ZeroDivisionError
 p=list(trim(p));r=[Q(0)]*max(0,len(p)-len(q)+1)
 while len(p)>=len(q):
  k=len(p)-len(q);c=p[-1]/q[-1];r[k]+=c
  for i,v in enumerate(q):p[i+k]-=c*v
  while p and not p[-1]:p.pop()
 return trim(r),trim(p)
def uq(p,q):
 a,b=ud(p,q)
 if b:raise ArithmeticError('inexact univariate division')
 return a
def ug(p,q):
 while q:p,q=q,ud(p,q)[1]
 return tuple(v/p[-1] for v in p) if p else ()
def uone():return (Q(1),)
def uni(p):
 out={}
 for m,c in p.items():
  if any(e for i,e in enumerate(m) if i!=A.LI):raise ValueError('denominator not univariate')
  out[m[A.LI]]=c
 r=[Q(0)]*(1+max(out,default=-1))
 for k,c in out.items():r[k]=c
 return trim(r)
def groups(p):
 g={}
 for m,c in p.items():
  mm=list(m);k=mm[A.LI];mm[A.LI]=0;mm=tuple(mm);g.setdefault(mm,{})[k]=c
 out={}
 for m,v in g.items():
  z=[Q(0)]*(1+max(v));
  for k,c in v.items():z[k]=c
  out[m]=trim(z)
 return out
def ungroup(g):
 out={}
 for m,u in g.items():
  for k,c in enumerate(u):
   if c:mm=list(m);mm[A.LI]=k;out[tuple(mm)]=c
 return out
def muluni(p,u):return ungroup({m:um(v,u) for m,v in groups(p).items()})
MODULUS=None
def set_modulus(p=None):
 global MODULUS
 MODULUS=trim(p) if p is not None else None
def uinv(a,m):
 r0,r1=m,a;s0,s1=(),(Q(1),)
 while r1:
  q,r=ud(r0,r1);r0,r1=r1,r;s0,s1=s1,ua(s0,tuple(-v for v in um(q,s1)))
 if len(r0)!=1:raise ZeroDivisionError('nonunit in coefficient quotient')
 return ud(tuple(v/r0[0] for v in s0),m)[1]
class G:
 def __init__(self,p=None,d=(Q(1),)):
  p={m:Q(c) for m,c in (p or {}).items() if c};d=trim(d)
  if not d:raise ZeroDivisionError
  if not p:self.p={};self.d=uone();return
  if MODULUS is not None:
   if d==uone() and all(m[A.LI]<len(MODULUS)-1 for m in p):self.p=p;self.d=uone();return
   d=ud(d,MODULUS)[1];inv=uinv(d,MODULUS)
   p=ungroup({m:ud(um(u,inv),MODULUS)[1] for m,u in groups(p).items()});self.p=p;self.d=uone();return
  lead=d[-1];p={m:c/lead for m,c in p.items()};d=tuple(c/lead for c in d)
  if len(d)>1:
   gs=groups(p);h=d
   for u in gs.values():
    h=ug(h,u)
    if len(h)==1:break
   if len(h)>1:p=ungroup({m:uq(u,h) for m,u in gs.items()});d=uq(d,h)
  self.p=p;self.d=d
 @staticmethod
 def c(x):return G({A.ZERO:Q(x)})
 def __neg__(self):return G({m:-c for m,c in self.p.items()},self.d)
 def __add__(self,o):
  if not isinstance(o,G):o=G.c(o)
  h=ug(self.d,o.d);l=uq(o.d,h);r=uq(self.d,h)
  return G(A.addp(muluni(self.p,l),muluni(o.p,r)),um(self.d,l))
 __radd__=__add__
 def __sub__(self,o):return self+-o
 def __rsub__(self,o):return -self+o
 def __mul__(self,o):
  if not isinstance(o,G):o=G.c(o)
  return G(A.mulp(self.p,o.p),um(self.d,o.d))
 __rmul__=__mul__
 def __truediv__(self,o):
  if not isinstance(o,G):o=G.c(o)
  return G(muluni(self.p,o.d),um(self.d,uni(o.p)))
 def __pow__(self,k):
  if k<0:return G.c(1)/(self**(-k))
  r=G.c(1);v=self
  while k:
   if k%2:r=r*v
   v=v*v;k//=2
  return r
 def iszero(self):return not self.p
 def key(self):return self.d,tuple(sorted(self.p.items()))
def var(name):
 m=list(A.ZERO);m[A.I[name]]=1;return G({tuple(m):Q(1)})
SY={v:var(v) for v in A.V};CACHE={}
def parse(s,env=None,alias=None):
 if s not in CACHE:CACHE[s]=ast.parse(s,mode='eval').body
 env=env or {};alias=alias or {}
 def ev(n):
  if isinstance(n,ast.Constant) and isinstance(n.value,int):return G.c(n.value)
  if isinstance(n,ast.Name):
   if n.id in env:return env[n.id]
   return alias[n.id] if n.id in alias else SY[n.id]
  if isinstance(n,ast.UnaryOp):return -ev(n.operand) if isinstance(n.op,ast.USub) else ev(n.operand)
  if isinstance(n,ast.BinOp):
   a=ev(n.left)
   if isinstance(n.op,ast.Pow):
    k=ast.literal_eval(n.right);assert isinstance(k,int);return a**k
   b=ev(n.right)
   if isinstance(n.op,ast.Add):return a+b
   if isinstance(n.op,ast.Sub):return a-b
   if isinstance(n.op,ast.Mult):return a*b
   if isinstance(n.op,ast.Div):return a/b
  raise ValueError(ast.dump(n))
 return ev(CACHE[s])
def keys(es):return {e.key() for e in es if not e.iszero()}
def supported_denominator(d,allowed):
 used=[];d=trim(d)
 for label,f in allowed:
  found=False
  while len(d)>1:
   q,r=ud(d,f)
   if r:break
   d=q;found=True
  if found:used.append(label)
 if len(d)>1:raise ValueError('uncovered denominator '+repr(d))
 return used

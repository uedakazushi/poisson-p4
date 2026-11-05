"""Prepared SOURCE-only packed5 profile engine, H-sliced canonical hash.

No coefficient matrix values: only existing-leaf support, QQ scaling, positive
integer majorants, and literal expression bindings. Future VALUE counts remain
five-coordinate coefficient maps; this engine's support packing does not apply
to their scalar ownership.
"""
from fractions import Fraction as Q
from math import gcd
import hashlib,json,re
WIDTH=256
STRIDES=tuple(256**j for j in range(5))
BIAS=64*sum(STRIDES)
def rawindex(e):return sum(e[j]*STRIDES[j]for j in range(5))
def pack(e):
 assert len(e)==5 and all(-64<=v<=191 for v in e)
 return BIAS+rawindex(e)
def unpack(v):
 assert 0<=v<256**5
 return tuple((v//STRIDES[j])%256-64 for j in range(5))
def count_numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,str):return 3 if re.fullmatch(r'[+-]?[0-9]+(?:/[0-9]+)?',v)else 0
 if isinstance(v,list):return sum(map(count_numeric,v))
 if isinstance(v,dict):return sum(map(count_numeric,v.values()))
 return 0
class Engine:
 def __init__(self,metadata_base=25000,guard=150000,checkpoint=None):
  self.base=metadata_base;self.guard=guard;self.live={};self.nodes={};self.ops=[]
  self.peak=metadata_base;self.meta=0;self.context={};self.checkpoint=checkpoint
  self.pair_count=0;self.hash_H_alias_peak=0
 def own(self,extra=256):
  n=self.base+self.meta+sum(len(x['S'])for x in self.live.values())+extra
  self.peak=max(self.peak,n)
  if n>self.guard:
   if self.checkpoint:self.checkpoint(self,n)
   raise RuntimeError('prospective SOURCE owned guard '+str(n))
  return n
 def drop(self,*names):
  for k in names:self.live.pop(k,None)
 @staticmethod
 def _lex_alias(S):
  a=list(S)
  def greater(x,y):
   for stride in STRIDES:
    u=x//stride%256;v=y//stride%256
    if u!=v:return u>v
   return False
  def down(root,end):
   while 2*root+1<=end:
    child=2*root+1
    if child+1<=end and greater(a[child+1],a[child]):child+=1
    if not greater(a[child],a[root]):break
    a[root],a[child]=a[child],a[root];root=child
  for i in range((len(a)-2)//2,-1,-1):down(i,len(a)-1)
  for end in range(len(a)-1,0,-1):a[0],a[end]=a[end],a[0];down(0,end-1)
  return a
 def supporthash(self,S):
  # H is the first coordinate in decoded tuple lex order. Each slice is
  # hashed then released; no full support alias list is constructed.
  self.own(768);counts={}
  for e in S:
   h=e%256;counts[h]=counts.get(h,0)+1
  nmax=max(counts.values(),default=0);self.own(nmax+768)
  self.hash_H_alias_peak=max(self.hash_H_alias_peak,nmax)
  hh=hashlib.sha256();hh.update(b'[');first=True
  for h in sorted(counts):
   subset=(e for e in S if e%256==h)
   # Construct exactly the forecast number of alias refs, then heapsort
   # that same list in _lex_alias. No simultaneous alias copy is allowed.
   a=list(subset);assert len(a)==counts[h]
   def greater(x,y):
    for stride in STRIDES:
     u=x//stride%256;v=y//stride%256
     if u!=v:return u>v
    return False
   def down(root,end):
    while 2*root+1<=end:
     child=2*root+1
     if child+1<=end and greater(a[child+1],a[child]):child+=1
     if not greater(a[child],a[root]):break
     a[root],a[child]=a[child],a[root];root=child
   for i in range((len(a)-2)//2,-1,-1):down(i,len(a)-1)
   for end in range(len(a)-1,0,-1):a[0],a[end]=a[end],a[0];down(0,end-1)
   for e in a:
    if not first:hh.update(b',')
    first=False;hh.update(json.dumps(list(unpack(e)),separators=(',',':')).encode())
   del a
  hh.update(b']\n');return hh.hexdigest(),counts
 def save(self,name,S,q,M,op):
  assert name not in self.nodes
  self.live[name]={'S':S,'q':Q(q),'M':int(M)};self.own()
  if S:box=[[min(unpack(e)[j]for e in S),max(unpack(e)[j]for e in S)]for j in range(5)]
  else:box=[[0,0]]*5
  h,counts=self.supporthash(S)
  node={'terms_upper':len(S),'known_QQ_content':str(Q(q)),'primitive_integer_abs_upper':str(M),'support_sha256':h,'exponent_box5':box,'H_slice_sizes':[[k-64,v]for k,v in sorted(counts.items())]}
  self.nodes[name]=node;self.ops.append({'output':name,**op})
  self.meta+=count_numeric(node)+count_numeric(op)+12
  self.own();return name
 def leaf(self,name,p,source=None,scalar=Q(1),shift=(0,0,0,0,0)):
  self.context={'op':'leaf','output':name};self.own(5*len(p['terms'])+768)
  S=set()
  for t in p['terms']:
   e=tuple(p['shift'][j]+int(t[j+1])+shift[j]for j in range(5));v=pack(e)
   if v not in S:self.own(len(S)+1+768);S.add(v)
  return self.save(name,S,Q(p['content'])*scalar,max([abs(int(t[0]))for t in p['terms']]+[0]),{'op':'input','source':source})
 def copy(self,name,A,scalar=Q(1),shift=(0,0,0,0,0)):
  x=self.live[A];self.context={'op':'copy','output':name};box=self.nodes[A]['exponent_box5']
  assert all(-64<=lo+shift[j]<=hi+shift[j]<=191 for j,(lo,hi)in enumerate(box))
  self.own(len(x['S'])+768);d=rawindex(shift)
  S={e+d for e in x['S']}
  return self.save(name,S,x['q']*scalar,x['M'],{'op':'copy','A':A,'scalar':str(scalar),'shift5':list(shift)})
 def mul(self,name,A,B,scalar=Q(1),shift=(0,0,0,0,0)):
  x=self.live[A];y=self.live[B];self.context={'op':'mul','output':name,'A':A,'B':B}
  bx=self.nodes[A]['exponent_box5'];by=self.nodes[B]['exponent_box5']
  box=[[bx[j][0]+by[j][0]+shift[j],bx[j][1]+by[j][1]+shift[j]]for j in range(5)]
  assert all(-64<=lo<=hi<=191 for lo,hi in box)
  d=rawindex(shift)-BIAS;S=set()
  for e in x['S']:
   for f in y['S']:
    g=e+f+d
    if g not in S:self.own(len(S)+1+768);S.add(g)
  n=min(len(x['S']),len(y['S']));self.pair_count+=len(x['S'])*len(y['S'])
  return self.save(name,S,x['q']*y['q']*scalar,n*x['M']*y['M'],{'op':'multiply','A':A,'B':B,'scalar':str(scalar),'shift5':list(shift),'pair_BOX_before_Minkowski':box,'signed_bias_correction':str(d)})
 def add(self,name,A,B,scalar=Q(1),consume=False):
  x=self.live[A];y=self.live[B];self.context={'op':'add','output':name,'A':A,'B':B}
  aa=x['q'];bb=y['q']*scalar;d=aa.denominator*bb.denominator//gcd(aa.denominator,bb.denominator)
  a=int(aa*d);b=int(bb*d);g=gcd(abs(a),abs(b));q=Q(g,d)if g else Q(0);wa=a//g if g else 0;wb=b//g if g else 0
  M=abs(wa)*x['M']+abs(wb)*y['M']
  if consume:
   S=x['S'];T=y['S'];self.drop(A,B)
   while T:
    self.own(len(S)+len(T)+768);S.add(T.pop())
  else:
   self.own(len(x['S'])+768);S=set(x['S'])
   for e in y['S']:
    if e not in S:self.own(len(S)+1+768);S.add(e)
  return self.save(name,S,q,M,{'op':'add','A':A,'B':B,'scalarB':str(scalar),'integer_weight_A':str(wa),'integer_weight_B':str(wb),'consume_lastuse':consume})

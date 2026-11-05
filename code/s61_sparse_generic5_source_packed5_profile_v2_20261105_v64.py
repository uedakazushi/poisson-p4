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
  q=Q(q);M=int(M)
  assert max(abs(q.numerator).bit_length(),q.denominator.bit_length(),abs(M).bit_length())<=4096
  self.live[name]={'S':S,'q':q,'M':M};self.own()
  # One linear pass obtains EXACT boxes and H-counts. Intermediates do not
  # allocate/sort a hash alias; boundary_hash does that only when requested.
  low=[191]*5;high=[-64]*5;counts={}
  self.own(1024)
  for e in S:
   for j in range(5):
    v=(e//STRIDES[j])%256-64;low[j]=min(low[j],v);high[j]=max(high[j],v)
   h=e%256;counts[h]=counts.get(h,0)+1
  box=[[low[j],high[j]]for j in range(5)]if S else[[0,0]]*5
  node={'terms_upper':len(S),'known_QQ_content':str(q),'primitive_integer_abs_upper':str(M),'support_sha256':None,'exponent_box5':box,'H_slice_sizes':[[k-64,v]for k,v in sorted(counts.items())],'hash_policy':'literal_DAG_and_exact_card_BOX_Hcounts; boundary canonical hash only'}
  self.nodes[name]=node;self.ops.append({'output':name,**op})
  self.meta+=count_numeric(node)+count_numeric(op)+12
  self.own();return name
 def boundary_hash(self,name,expected=None):
  h,counts=self.supporthash(self.live[name]['S'])
  assert [[k-64,v]for k,v in sorted(counts.items())]==self.nodes[name]['H_slice_sizes']
  if expected is not None:assert h==expected
  self.nodes[name]['support_sha256']=h;self.nodes[name]['hash_policy']='verified canonical decoded tuplelex at required boundary'
  return h
 def snapshot(self,path):
  # Every entry is streamed twice in the SAME unmodified set iteration order:
  # writer then exact readback. No whole-live copy or alias list is made.
  from pathlib import Path
  path=Path(path);assert not path.exists()
  names=sorted(self.live);total=sum(len(self.live[k]['S'])for k in names)
  assert len(names)<=64 and total<=125000 and all(len(k)<=80 for k in names)
  self.own(8192)
  header={'format':'Q5SUPP1','packed_base':256,'bias_coordinate':64,'maps':[]}
  for k in names:
   v=self.live[k];q=v['q'];M=v['M'];assert max(abs(q.numerator).bit_length(),q.denominator.bit_length(),abs(M).bit_length())<=4096
   header['maps'].append({'id':k,'q':str(q),'M':str(M),'terms':len(v['S'])})
  b=json.dumps(header,separators=(',',':')).encode();assert len(b)<=4096*len(names)+256
  cap=12+4096*64+256+5*125000;tmp=path.with_suffix(path.suffix+'.partial');assert not tmp.exists()
  with tmp.open('wb')as f:
   f.write(b'Q5SUPP1\n');f.write(len(b).to_bytes(4,'big'));f.write(b)
   for k in names:
    for e in self.live[k]['S']:assert 0<=e<2**40;f.write(e.to_bytes(5,'big'))
  assert tmp.stat().st_size<=cap
  with tmp.open('rb')as f:
   assert f.read(8)==b'Q5SUPP1\n';n=int.from_bytes(f.read(4),'big');assert n==len(b)and f.read(n)==b
   for k in names:
    for e in self.live[k]['S']:assert f.read(5)==e.to_bytes(5,'big')
   assert not f.read(1)
  tmp.rename(path);hh=hashlib.sha256()
  with path.open('rb')as f:
   for b in iter(lambda:f.read(65536),b''):hh.update(b)
  return{'file':path.name,'sha256':hh.hexdigest(),'bytes':path.stat().st_size,'maps':len(names),'total_keys':total,'all_keys_exact_stream_roundtrip':True,'no_coefficient_VALUE':True}
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

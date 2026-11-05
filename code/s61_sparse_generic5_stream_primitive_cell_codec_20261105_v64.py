"""Canonical Q5CELL1 primitive cells with alias-key heapsort and stream readback."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import hashlib
MAGIC=b'Q5CELL1\n'
def uint(n):
 assert isinstance(n,int)and n>=0
 b=bytearray()
 while n>=128:b.append((n&127)|128);n>>=7
 b.append(n);return bytes(b)
def sint(n):
 assert isinstance(n,int);m=abs(n);b=m.to_bytes((m.bit_length()+7)//8,'big');return uint(2*len(b)+(n<0))+b
def _greater(a,b):
 x=sum(a);y=sum(b);return x>y or(x==y and a>b)
def descending_alias_keys(values):
 # Only this list of existing key references is O(n). No decoration or key copy.
 k=list(values);n=len(k)
 def down(start,end):
  root=start
  while 2*root+1<=end:
   child=2*root+1
   if child+1<=end and _greater(k[child+1],k[child]):child+=1
   if not _greater(k[child],k[root]):break
   k[root],k[child]=k[child],k[root];root=child
 for start in range((n-2)//2,-1,-1):down(start,n-1)
 for end in range(n-1,0,-1):k[0],k[end]=k[end],k[0];down(0,end-1)
 k.reverse();return k
def _vh_start(n):
 h=hashlib.sha256();h.update((str(n)+'\n').encode());return h
def _vh_term(h,q,e):h.update((str(q)+' '+' '.join(map(str,e))+'\n').encode())
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def write_cell(path,q,values,shift=None):
 path=Path(path);assert not path.exists();tmp=path.with_suffix(path.suffix+'.partial');assert not tmp.exists();q=Q(q);n=len(values);h=_vh_start(n)
 if not n:
  with tmp.open('wb')as f:f.write(MAGIC+b'\x00')
  tmp.rename(path);return{'terms':0,'G':1,'content':'0','shift':[0]*5,'value_sha256':h.hexdigest(),'sha256':sha(path),'bytes':path.stat().st_size}
 assert q and all(isinstance(c,int)and c for c in values.values());G=0
 for c in values.values():G=gcd(G,abs(c))
 if shift is None:shift=[min(e[j]for e in values)for j in range(5)]
 assert len(shift)==5 and all(len(e)==5 and all(e[j]>=shift[j]for j in range(5))for e in values)
 keys=descending_alias_keys(values);qq=q*G
 with tmp.open('wb')as f:
  f.write(MAGIC+b'\x01');f.write(sint(qq.numerator));f.write(sint(qq.denominator))
  for x in shift:f.write(uint(2*x if x>=0 else-2*x-1))
  f.write(uint(n));previous=None
  for e in keys:
   assert previous is None or _greater(previous,e);previous=e;c=values[e];assert c%G==0;f.write(sint(c//G))
   for j in range(5):f.write(uint(e[j]-shift[j]))
   _vh_term(h,q*c,e)
 del keys;tmp.rename(path);return{'terms':n,'G':G,'content':str(qq),'shift':list(shift),'value_sha256':h.hexdigest(),'sha256':sha(path),'bytes':path.stat().st_size}
def validate_cell(path,q,values,record):
 path=Path(path);q=Q(q);G=record['G'];h=_vh_start(len(values))
 with path.open('rb')as f:
  assert f.read(8)==MAGIC;flag=f.read(1)
  if not values:
   assert flag==b'\x00'and not f.read(1)and record['terms']==0
  else:
   assert flag==b'\x01'
   def u():
    n=0;s=0;b=bytearray()
    while True:
     c=f.read(1);assert c;b.extend(c);v=c[0];n|=(v&127)<<s
     if not v&128:break
     s+=7;assert s<100000
    assert bytes(b)==uint(n);return n
   def i():
    tag=u();b=f.read(tag//2);assert len(b)==tag//2 and(not b or b[0]);n=int.from_bytes(b,'big');assert n or not tag&1;v=-n if tag&1 else n;assert sint(v)==uint(tag)+b;return v
   num=i();den=i();assert den>0 and Q(num,den)==q*G and(Q(num,den).numerator,Q(num,den).denominator)==(num,den);shift=[]
   for _ in range(5):v=u();shift.append(v//2 if not v&1 else-(v//2)-1)
   assert shift==record['shift'];n=u();assert n==len(values)==record['terms'];previous=None;gc=0
   for _ in range(n):
    c=i();assert c;e=tuple(shift[j]+u()for j in range(5));assert previous is None or _greater(previous,e);previous=e;assert values.get(e)==c*G;gc=gcd(gc,abs(c));_vh_term(h,q*c*G,e)
   assert gc==1 and not f.read(1)
 assert h.hexdigest()==record['value_sha256']and sha(path)==record['sha256'];return True

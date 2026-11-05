"""Prepared SOURCE-only immutable support views and fused union-product."""
from fractions import Fraction as Q
from math import gcd
from s61_sparse_generic5_source_packed5_profile_v2_20261105_v64 import Engine as Base,STRIDES,BIAS,rawindex,unpack,count_numeric
class Engine(Base):
 def own(self,extra=256):
  assert len(self.live)<=64
  seen={id(v['S']):len(v['S'])for v in self.live.values()}
  n=self.base+self.meta+sum(seen.values())+4*len(self.live)+4*len(seen)+4096+extra
  self.peak=max(self.peak,n)
  if n>self.guard:
   if self.checkpoint:self.checkpoint(self,n)
   raise RuntimeError('prospective SOURCE unique-backing guard '+str(n))
  return n
 def owners(self,S):return sum(v['S']is S for v in self.live.values())
 @staticmethod
 def fields(q,M):
  assert abs(Q(q).numerator).bit_length()<=128 and Q(q).denominator.bit_length()<=128
  assert int(M)>=0 and int(M).bit_length()<=512
 def _record(self,name,S,q,M,off,op,box=None,hs=None):
  assert name not in self.nodes
  q=Q(q);M=int(M);self.fields(q,M);self.own(1024)
  self.live[name]={'S':S,'q':q,'M':M,'off':off};self.own(1024)
  if box is None:
   lo=[191]*5;hi=[-64]*5;counts={}
   for e in S:
    true=e+off;assert 0<=true<2**40
    for j in range(5):
     v=true//STRIDES[j]%256-64;lo[j]=min(lo[j],v);hi[j]=max(hi[j],v)
    h=true%256-64;counts[h]=counts.get(h,0)+1
   box=[[lo[j],hi[j]]for j in range(5)]if S else[[0,0]]*5;hs=sorted(counts.items())
  n={'terms_upper':len(S),'known_QQ_content':str(q),'primitive_integer_abs_upper':str(M),'support_sha256':None,'exponent_box5':box,'H_slice_sizes':[[h,c]for h,c in hs],'view_offset_integer':off,'support_profile':'exact materialized union or immutable view'}
  self.nodes[name]=n;self.ops.append({'output':name,**op});self.meta+=count_numeric(n)+count_numeric(op)+12;self.own();return name
 def save(self,name,S,q,M,op):return self._record(name,S,q,M,0,op)
 def copy(self,name,A,scalar=Q(1),shift=(0,0,0,0,0)):
  x=self.live[A];self.context={'op':'SOURCE_immutable_shift_view','output':name,'A':A}
  box=[[a+shift[j],b+shift[j]]for j,(a,b)in enumerate(self.nodes[A]['exponent_box5'])]
  assert all(-64<=a<=b<=191 for a,b in box)
  hs=[(h+shift[0],n)for h,n in self.nodes[A]['H_slice_sizes']]
  if not x['S']:box=[[0,0]]*5
  return self._record(name,x['S'],x['q']*scalar,x['M'],x['off']+rawindex(shift),{'op':'immutable_SOURCE_view','A':A,'scalar':str(scalar),'shift5':list(shift),'no_set_copy':True},box,hs)
 def mul(self,name,A,B,scalar=Q(1),shift=(0,0,0,0,0)):
  x=self.live[A];y=self.live[B];self.context={'op':'multiply','output':name,'A':A,'B':B}
  bx=self.nodes[A]['exponent_box5'];by=self.nodes[B]['exponent_box5']
  box=[[bx[j][0]+by[j][0]+shift[j],bx[j][1]+by[j][1]+shift[j]]for j in range(5)]
  assert all(-64<=a<=b<=191 for a,b in box)
  q=x['q']*y['q']*scalar;M=min(len(x['S']),len(y['S']))*x['M']*y['M'];self.fields(q,M)
  S=set();d=x['off']+y['off']+rawindex(shift)-BIAS
  for e in x['S']:
   for f in y['S']:
    v=e+f+d
    if v not in S:self.own(len(S)+1+1024);S.add(v)
  self.pair_count+=len(x['S'])*len(y['S'])
  return self.save(name,S,q,M,{'op':'multiply','A':A,'B':B,'scalar':str(scalar),'shift5':list(shift),'pair_BOX_before_Minkowski':box,'true_bias_offset':str(d)})
 def add(self,*args,**kwargs):raise RuntimeError('Fixed resume uses SOURCE fused_union_product; standalone add not implemented')
 def fused_union_product(self,name,C,A,B,dotname):
  x=self.live[A];y=self.live[B];self.context={'op':'SOURCE_fused_union_product','output':name,'C':C,'A':A,'B':B,'literal_dot':dotname}
  self.own(2048)
  bx=self.nodes[A]['exponent_box5'];by=self.nodes[B]['exponent_box5']
  db=[[bx[j][0]+by[j][0],bx[j][1]+by[j][1]]for j in range(5)];assert all(-64<=a<=b<=191 for a,b in db)
  dq=x['q']*y['q'];dM=min(len(x['S']),len(y['S']))*x['M']*y['M'];self.fields(dq,dM)
  hupper={}
  for a,na in self.nodes[A]['H_slice_sizes']:
   for b,nb in self.nodes[B]['H_slice_sizes']:hupper[a+b]=hupper.get(a+b,0)+na*nb
  if C is None:
   S=set();q=dq;M=dM;wa=0;wb=1;lo=[191]*5;hi=[-64]*5;counts={}
   self.live[name]={'S':S,'q':q,'M':M,'off':0}
  else:
   c=self.live[C];assert C not in(A,B)and c['S']is not x['S']and c['S']is not y['S']
   assert c['off']==0 and self.owners(c['S'])==1
   den=c['q'].denominator*dq.denominator//gcd(c['q'].denominator,dq.denominator);aa=int(c['q']*den);bb=int(dq*den);g=gcd(abs(aa),abs(bb));wa=aa//g if g else 0;wb=bb//g if g else 0
   q=Q(g,den)if g else Q(0);M=abs(wa)*c['M']+abs(wb)*dM;self.fields(q,M)
   S=c['S'];oldbox=self.nodes[C]['exponent_box5'];lo=[b[0]for b in oldbox];hi=[b[1]for b in oldbox];counts=dict(self.nodes[C]['H_slice_sizes'])
   if not S:lo=[191]*5;hi=[-64]*5
   self.live.pop(C);self.live[name]={'S':S,'q':q,'M':M,'off':0}
  d=x['off']+y['off']-BIAS
  for e in x['S']:
   for f in y['S']:
    v=e+f+d
    if v not in S:
     self.own(1+2048);S.add(v);coords=unpack(v)
     for j,a in enumerate(coords):lo[j]=min(lo[j],a);hi[j]=max(hi[j],a)
     counts[coords[0]]=counts.get(coords[0],0)+1
  self.pair_count+=len(x['S'])*len(y['S'])
  if C is not None:
   dot={'terms_upper':len(x['S'])*len(y['S']),'known_QQ_content':str(dq),'primitive_integer_abs_upper':str(dM),'support_sha256':None,'exponent_box5':db,'H_slice_sizes':[[h,n]for h,n in sorted(hupper.items())],'support_profile':'UNMATERIALIZED dot: cardinality and Hcounts are SAFE UPPERS, not exact'}
   self.nodes[dotname]=dot;op={'output':dotname,'op':'unmaterialized_multiply_for_SOURCE_union','A':A,'B':B,'qM_and_BOX_bound':True}
   self.ops.append(op);self.meta+=count_numeric(dot)+count_numeric(op)+12
  # _record must not regard the moved, already live root as a second set.
  return self._record(name,S,q,M,0,{'op':'SOURCE_fused_union_product','C':C,'A':A,'B':B,'literal_dot':dotname,'integer_weight_old':str(wa),'integer_weight_dot':str(wb),'no_dot_support_allocated':True},[[lo[j],hi[j]]for j in range(5)]if S else[[0,0]]*5,sorted(counts.items()))
 def boundary_hash(self,name,expected=None):
  assert self.live[name]['off']==0
  return super().boundary_hash(name,expected)
 def snapshot(self,path):
  assert all(v['off']==0 for v in self.live.values())
  assert len({id(v['S'])for v in self.live.values()})==len(self.live)
  return super().snapshot(path)
 def failure_view_snapshot(self,path):
  """Lossless diagnostic state, explicitly NOT a completed pending operation."""
  from pathlib import Path
  import hashlib,json
  path=Path(path);assert not path.exists()
  # The always-reserved4096 SOURCE scalar slots pay for this failure-only
  # header/index (<=2048), even while caller scratch remains live. Do not
  # recurse through own/checkpoint while consuming that existing reserve.
  names=sorted(self.live);assert len(names)<=64
  blocks=[];index={};maps=[]
  for k in names:
   v=self.live[k];sid=id(v['S'])
   if sid not in index:index[sid]=len(blocks);blocks.append(v['S'])
   self.fields(v['q'],v['M']);assert len(k)<=80
   maps.append({'id':k,'b':index[sid],'o':v['off'],'q':str(v['q']),'M':str(v['M']),'n':len(v['S'])})
  total=sum(map(len,blocks));assert total<=125000
  header={'format':'Q5VIEW1','pending_operation_not_complete':True,'context':self.context,'maps':maps,'blocks':[len(S)for S in blocks]}
  assert count_numeric(header)+6*len(blocks)+4*len(names)<=2048
  b=json.dumps(header,separators=(',',':')).encode();assert len(b)<=1024*64+4096
  tmp=path.with_suffix(path.suffix+'.partial');assert not tmp.exists()
  with tmp.open('wb')as f:
   f.write(b'Q5VIEW1\n');f.write(len(b).to_bytes(4,'big'));f.write(b)
   for S in blocks:
    for e in S:assert 0<=e<2**40;f.write(e.to_bytes(5,'big'))
  with tmp.open('rb')as f:
   assert f.read(8)==b'Q5VIEW1\n';n=int.from_bytes(f.read(4),'big');assert n==len(b)and f.read(n)==b
   for S in blocks:
    for e in S:assert f.read(5)==e.to_bytes(5,'big')
   assert not f.read(1)
  tmp.rename(path);h=hashlib.sha256()
  with path.open('rb')as f:
   for b in iter(lambda:f.read(65536),b''):h.update(b)
  return{'file':path.name,'sha256':h.hexdigest(),'bytes':path.stat().st_size,'maps':len(names),'backings':len(blocks),'unique_keys':total,'all_backing_keys_header_and_view_offsets_stream_roundtrip':True,'pending_operation_not_completed':True,'not_VALUE_data':True}

#!/usr/bin/env python3
"""ONE source-faithful tail3 and COMBINED D12*tail3 in the same Haxis quotient."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1];P=32003;OWN=12000;RAW=8388608;PRODUCTS=100000
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(custody,custody_sha,out):
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60)
 start=time.monotonic();out.mkdir(exist_ok=False,parents=True);pins={};readbytes=0;peak=0;bits=0
 s=dict(version='J3HAX001',stage='pins',p=P,products=0,pending=0,cursor={},work={},nodes=[],powers={})
 def pin(p,h):
  nonlocal readbytes
  assert sha(p)==h;readbytes+=p.stat().st_size;pins[str(p.relative_to(ROOT))]=h
 def guard():
  nonlocal peak
  n=6000+4*len(s['nodes'])+5*len(s['powers']);peak=max(peak,n);assert n<=OWN
  assert time.monotonic()-start<45,'soft_clock'
 def mul(a,b):
  nonlocal bits
  s['pending']=a*b;bits=max(bits,s['pending'].bit_length());s['products']+=1
  assert s['products']<=PRODUCTS and s['pending'].bit_length()<=30
  if s['products']%512==0:guard()
  v=s['pending']%P;s['pending']=0;return v
 def inv(a):
  assert a%P;r=1;b=a%P;n=P-2
  while n:
   s['inverse_cursor']=dict(a=a,r=r,b=b,n=n)
   if n&1:r=mul(r,b)
   b=mul(b,b);n>>=1
  assert mul(a%P,r)==1;s['inverse_cursor']={};return r
 def rm(a,b):
  c=[0]*7;s['work']=dict(left=a,right=b,current=c,phase='convolution')
  for i in range(4):
   for j in range(4):s['work'].update(i=i,j=j);c[i+j]=(c[i+j]+mul(a[i],b[j]))%P
  for k in range(6,3,-1):
   t=c[k];s['work'].update(phase='monic_reduction',k=k,t=t)
   for j in range(4):s['work']['j']=j;c[k-4+j]=(c[k-4+j]-mul(t,s['h'][j]))%P
   c[k]=0
  z=c[:4];s['work']={};return z
 def dump(n,obj):
  with (out/n).open('w') as f:
   for chunk in json.JSONEncoder(separators=(',',':')).iterencode(obj):
    assert readbytes+sum(p.stat().st_size for p in out.iterdir() if p.is_file())+len(chunk.encode())+65536<=RAW
    f.write(chunk)
   f.write('\n')
 try:
  pin(custody,custody_sha);c=json.loads(custody.read_text());assert c['prepared_only']
  pin(Path(__file__),c['worker_sha256'])
  for p,h in c['runtime_input_sha256'].items():pin(ROOT/p,h)
  s['h']=c['monic_H_axis'];assert s['h']==[17430,11013,30453,27809,1]
  s['powers'][0]=[1,0,0,0];s['powers'][1]=[0,1,0,0]
  v=inv(s['h'][0]);ai=[mul((-s['h'][k])%P,v) for k in (1,2,3)]+[(-v)%P]
  assert rm(s['powers'][1],ai)==[1,0,0,0];s['powers'][-1]=ai
  for k in range(2,c['maximum_a_exponent']+1):s['powers'][k]=rm(s['powers'][k-1],s['powers'][1])
  for k in range(-2,c['minimum_a_exponent']-1,-1):s['powers'][k]=rm(s['powers'][k+1],ai)
  s['stage']='stream_literal_atoms';s['nodes']=[[0]*4 for _ in c['instructions']]
  with (ROOT/c['native_input']).open() as f:
   assert f.readline().strip()=='J3HAX001';assert int(f.readline())==c['literal_atom_count']
   for atom in range(c['literal_atom_count']):
    node,D,n=map(int,f.readline().split());assert D%P and n>=0
    z=[0]*4;s['current_atom']=z;s['cursor']=dict(atom=atom,node=node,D=D,terms=n,term=0)
    for t in range(n):
     a,e,d,v=map(int,f.readline().split());s['cursor']['term']=t
     assert c['minimum_a_exponent']<=a<=c['maximum_a_exponent'] and e>=0 and d>=0
     if e==d==0:
      for j in range(4):z[j]=(z[j]+mul(v%P,s['powers'][a][j]))%P
    dv=inv(D%P);s['nodes'][node]=[mul(x,dv) for x in z]
   assert not f.read().strip()
  s['stage']='complete_fixed_typed_DAG'
  # Column closure plus cached4; four dot products and their sum are appended
  # by the same source plan and are active even though outside that prefix.
  active=set(c['active_column_and_cached_nodes'])|set(range(155,len(c['instructions'])))
  for i,(op,a,b) in enumerate(c['instructions']):
   if op==0 or i not in active:continue
   assert 0<=a<i and 0<=b<i;s['cursor']=dict(node=i,op=op,a=a,b=b)
   if op==3:s['nodes'][i]=rm(s['nodes'][a],s['nodes'][b])
   else:s['nodes'][i]=[(x+y if op==1 else x-y)%P for x,y in zip(s['nodes'][a],s['nodes'][b])]
  s['T3_remainder']=s['nodes'][c['output_root']]
  # A zero-divisor ring requires this single combined test. Separate nonzero
  # factors would not establish a common component/open intersection.
  r=json.loads((ROOT/c['D12_report']).read_text());assert r['all_passed']
  s['D12_remainder']=r['D12_remainder'];s['combined_remainder']=rm(s['D12_remainder'],s['T3_remainder'])
  s['stage']='completed';guard();positive=any(s['combined_remainder']);dump('complete-fixed-source-ring-certificate.json',s)
  r=dict(id=out.name,all_passed=True,computations_completed=True,combined_D12_times_T3_NONZERO_PRIMARY=positive,
   separate_nonzero_factors_NOT_used=True,zero_combined_remainder_inconclusive=not positive,
   original32_target=3,current15_pair_rows=[2,4,5,8],source_q2_roots=[2551,2761,2866,3181],
   genuine_global_minor_identity='A_anchor*D15_3=-T3*D12, only six named FULLU units',
   no_kappa_or_oldOmega_inverse=True,T3_remainder=s['T3_remainder'],combined_remainder=s['combined_remainder'],
   products=s['products'],peak_owned=peak,max_field_temporary_bits=bits,input_sha256=pins,
   evidence_sha256={'complete-fixed-source-ring-certificate.json':sha(out/'complete-fixed-source-ring-certificate.json')},
   independent_value_peer_required=True,no_full_H_coverage_Fitting_or_bad_rank_claim=True,seconds=time.monotonic()-start)
  dump('report.json',r);print(json.dumps(r))
 except BaseException as e:
  s['failure']=repr(e);dump('failure-complete-prefix.json',dict(state=s,input_sha256=pins,no_retry=True));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--custody',type=Path,required=True);p.add_argument('--custody-sha',required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();main(a.custody.resolve(),a.custody_sha,a.output.resolve())

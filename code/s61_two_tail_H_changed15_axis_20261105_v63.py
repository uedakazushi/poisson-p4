#!/usr/bin/env python3
"""ONE changed genuine15 minor via pair4 and combined D12 over fixed Haxis ring."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1];P=32003;CAP=10000;OWN=4096;RAW=4194304
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(custody,custody_sha,out):
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60)
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={};readbytes=0;bits=0;peak=0
 s=dict(version='HCHG1501',stage='pins',p=P,products=0,pending=0,cursor={},work={},matrix=[],powers={},dp=[0]*64,current=[0]*4)
 def pin(p,h):
  nonlocal readbytes
  assert sha(p)==h;readbytes+=p.stat().st_size;pins[str(p.relative_to(ROOT))]=h
 def guard():
  nonlocal peak
  n=2048+64+4*len(s['matrix'])+5*len(s['powers']);peak=max(peak,n);assert n<=OWN
  assert time.monotonic()-start<45,'soft_clock'
 def mul(a,b):
  nonlocal bits
  s['pending']=a*b;s['products']+=1;bits=max(bits,s['pending'].bit_length())
  assert s['products']<=CAP and s['pending'].bit_length()<=30
  if s['products']%128==0:guard()
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
  v=c[:4];s['work']={};return v
 def dump(name,obj):
  with (out/name).open('w') as f:
   for chunk in json.JSONEncoder(separators=(',',':')).iterencode(obj):
    assert readbytes+sum(p.stat().st_size for p in out.iterdir() if p.is_file())+len(chunk.encode())+65536<=RAW
    f.write(chunk)
   f.write('\n')
 try:
  pin(custody,custody_sha);c=json.loads(custody.read_text());assert c['prepared_only']
  pin(Path(__file__),c['worker_sha256'])
  for p,h in c['runtime_input_sha256'].items():pin(ROOT/p,h)
  s['h']=[17430,11013,30453,27809,1];s['powers'][0]=[1,0,0,0];s['powers'][1]=[0,1,0,0]
  v=inv(s['h'][0]);ai=[mul((-s['h'][k])%P,v) for k in (1,2,3)]+[(-v)%P]
  assert rm(s['powers'][1],ai)==[1,0,0,0];s['powers'][-1]=ai
  for k in range(2,c['maximum_a_exponent']+1):s['powers'][k]=rm(s['powers'][k-1],s['powers'][1])
  for k in range(-2,c['minimum_a_exponent']-1,-1):s['powers'][k]=rm(s['powers'][k+1],ai)
  s['stage']='sixteen_saved_axis_polynomials'
  for index,q in enumerate(c['axis_input_polynomials']):
   D=int(q['D']);assert D%P;z=[0]*4;s['current']=z
   for t,(a,v) in enumerate(q['coefficients']):
    s['cursor']=dict(entry=index,term=t,a=a,integer=v,component=0)
    for j in range(4):s['cursor']['component']=j;z[j]=(z[j]+mul(int(v)%P,s['powers'][a][j]))%P
   dv=inv(D%P);s['matrix'].append([mul(x,dv) for x in z])
  s['stage']='fixed_pair4_determinant';s['dp'][0]=1
  for mask in range(1,16):
   row=mask.bit_count()-1;v=[0]*4;s['current']=v
   for j in range(4):
    if not mask>>j&1:continue
    old=mask^(1<<j);negative=(row+(mask&((1<<j)-1)).bit_count())&1;s['cursor']=dict(mask=mask,row=row,column=j,previous=old,negative=negative)
    z=rm(s['dp'][4*old:4*old+4],s['matrix'][4*row+j])
    for k in range(4):v[k]=(v[k]-z[k] if negative else v[k]+z[k])%P
   s['dp'][4*mask:4*mask+4]=v
  s['detPair_remainder']=s['dp'][-4:];s['stage']='COMBINED_D12_times_pair'
  peer=json.loads((ROOT/c['D12_peer']).read_text());assert peer['all_passed'] and peer['H_nondivides_D12_in_BASE6_independently_proved']
  s['D12_remainder']=peer['D12_remainder'];s['combined_remainder']=rm(s['D12_remainder'],s['detPair_remainder'])
  s['stage']='completed';guard();positive=any(s['combined_remainder']);dump('fixed-pair4-and-combined-certificate.json',s)
  r=dict(id=out.name,all_passed=True,computations_completed=True,combined_D12_times_changed_pair_NONZERO_PRIMARY=positive,
   zero_combined_inconclusive=not positive,no_field_or_etale_assumption=True,separate_nonzero_factors_NOT_used=True,
   source_identity='A*D15_new=+detPair4*D12',no_anchor_kappa_oldOmega_inverse=True,
   detPair_remainder=s['detPair_remainder'],combined_remainder=s['combined_remainder'],products=s['products'],peak_owned=peak,max_field_temporary_bits=bits,
   input_sha256=pins,evidence_sha256={'fixed-pair4-and-combined-certificate.json':sha(out/'fixed-pair4-and-combined-certificate.json')},
   independent_peer_required=True,no_whole_H_coverage_or_Fitting_or_bad_rank_claim=True,seconds=time.monotonic()-start)
  dump('report.json',r);print(json.dumps(r))
 except BaseException as e:
  s['failure']=repr(e);dump('failure-complete-prefix.json',dict(state=s,input_sha256=pins,no_retry=True));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--custody',type=Path,required=True);p.add_argument('--custody-sha',required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();main(a.custody.resolve(),a.custody_sha,a.output.resolve())

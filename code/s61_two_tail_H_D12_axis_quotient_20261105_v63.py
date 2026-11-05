#!/usr/bin/env python3
"""ONE fixed D12-axis determinant modulo H_axis over F32003; help only until GO."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1];P=32003;OWN=28000;RAW=1048576;CAP=750000
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(custody,custody_sha,out):
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60)
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={};readbytes=0
 s=dict(version='HAXD1201',stage='pins',p=P,products=0,pending=0,cursor={},work={},matrix=[],powers=[],dp=[],current=[0]*4)
 peak=0;bits=0
 def pin(p,h):
  nonlocal readbytes
  assert sha(p)==h;readbytes+=p.stat().st_size;pins[str(p.relative_to(ROOT))]=h
 def guard():
  nonlocal peak
  # Charge full allocated DP plus all fixed inputs/control, regardless of zeros.
  owned=10000+len(s['dp'])+4*len(s['matrix'])+4*len(s['powers'])
  peak=max(peak,owned);assert owned<=OWN
  assert time.monotonic()-start<45,'soft_clock'
 def mul(a,b):
  nonlocal bits
  s['pending']=a*b;s['products']+=1;bits=max(bits,abs(s['pending']).bit_length())
  assert s['products']<=CAP and abs(s['pending']).bit_length()<=30
  if s['products']%2048==0:guard()
  z=s['pending']%P;s['pending']=0;return z
 def inv(a):
  assert a%P;r=1;b=a%P;n=P-2
  while n:
   s['inverse_cursor']=dict(input=a,result=r,base=b,exponent=n)
   if n&1:r=mul(r,b)
   b=mul(b,b);n>>=1
  assert mul(a%P,r)==1;s['inverse_cursor']={};return r
 def rm(a,b):
  assert len(a)==len(b)==4;c=[0]*7
  s['work']=dict(left=a,right=b,current=c,phase='convolution',i=0,j=0)
  for i in range(4):
   for j in range(4):
    s['work'].update(i=i,j=j);c[i+j]=(c[i+j]+mul(a[i],b[j]))%P
  for k in range(6,3,-1):
   t=c[k];s['work'].update(phase='monic_reduction',k=k,t=t)
   for j in range(4):
    s['work']['j']=j;c[k-4+j]=(c[k-4+j]-mul(t,s['h'][j]))%P
   c[k]=0
  z=c[:4];s['work']={};return z
 def dump(name,obj):
  p=out/name
  with p.open('w') as f:
   for chunk in json.JSONEncoder(separators=(',',':')).iterencode(obj):
    used=readbytes+sum(x.stat().st_size for x in out.iterdir() if x.is_file())+len(chunk.encode())+65536
    assert used<=RAW,('RAW',used);f.write(chunk)
   f.write('\n')
 try:
  pin(custody,custody_sha);c=json.loads(custody.read_text());assert c['preparation_only']
  pin(Path(__file__),c['worker_sha256'])
  for p,h in c['runtime_input_sha256'].items():pin(ROOT/p,h)
  mpath=ROOT/c['matrix_input'];m=json.loads(mpath.read_text());M=m['matrix'];assert len(M)==12 and all(len(r)==12 for r in M)
  s['stage']='fixed_axis_modulus';s['H_axis_integer']=c['H_axis_integer'];h=[x%P for x in c['H_axis_integer']]
  assert len(h)==5 and h[-1] and h[0];v=inv(h[-1]);s['h']=[mul(x,v) for x in h];assert s['h'][-1]==1
  lcd=int(m['LCD_product']);assert lcd%P;vinv=inv(lcd%P);s['LCD_mod_p']=lcd%P;s['LCD_inverse']=vinv
  s['powers']=[[1,0,0,0]]
  for _ in range(8):s['powers'].append(rm(s['powers'][-1],[0,1,0,0]))
  s['stage']='input_reduction'
  for i in range(12):
   for j in range(12):
    v=[0]*4;s['current']=v
    for k,x in M[i][j]:
     s['cursor']=dict(row=i,column=j,source_degree=k);x=int(x)%P;assert 0<=k<=8
     for t in range(4):v[t]=(v[t]+mul(x,s['powers'][k][t]))%P
    s['matrix'].append(v)
  del M,m;guard();s['stage']='subset_determinant';s['dp']=[0]*(4096*4);s['dp'][0]=1
  for mask in range(1,4096):
   row=mask.bit_count()-1;v=[0]*4;s['current']=v
   for j in range(12):
    if not mask>>j&1:continue
    prev=mask^(1<<j);negative=(row+(mask&((1<<j)-1)).bit_count())&1
    s['cursor']=dict(mask=mask,row=row,column=j,previous=prev,negative=negative)
    z=rm(s['dp'][4*prev:4*prev+4],s['matrix'][12*row+j])
    for k in range(4):v[k]=(v[k]-z[k] if negative else v[k]+z[k])%P
   s['dp'][4*mask:4*mask+4]=v
   if mask%32==0:guard()
  s['detC_remainder']=s['dp'][-4:];s['stage']='a5_and_LCD_normalization'
  rem=rm(s['detC_remainder'],s['powers'][5]);s['D12_remainder']=[mul(x,vinv) for x in rem]
  # Explicit a-unit identity in the quotient, with no field assumption on S.
  hi=inv(s['h'][0]);ainv=[mul((-s['h'][k])%P,hi) for k in (1,2,3)]+[(-hi)%P]
  assert rm([0,1,0,0],ainv)==[1,0,0,0];s['a_inverse']=ainv
  s['stage']='completed';guard();positive=any(s['D12_remainder'])
  dump('complete-subset-and-ring-certificate.json',s)
  r=dict(id=out.name,all_passed=True,computations_completed=True,H_nondivides_D12_PRIMARY=positive,
   zero_remainder_is_inconclusive=not positive,no_field_or_etale_assumption=True,existing_p=P,
   no_new_beta_point_or_source_or_other_axis=True,products=s['products'],peak_owned=peak,max_field_temporary_bits=bits,
   detC_remainder=s['detC_remainder'],D12_remainder=s['D12_remainder'],input_sha256=pins,
   evidence_sha256={'complete-subset-and-ring-certificate.json':sha(out/'complete-subset-and-ring-certificate.json')},
   independent_peer_required_before_nondivision_adoption=True,no_actual_Fitting_or_fullrank_claim=True,seconds=time.monotonic()-start)
  dump('report.json',r);print(json.dumps(r))
 except BaseException as e:
  s['failure']=repr(e);dump('failure-complete-prefix.json',dict(state=s,input_sha256=pins,no_retry=True));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--custody',type=Path,required=True);p.add_argument('--custody-sha',required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();main(a.custody.resolve(),a.custody_sha,a.output.resolve())

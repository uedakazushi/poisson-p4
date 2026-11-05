#!/usr/bin/env python3
"""ONE fixed H(e=0,d=1) Rabin test over F32003; --help performs no values."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
P=32003;MODCAP=14000;OWNCAP=8192;RAWCAP=1048576
FACTOR='reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63/factor-sparse.txt'
FACTOR_SHA='241488fafa9b45e5353f0b707f75902bc5bc32e85ee093b9f2e0a5f07c7d571c'
PEER='reports/S61_TWO_TAIL_CRT_ROOT_DENSE_REVIEW_20261105_v63/report.json'
PEER_SHA='6d7f852148000830f26d4f9ed709a7e7a5bc205e8c2678e4613648ea82b23bc7'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def numeric(x):
 if isinstance(x,bool) or x is None:return 0
 if isinstance(x,int):return 1
 if isinstance(x,(list,tuple)):return sum(numeric(v) for v in x)
 if isinstance(x,dict):return len(x)+sum(numeric(v) for v in x.values())
 return 0
def main(custody,custody_sha,out):
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60)
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 s=dict(version='H6RAB001',stage='pins',p=P,modular_products=0,peak_owned=0,pending=0,
        factor=[],specialized=[],monic=[],frobenius=[[0,1]],euclid=[],cursor={},work={},scalar_inverse={},gcd_work={})
 pins={};maxbits=0;rawread=0
 def pin(p,expected):
  nonlocal rawread
  h=sha(p);assert h==expected,(str(p),h);pins[str(p.relative_to(ROOT))]=h;rawread+=p.stat().st_size
 def guard():
  owned=2048+numeric(s);s['peak_owned']=max(s['peak_owned'],owned);assert owned<=OWNCAP,('owned',owned)
  assert time.monotonic()-started<45,'soft_clock'
 def dump(n,x):
  b=(json.dumps(x,separators=(',',':'))+'\n').encode()
  # Include all external reads, all already-created files and this complete
  # new payload; no atomic second copy or omitted failure-prefix reserve.
  used=rawread+sum(p.stat().st_size for p in out.iterdir() if p.is_file())+len(b)
  assert used<=RAWCAP,('RAW',used);(out/n).write_bytes(b)
 def mul(a,b):
  nonlocal maxbits
  s['pending']=a*b;s['modular_products']+=1
  maxbits=max(maxbits,abs(s['pending']).bit_length())
  assert s['modular_products']<=MODCAP and abs(s['pending']).bit_length()<=30
  if s['modular_products']%128==0:guard()
  v=s['pending']%P;s['pending']=0;return v
 def inv(a):
  assert a%P
  r=1;b=a%P;n=P-2
  while n:
   s['scalar_inverse']=dict(input=a,result=r,base=b,pending_exponent=n)
   if n&1:r=mul(r,b)
   s['scalar_inverse']['result']=r;b=mul(b,b);n>>=1
  assert mul(a%P,r)==1
  s['scalar_inverse']={}
  return r
 def trim(a):
  while a and a[-1]==0:a.pop()
  return a
 def pmul(a,b,f):
  assert len(a)<=6 and len(b)<=6 and len(f)==7 and f[-1]==1
  c=[0]*11
  s['work']=dict(kind='pmul',left=a,right=b,modulus=f,current=c,cursor={})
  for i,x in enumerate(a):
   for j,y in enumerate(b):
    s['work']['cursor']=dict(phase='convolution',i=i,j=j);c[i+j]=(c[i+j]+mul(x,y))%P
  for k in range(10,5,-1):
   t=c[k]
   for j in range(6):
    s['work']['cursor']=dict(phase='monic_reduction',k=k,j=j);c[k-6+j]=(c[k-6+j]-mul(t,f[j]))%P
   c[k]=0
  result=trim(c[:6]);s['work']={};return result
 def ppow(a,f,layer):
  r=[1];b=a[:];n=P
  while n:
   s['cursor']=dict(layer=layer,pending_exponent=n,result=r,base=b);guard()
   if n&1:r=pmul(r,b,f)
   s['cursor']['result']=r;b=pmul(b,b,f);n>>=1
  s['cursor']={};return r
 def pdiv(a,b):
  assert b and len(a)<=7 and len(b)<=7
  r=a[:];q=[0]*max(0,len(a)-len(b)+1);steps=0
  s['work']=dict(kind='pdiv',left=a,right=b,current_remainder=r,current_quotient=q,cursor={})
  v=inv(b[-1])
  while r and len(r)>=len(b):
   k=len(r)-len(b);t=mul(r[-1],v);q[k]=t
   for j,c in enumerate(b):
    s['work']['cursor']=dict(shift=k,j=j,t=t);r[k+j]=(r[k+j]-mul(t,c))%P
   trim(r);steps+=1;assert steps<=7
  result=trim(q);s['work']={};return result,r
 def pgcd(a,b,which):
  chain=[];s['euclid'].append(dict(which=which,steps=chain));steps=0
  while b:
   s['gcd_work']=dict(which=which,left=a,right=b)
   guard();q,r=pdiv(a,b);chain.append(dict(left=a,right=b,quotient=q,remainder=r))
   a,b=b,r;steps+=1;assert steps<=7
  v=inv(a[-1]);result=[mul(c,v) for c in a];s['gcd_work']={};return result
 try:
  pin(custody,custody_sha);c=json.loads(custody.read_text())
  assert c['preparation_only'] and c['no_target_values_or_Rabin_executed']
  pin(Path(__file__),c['worker_sha256']);pin(ROOT/FACTOR,FACTOR_SHA);pin(ROOT/PEER,PEER_SHA)
  peer=json.loads((ROOT/PEER).read_text());assert peer['all_passed'] and peer['exact_char0_gcd_proved_up_to_QQ_scalar']
  with (ROOT/FACTOR).open() as f:
   assert int(f.readline())==94
   s['factor']=[tuple(map(int,line.split())) for line in f]
  assert len(s['factor'])==94 and all(abs(v).bit_length()<=40 for a,e,d,v in s['factor'])
  s['stage']='fixed_e0_d1_specialization';f=[0]*7;s['specialized']=f
  for a,e,d,v in s['factor']:
   assert 0<=a<=6 and 0<=e<=9 and 0<=d<=4
   if e==0:f[a]=(f[a]+v)%P
  s['specialized']=f[:];assert f[6]!=0,'fixed_specialization_degree_dropped'
  v=inv(f[6]);f=[mul(x,v) for x in f];assert f[-1]==1;s['monic']=f[:]
  s['stage']='six_fixed_frobenius_powers'
  for layer in range(1,7):s['frobenius'].append(ppow(s['frobenius'][-1],f,layer));guard()
  s['stage']='two_fixed_Rabin_gcds'
  gcds={}
  for k in (3,2):
   b=s['frobenius'][k][:]
   while len(b)<2:b.append(0)
   b[1]=(b[1]-1)%P;gcds[str(k)]=pgcd(f[:],trim(b),k)
  s['stage']='completed'
  irreducible=s['frobenius'][6]==[0,1] and all(v==[1] for v in gcds.values());guard()
  evidence=dict(state=s,Rabin_gcds=gcds,x_p6_equals_x=s['frobenius'][6]==[0,1],
   fixed_specialization_degree=6,one_fixed_specialization_only=True,max_field_temporary_bits=maxbits)
  dump('fixed6-Rabin-certificate.json',evidence)
  rep=dict(id=out.name,all_passed=irreducible,computations_completed=True,irreducible_over_QQ_a_e_d_verified=irreducible,
   fixed_p=P,fixed_e=0,fixed_d=1,no_reset_extra_prime_point_or_factor_family=True,
   modular_products=s['modular_products'],peak_owned=s['peak_owned'],max_field_temporary_bits=maxbits,
   input_sha256=pins,evidence_sha256={'fixed6-Rabin-certificate.json':sha(out/'fixed6-Rabin-certificate.json')},
   no_Fitting_or_fullrank_or_exhaustion_claim=True,seconds=time.monotonic()-started)
  dump('report.json',rep);print(json.dumps(rep))
 except BaseException as e:
  s['guard_failure']=repr(e)
  # Micro-state is bounded by the same source reserve; no computation retry.
  dump('failure-full-state.json',dict(state=s,input_sha256=pins,max_field_temporary_bits=maxbits,no_retry=True))
  raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--custody',type=Path,required=True)
 p.add_argument('--custody-sha',required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();main(a.custody.resolve(),a.custody_sha,a.output.resolve())

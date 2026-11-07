#!/usr/bin/env python3
"""ONE direct3var factor call plus conditional384bit exact multiplyback; no automatic retry."""
from pathlib import Path
import argparse,hashlib,json,math,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1];RAW=4194304;OWN=50000;PAIRS=1500000
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def numeric(x):
 if isinstance(x,bool) or x is None:return 0
 if isinstance(x,(int,float)):return 1
 if isinstance(x,str):return int(x.lstrip('-').isdigit())
 if isinstance(x,(tuple,list)):return sum(numeric(v) for v in x)
 if isinstance(x,dict):return len(x)+sum(numeric(k)+numeric(v) for k,v in x.items())
 return 0
def main(custody,custody_sha,out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 out.mkdir(exist_ok=False,parents=True);pins={};readbytes=0;peak=0;maxbits=0
 s=dict(version='HQQ38401',phase='pins',factors=[],original={},current={0:1},next={},pending_product=0,pending_sum=0,cursor={},pairs=0,stages=[])
 def pin(p,h):
  nonlocal readbytes
  assert sha(p)==h;readbytes+=p.stat().st_size;pins[str(p.relative_to(ROOT))]=h
 def dump(n,x):
  other=readbytes+sum(p.stat().st_size for p in out.iterdir() if p.is_file());written=0
  with (out/n).open('w') as f:
   for chunk in json.JSONEncoder(separators=(',',':')).iterencode(x):
    written+=len(chunk.encode());assert other+written+65536<=RAW,('RAW',other+written);f.write(chunk)
   f.write('\n')
 def guard(extra=0):
  nonlocal peak
  n=12000+numeric(s)+extra;peak=max(peak,n);assert n<=OWN,('owned',n)
  assert time.monotonic()-start<45,'soft_clock'
 def finish_frontier(reason):
  s['frontier_reason']=reason;dump('verification-complete-prefix.json',s)
  r=dict(id=out.name,all_passed=False,frontier_only=True,no_verified_factorization_claim=True,reason=reason,input_sha256=pins,
   library_internal_growth_opaque=True,no_retry_or_capraise=True,evidence_sha256={'verification-complete-prefix.json':sha(out/'verification-complete-prefix.json')})
  dump('report.json',r);print(json.dumps(r))
 def shape(poly):
  keys=list(poly);return [max(k//50 for k in keys),max((k%50)//5 for k in keys),max(k%5 for k in keys)]
 try:
  pin(custody,custody_sha);c=json.loads(custody.read_text());assert c['prepared_only']
  for p,h in c['runtime_input_sha256'].items():pin(ROOT/p,h)
  pin(Path(__file__),c['driver_sha256']);pin(ROOT/c['native_executable'],c['native_executable_sha256'])
  with (ROOT/c['fixed_H_input']).open() as f:
   assert int(f.readline())==94
   for line in f:
    a,e,d,v=map(int,line.split());assert abs(v).bit_length()<=40 and 0<=a<=6 and 0<=e<=9 and 0<=d<=4
    k=50*a+5*e+d;assert k not in s['original'];s['original'][k]=v
  assert len(s['original'])==94 and shape(s['original'])==[6,9,4]
  s['phase']='ONE opaque native3var factor call';dump('phase-before-native.json',dict(input_sha256=pins,opaque_API_internal_bits_not_384_certified=True,no_factor_values_yet=True))
  sec=45-(time.monotonic()-start);assert sec>1
  cmd=[str(ROOT/c['native_executable']),'--run',str(ROOT/c['fixed_H_input']),str(out),str(readbytes),str(sec)]
  try:
   with (out/'native.log').open('w') as log:p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,58-(time.monotonic()-start)))
  except subprocess.TimeoutExpired:finish_frontier('opaque native library timeout; inputs and logs preserved, internal library state unavailable');return
  if p.returncode!=0:finish_frontier('native API/canonicalization stopped; no verified factorization');return
  status=json.loads((out/'native-status.json').read_text());assert status['candidate_list_completed'] and status['factor_API_status']==1
  doc=json.loads((out/'canonical-factor-candidates.json').read_text());assert doc['factor_API_status']==1 and 1<=len(doc['factors'])<=11
  s['phase']='canonical candidate parsing'
  for f in doc['factors']:
   poly={};g=0
   for (a,e,d),v in f['terms']:
    v=int(v);assert v and abs(v).bit_length()<=357 and 0<=a<=6 and 0<=e<=9 and 0<=d<=4
    k=50*a+5*e+d;assert k not in poly;poly[k]=v;g=math.gcd(g,abs(v))
   assert len(poly)<=350 and g==1 and poly[max(poly)]>0
   s['factors'].append(dict(multiplicity=f['multiplicity'],polynomial=poly,maxbits=max(abs(v).bit_length() for v in poly.values()),degrees=shape(poly)))
  # Count even the transient original JSON plus converted integer candidates.
  guard(numeric(doc));del doc
  assert [sum(f['multiplicity']*f['degrees'][j] for f in s['factors']) for j in range(3)]==[6,9,4]
  s['phase']='conditional stagewise384bit exact multiplyback'
  for fi,f in enumerate(s['factors']):
   for repeat in range(f['multiplicity']):
    A=s['current'];B=f['polynomial'];ab=max(abs(v).bit_length() for v in A.values());bb=f['maxbits'];upper=ab+bb+18
    if upper>384:finish_frontier('conditional exact-product source bound exceeds384 before any product of this stage');return
    degA=shape(A);assert all(a+b<=u for a,b,u in zip(degA,f['degrees'],[6,9,4]))
    s['next']={};s['stages'].append(dict(factor=fi,repeat=repeat,partial_bits=ab,factor_bits=bb,prospective_bound=upper))
    for ka,x in A.items():
     for kb,y in B.items():
      k=ka+kb;assert 0<=k<350;s['cursor']=dict(factor=fi,repeat=repeat,left=ka,right=kb,result=k,step='before_product')
      s['pairs']+=1;assert s['pairs']<=PAIRS
      # BEFORE this multiply the bound is proven from actual operand maxima
      # and all<=350² possible contributions; no cancellation is assumed.
      s['pending_product']=x*y;assert abs(s['pending_product']).bit_length()<=384
      s['cursor']['step']='before_add';s['pending_sum']=s['next'].get(k,0)+s['pending_product'];assert abs(s['pending_sum']).bit_length()<=upper<=384
      maxbits=max(maxbits,abs(s['pending_product']).bit_length(),abs(s['pending_sum']).bit_length())
      if k not in s['next']:guard(3)
      if s['pending_sum']:s['next'][k]=s['pending_sum']
      else:s['next'].pop(k,None)
      s['pending_product']=s['pending_sum']=0;s['cursor']['step']='committed'
      if s['pairs']%4096==0:guard()
    s['current']=s['next'];s['next']={};assert len(s['current'])<=350;guard()
  assert s['current']==s['original'],'ALL exact multivariate coefficient product mismatch'
  assert all(f['multiplicity']==1 for f in s['factors']),'conflicts with accepted squarefree input'
  s['phase']='complete exact factorization PRIMARY';dump('canonical-factors-and-full-product-certificate.json',s)
  r=dict(id=out.name,all_passed=True,computations_completed=True,exact_QQ_factorization_PRIMARY=True,
   primitive_positive_factor_count=len(s['factors']),QQ_reducible_PRIMARY=len(s['factors'])>1,QQ_irreducible_API_contract_PRIMARY=len(s['factors'])==1,
   factor_irreducibility_uses_COMPLETE_FLINT_API_contract=True,geometric_irreducibility_NOT_claimed=True,
   all_original94_coefficients_and_other350box_zeros_exact=True,canonical_output_bound357=True,
   verified_identity_hard_bits=384,observed_verification_bits=maxbits,verification_pairs=s['pairs'],Python_peak_owned=peak,native_accessible_peak=status['accessible_peak'],
   native_API_internal_growth_OPAQUE_under_RAM_time=True,input_sha256=pins,
   evidence_sha256={n:sha(out/n) for n in ['canonical-factor-candidates.json','native-status.json','canonical-factors-and-full-product-certificate.json']},
   independent_factor_identity_and_irreducibility_peer_REQUIRED=True,no_bad_Fitting_fullrank_or_exhaustion_claim=True,seconds=time.monotonic()-start)
  dump('report.json',r);print(json.dumps(r))
 except BaseException as e:
  finish_frontier(type(e).__name__+': '+str(e));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--custody',type=Path,required=True);p.add_argument('--custody-sha',required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();main(a.custody.resolve(),a.custody_sha,a.output.resolve())

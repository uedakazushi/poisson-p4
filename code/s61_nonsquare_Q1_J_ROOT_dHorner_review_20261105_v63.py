#!/usr/bin/env python3
"""Distinct d-Horner replay of the first J layer with symbolic a^18 denominator."""
import hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PRIMARY=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_J_HASSE_ENTRY_20261105_v63'
SOURCE=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_J_DHORNER_SOURCE_20261105_v63'
ORIGINAL=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt'
OUT=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_J_ROOT_DHORNER_REVIEW_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));OUT.mkdir(exist_ok=False);pins={}
 for p,want in [(PRIMARY/'report.json','0e25e61b3dda79cfba5d6b51aedabd22b4961dd4509b705a87e78a0e3c59e2fa'),(SOURCE/'report.json','2d7e00ba86ec0aef44da55c17d65094a78e32b0206059bad6d58c2dfbd608ba0')]:
  assert sha(p)==want;d=json.loads(p.read_text());assert d['all_passed'];pins[str(p.relative_to(ROOT))]=want
  for name,h in d['evidence_sha256'].items():assert sha(p.parent/name)==h;pins[str((p.parent/name).relative_to(ROOT))]=h
 assert sha(ORIGINAL)=='d88cfe290ecbc71b045a382a333516815967ce6117576697c7cce78f2b00a6f7';pins[str(ORIGINAL.relative_to(ROOT))]=sha(ORIGINAL)
 groups=[{}for _ in range(19)];l1=0
 with ORIGINAL.open()as f:
  assert f.readline().strip()=='T1C4QN001';D=int(f.readline());assert int(f.readline())==4811
  for _ in range(4811):
   key,c=map(int,f.readline().split());a,rem=divmod(key,1048576);e,d=divmod(rem,1024)
   assert 0<=a<=16 and 0<=e<=33 and 0<=d<=18 and a+e+d<=50 and c
   groups[d][1024*a+e]=c;l1+=abs(c)
  assert not f.read().strip()
 assert l1<2**301 and D.bit_length()==197
 assert 2*(4811+2*1667)+16000==32290 and 2**301*10**18<2**361
 assert 18*2*1667+4811<65000
 (OUT/'pre-arithmetic-source-and-budget.json').write_text(json.dumps(dict(
  recurrence='H_l=(4+6e)H_(l+1)+(5a)^(18-l)I_l',source_bits361=361,preoperation_bits364=364,hard518=518,
  whole_owned32290=32290,products65000=65000,RAW_hard4000000=4000000,
  denominator='numeric D*5^18, symbolic a^18',no_binomial_Rpower_or_producer_engine=True),indent=2)+'\n')
 H=dict(groups[18]);N={};maximum=0;peak=0;products=0;cursor={}
 def account():
  nonlocal peak
  n=2*(4811+len(H)+len(N))+16000;peak=max(peak,n);assert n<=32290
 def add(key,v):
  nonlocal maximum
  a,e=divmod(key,1024);assert 0<=a<=34 and 0<=e<=51 and a+e<=68
  maximum=max(maximum,abs(v).bit_length());assert abs(v).bit_length()<=364
  value=N.get(key,0)+v;maximum=max(maximum,abs(value).bit_length())
  if abs(value).bit_length()>364:
   (OUT/'guard-prefix.json').write_text(json.dumps(dict(cursor=cursor,uncommitted=str(value),H={str(k):str(v)for k,v in H.items()},N={str(k):str(v)for k,v in N.items()},input_sha256=pins,no_auto_retry=True),separators=(',',':'))+'\n');raise RuntimeError('source364guard')
  if value:N[key]=value
  else:N.pop(key,None)
  account()
 for l in range(17,-1,-1):
  N={}
  for key,c in H.items():
   cursor={'phase':'literal_root_times_H','d':l,'key':key};add(key,4*c);add(key+1,6*c);products+=2
  scale=5**(18-l);shift=1024*(18-l)
  for key,c in groups[l].items():
   cursor={'phase':'scaled_original_input','d':l,'key':key};add(key+shift,scale*c);products+=1
  H=N;N={};account()
 assert len(H)==606 and products<65000 and not(PRIMARY/'completed-zero-orders.txt').read_bytes()
 with(PRIMARY/'first-nonzero-J-Hasse.txt').open()as f:
  assert f.readline().strip()=='Q1JHV001';k,den,external_a=map(int,f.readline().split());assert k==0 and den==D*5**18 and external_a==18
  assert int(f.readline())==606
  for key in sorted(H):
   a,e,c=map(int,f.readline().split());assert key==1024*a+e and H[key]==c
  assert not f.read().strip()
 facts=dict(J_valuation_exact0=0,all606_integer_coefficients_independently_equal=True,numeric_denominator_equal=True,symbolic_a18_denominator_equal=True,
  zero_orders_empty=True,method='descending-d integer Horner on original4811 coefficients, literal two-term R and shifted powers5a',
  no_producer_imports_binomial_or_root_power_engine=True,products=products,maximum_integer_bits=maximum,owned_upper=peak,no_J_quotient_or_GCD=True)
 (OUT/'verification.json').write_text(json.dumps(facts,indent=2)+'\n');pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 raw=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in OUT.iterdir())+20000;assert raw<4000000
 result=dict(id=OUT.name,all_passed=True,frozen=True,independent_actual_layer_reconstruction=True,input_sha256=pins,
  evidence_sha256={p.name:sha(p)for p in OUT.iterdir()},RAW_total_upper=raw,seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(facts,indent=2));print('report_sha256',sha(OUT/'report.json'))
if __name__=='__main__':main()

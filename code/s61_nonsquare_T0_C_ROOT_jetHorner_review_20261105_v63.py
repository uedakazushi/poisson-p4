#!/usr/bin/env python3
"""Independent five-jet integer Horner replay of T0/C orders0..4."""
import hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C_HASSE_WORKER_CUSTODY_20261105_v63'
PRIMARY=ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C_HASSE_ENTRY_20261105_v63'
OUT=ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C_ROOT_JET_HORNER_REVIEW_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));OUT.mkdir(exist_ok=False)
 pins={};p=PRIMARY/'report.json';assert sha(p)=='05172383bada31e9354a618d5236bb3c01fe3f313f11e27c37e04e397ce6469a'
 r=json.loads(p.read_text());assert r['all_passed'] and r['result']['C_valuation']==4;pins[str(p.relative_to(ROOT))]=sha(p)
 for name,want in r['evidence_sha256'].items():assert sha(PRIMARY/name)==want;pins[str((PRIMARY/name).relative_to(ROOT))]=want
 p=BASE/'report.json';assert sha(p)=='d96426dfc0506657a31817eaba3f50672715effa9630258b20a1b456b58deff4';pins[str(p.relative_to(ROOT))]=sha(p)
 c=json.loads(p.read_text());p=BASE/'native-input.txt';assert sha(p)==c['evidence_sha256'][p.name];pins[str(p.relative_to(ROOT))]=sha(p)
 groups=[{}for _ in range(21)];l1=0;n=0
 with p.open()as f:
  assert f.readline().strip()=='T0CH001 7208 20 0';D=int(f.readline())
  for line in f:
   key,value=map(int,line.split());a,rem=divmod(key,1048576);e,d=divmod(rem,1024)
   assert 0<=a<=20 and 0<=e<=40 and 0<=d<=18 and a+e+d<=57 and value
   assert 64*e+d not in groups[a];groups[a][64*e+d]=value;l1+=abs(value);n+=1
 assert n==7208 and l1<2**326 and D<2**197
 # V_i=(R+100z)V_(i+1)+100^(20-i)I_i, truncatedmodz^5.
 # This is a coefficient-major Horner engine, no source binomial or R^i.
 assert 2*(7208+2*5*1539)+16000==61196
 assert 2**326*129**20<2**481
 root_products=5*20*1539*3;z_products=4*20*1539;scalar_products=7208
 assert root_products+z_products+scalar_products<600000
 (OUT/'pre-arithmetic-source-and-budget.json').write_text(json.dumps(dict(
  recurrence='V_i=(2+9e-18e^2+100z)*V_(i+1)+100^(20-i)I_i mod z^5',
  integer_partial_bits481=481,integer_hard518=518,whole_owned61196=61196,
  all_products_upper600000=600000,RAW_hard4000000=4000000,
  fixed_input_already_T0_div_a4=True,no_second_shift=True,no_new_factor_or_quotient=True),indent=2)+'\n')
 H=[dict(groups[20])]+[{}for _ in range(4)];N=[{}for _ in range(5)];maximum=0;peak=0;products=0;cursor={}
 def account():
  nonlocal peak
  size=2*(7208+sum(map(len,H))+sum(map(len,N)))+16000;peak=max(peak,size);assert size<=61196
 def add(j,key,v):
  nonlocal maximum
  e,d=divmod(key,64);assert 0<=e<=80 and 0<=d<=18
  maximum=max(maximum,abs(v).bit_length());assert abs(v).bit_length()<=481
  value=N[j].get(key,0)+v;maximum=max(maximum,abs(value).bit_length())
  if abs(value).bit_length()>481:
   (OUT/'guard-prefix.json').write_text(json.dumps(dict(cursor=cursor,uncommitted=str(value),H=[{str(k):str(v)for k,v in p.items()}for p in H],N=[{str(k):str(v)for k,v in p.items()}for p in N],input_sha256=pins,no_auto_retry=True),separators=(',',':'))+'\n');raise RuntimeError('source481guard')
  if value:N[j][key]=value
  else:N[j].pop(key,None)
  account()
 for i in range(19,-1,-1):
  N=[{}for _ in range(5)]
  for j in range(5):
   for key,value in H[j].items():
    cursor={'phase':'root-times-oldjet','a':i,'jet':j,'key':key}
    for shift,scalar in [(0,2),(64,9),(128,-18)]:add(j,key+shift,value*scalar);products+=1
   if j:
    for key,value in H[j-1].items():
     cursor={'phase':'z-times-oldjet','a':i,'jet':j,'key':key};add(j,key,100*value);products+=1
  scalar=100**(20-i)
  for key,value in groups[i].items():
   cursor={'phase':'scaled-input','a':i,'key':key};add(0,key,scalar*value);products+=1
  H=N;N=[{}for _ in range(5)];account()
 assert not any(H[:4]) and len(H[4])==646 and products<600000
 assert (PRIMARY/'completed-zero-orders.txt').read_text().splitlines()==[f'{i} ZERO'for i in range(4)]
 with (PRIMARY/'first-nonzero-C-Hasse.txt').open()as f:
  assert f.readline().strip()=='T0CHV001';k,den=map(int,f.readline().split());assert k==4 and den==D*100**20
  assert int(f.readline())==646
  for key in sorted(H[4]):
   e,d,c=map(int,f.readline().split());assert 64*e+d==key and H[4][key]==c
  assert not f.read().strip()
 facts=dict(C_valuation_exact4=4,all_four_zero_layers_independently_rebuilt=True,all646_nonzero_layer_coefficients_equal=True,
  denominator_equal=True,method='truncated five-jet descending-a integer Horner with literal R+100z',
  no_producer_imports_or_binomial_or_R_power_engine=True,maximum_integer_bits=maximum,products=products,owned_upper=peak,
  no_quotient_or_GCD_evaluated=True)
 (OUT/'verification.json').write_text(json.dumps(facts,indent=2)+'\n');pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 raw=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in OUT.iterdir())+20000;assert raw<4000000
 result=dict(id=OUT.name,all_passed=True,frozen=True,independent_actual_value_replay=True,input_sha256=pins,
  evidence_sha256={p.name:sha(p)for p in OUT.iterdir()},RAW_total_upper=raw,seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(facts,indent=2));print('report_sha256',sha(OUT/'report.json'))
if __name__=='__main__':main()

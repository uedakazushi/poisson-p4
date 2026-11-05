#!/usr/bin/env python3
"""Distinct integer Horner reconstruction of the first F layer, no binomial/root powers."""
import hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_F_HORNER_REVIEW_SOURCE_20261105_v63'
PRIMARY=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_F_HASSE_ENTRY_20261105_v63'
BASE=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_F_HASSE_WORKER_CUSTODY_20261105_v63'
OUT=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_F_ROOT_HORNER_REVIEW_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 OUT.mkdir(exist_ok=False);pins={};products=0;maximum=0;owned=0;cursor={};H={};N={};groups=[{}for _ in range(17)]
 for p,h in [(SOURCE/'report.json','8f7de1c559889b6d9e10c1509bed75bcba3923fc6f7cb88f7c5a0da15bd07add'),(PRIMARY/'report.json','76e6e8eeb4cf89410ec5885caf7c7a935bee143308108f076002ced96f37f153')]:
  assert sha(p)==h;d=json.loads(p.read_text());assert d['all_passed'];pins[str(p.relative_to(ROOT))]=h
  for name,want in d.get('evidence_sha256',{}).items():
   q=p.parent/name;assert sha(q)==want;pins[str(q.relative_to(ROOT))]=want
 p=BASE/'native-input.txt';pins[str(p.relative_to(ROOT))]=sha(p)
 with p.open()as f:
  assert f.readline().strip()=='Q1FH001 4811 16 0';D=int(f.readline());previous=-1;count=0;l1=0
  for line in f:
   key,c=map(int,line.split());assert previous<key and c;previous=key;a,rem=divmod(key,1048576);e,dd=divmod(rem,1024)
   assert 0<=a<=16 and 0<=e<=33 and 0<=dd<=18 and a+e+dd<=50
   groups[a][e*64+dd]=c;count+=1;l1+=abs(c)
 assert count==4811 and l1<2**301 and 0<D<2**197
 assert 2*(4811+1254+1254)+16000==30638
 assert 2**301*500**16<2**445
 def account():
  nonlocal owned
  n=2*(sum(map(len,groups))+len(H)+len(N))+16000;owned=max(owned,n);assert n<=30638
 def add(key,v):
  nonlocal maximum
  e,dd=divmod(key,64);assert 0<=e<=65 and 0<=dd<=18
  b=abs(v).bit_length();maximum=max(maximum,b);assert b<=445
  value=N.get(key,0)+v;b=abs(value).bit_length();maximum=max(maximum,b)
  if b>445:
   receipt=dict(reason='445bit source guard',cursor=cursor,uncommitted=str(value),H={str(k):str(v)for k,v in H.items()},N={str(k):str(v)for k,v in N.items()},input_sha256=pins,no_auto_retry=True)
   (OUT/'guard-prefix.json').write_text(json.dumps(receipt,separators=(',',':'))+'\n');raise RuntimeError('guard')
  if value:N[key]=value
  else:N.pop(key,None)
  account()
 # Descending-a Horner uses a fresh accumulator per step, no expansion of R^i.
 H=dict(groups[16]);account()
 for i in range(15,-1,-1):
  N={}
  for key,c in H.items():
   cursor={'phase':'Horner-root-multiply','a':i,'key':key}
   for shift,scale in [(0,24),(64,-3),(128,-90)]:add(key+shift,c*scale);products+=1
  scalar=500**(16-i)
  for key,c in groups[i].items():
   cursor={'phase':'Horner-scaled-input','a':i,'key':key};add(key,c*scalar);products+=1
  H=N;N={};account()
 assert products<=65020 and len(H)==665
 del groups
 assert (PRIMARY/'completed-zero-orders.txt').read_bytes()==b''
 pins[str((PRIMARY/'first-nonzero-F-Hasse.txt').relative_to(ROOT))]=sha(PRIMARY/'first-nonzero-F-Hasse.txt')
 with (PRIMARY/'first-nonzero-F-Hasse.txt').open()as f:
  assert f.readline().strip()=='Q1FHV001';k,den=map(int,f.readline().split());assert k==0 and den==D*500**16
  n=int(f.readline());assert n==len(H);previous=None
  for key in sorted(H):
   e,dd,c=map(int,f.readline().split());assert key==64*e+dd and H[key]==c
  assert not f.read().strip()
 facts=dict(all665_coefficients_independently_equal=True,denominator_equal=True,zero_orders_empty=True,F_valuation_exact=0,
  method='descending-a integer Horner, 16 multiplication steps by literal 3-term R, with scaled input groups',
  no_producer_imports_or_binomial_or_R_power_engine=True,products=products,max_integer_bits=maximum,source_bits445=445,hard_bits518=518,
  whole_owned_upper=owned,whole_owned_hard30638=30638,only_F0_layer_evaluated=True,no_F_division_or_GCD=True)
 (OUT/'verification.json').write_text(json.dumps(facts,indent=2)+'\n')
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 raw=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in OUT.iterdir())+20000;assert raw<4000000
 result=dict(id=OUT.name,all_passed=True,frozen=True,independent_actual_value_reconstruction=True,
  input_sha256=pins,evidence_sha256={'verification.json':sha(OUT/'verification.json')},RAW_total_upper=raw,seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(facts,indent=2))
if __name__=='__main__':main()

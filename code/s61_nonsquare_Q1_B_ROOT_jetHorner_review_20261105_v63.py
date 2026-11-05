#!/usr/bin/env python3
"""Independent five-jet d-Horner with literal R+b1*z; no native S_h/binomial engine."""
import hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PRIMARY=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B_HASSE_ENTRY_20261105_v63'
SOURCE=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B_NEXT_SOURCE_20261105_v63'
ORIGINAL=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt'
OUT=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B_ROOT_JET_HORNER_REVIEW_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));OUT.mkdir(exist_ok=False);pins={}
 for p,want in [(PRIMARY/'report.json','b573db2d9d165cd311786faadf3f2dd1ca3488276604086bdd15e396e7d509a7'),(SOURCE/'report.json','14814743dfd4f4be1b320cf4b378c5f5ecad6700747bfe16f5f2c1c538a1bc4e')]:
  assert sha(p)==want;doc=json.loads(p.read_text());assert doc['all_passed'];pins[str(p.relative_to(ROOT))]=want
  for name,h in doc['evidence_sha256'].items():assert sha(p.parent/name)==h;pins[str((p.parent/name).relative_to(ROOT))]=h
 assert sha(ORIGINAL)=='d88cfe290ecbc71b045a382a333516815967ce6117576697c7cce78f2b00a6f7';pins[str(ORIGINAL.relative_to(ROOT))]=sha(ORIGINAL)
 groups=[{}for _ in range(19)];l1=0
 with ORIGINAL.open()as f:
  assert f.readline().strip()=='T1C4QN001';D=int(f.readline());assert int(f.readline())==4811
  for _ in range(4811):
   key,c=map(int,f.readline().split());a,rem=divmod(key,1048576);e,d=divmod(rem,1024)
   assert 0<=a<=16 and 0<=e<=33 and 0<=d<=18 and a+e+d<=50 and c
   groups[d][128*a+e]=c;l1+=abs(c)
  assert not f.read().strip()
 assert l1<2**301 and D.bit_length()==197
 # Two maps of 5-component jets charge six numeric cells per nonzero key.
 # Current/next b1 powers are charged separately. Allzero entries count too.
 assert 2*4811+2*2297*6+2*2*361+16000==54630<60000
 assert 2**301*179**18<2**436
 product_upper=18*5*2297*2+18*4*2297*4+4811*361+18*361*4
 assert product_upper<3000000
 (OUT/'pre-arithmetic-source-and-budget.json').write_text(json.dumps(dict(
  recurrence='V_l=(R+b1*z)V_(l+1)+b1^(18-l)I_l mod z^5',R='40-60e',b1='50a+18e^2-9e-2',
  source_bits436=436,preoperation437=437,hard518=518,whole_owned60000=60000,
  all_products_upper=product_upper,RAW_hard8000000=8000000,denominator='D197 times symbolic b1^18',
  no_binomial_native_Sh_cache_or_producer_import=True),indent=2)+'\n')
 H={key:[c,0,0,0,0]for key,c in groups[18].items()};N={};P={0:1};PN={};maximum=0;products=0;peak=0;cursor={}
 b1={128:50,2:18,1:-9,0:-2}
 def account():
  nonlocal peak
  n=2*4811+6*(len(H)+len(N))+2*(len(P)+len(PN))+16000;peak=max(peak,n);assert n<60000
 def guard(v):
  nonlocal maximum
  maximum=max(maximum,abs(v).bit_length())
  if abs(v).bit_length()>437:
   state=dict(cursor=cursor,uncommitted=str(v),H={str(k):list(map(str,x))for k,x in H.items()},N={str(k):list(map(str,x))for k,x in N.items()},P={str(k):str(x)for k,x in P.items()},input_sha256=pins,no_auto_retry=True)
   (OUT/'guard-prefix.json').write_text(json.dumps(state,separators=(',',':'))+'\n');raise RuntimeError('source437guard')
  return v
 def add(key,j,value):
  a,e=divmod(key,128);assert 0<=a<=34 and 0<=e<=69 and a+e<=86
  value=guard(value);old=N.get(key);q=guard((old[j]if old else 0)+value)
  if not old:
   if not q:return
   old=[0]*5;N[key]=old
  old[j]=q
  if not any(old):N.pop(key)
  account()
 for l in range(17,-1,-1):
  PN={}
  for key,c in P.items():
   for shift,x in b1.items():
    cursor={'phase':'current_b1_power','d':l,'key':key,'shift':shift};target=key+shift;v=guard(c*x);PN[target]=guard(PN.get(target,0)+v);products+=1
  PN={k:v for k,v in PN.items()if v};account();P=PN;PN={};N={}
  for key,vec in H.items():
   for j,c in enumerate(vec):
    if not c:continue
    cursor={'phase':'R_times_oldjet','d':l,'key':key,'jet':j};add(key,j,40*c);add(key+1,j,-60*c);products+=2
    if j<4:
     for shift,x in b1.items():
      cursor={'phase':'b1_z_times_oldjet','d':l,'key':key,'jet':j+1,'shift':shift};add(key+shift,j+1,c*x);products+=1
  for key,c in groups[l].items():
   for shift,x in P.items():
    cursor={'phase':'scaled_original_input','d':l,'key':key,'b1_power_key':shift};add(key+shift,0,c*x);products+=1
  H=N;N={};account()
 assert len(H)==1107 and all(v[:4]==[0,0,0,0]for v in H.values()) and products<3000000
 assert(PRIMARY/'completed-zero-orders.txt').read_text().splitlines()==[f'{i} ZERO'for i in range(4)]
 with(PRIMARY/'first-nonzero-B-Hasse.txt').open()as f:
  assert f.readline().strip()=='Q1BHV001';k,den,symbolic_b1=map(int,f.readline().split());assert k==4 and den==D and symbolic_b1==18
  assert int(f.readline())==1107
  for key in sorted(H):
   a,e,c=map(int,f.readline().split());assert key==128*a+e and H[key][4]==c
  assert not f.read().strip()
 facts=dict(B_valuation_exact4=4,all_four_zero_layers_independently_reconstructed=True,all1107_integer_coefficients_equal=True,
  numeric_D197_and_symbolic_b1_18_equal=True,method='descending-d 5-jet integer Horner, literal R+b1*z and one current b1 power',
  no_native_Sh_cache_binomial_or_producer_imports=True,products=products,maximum_integer_bits=maximum,owned_upper=peak,
  no_B_quotient_or_GCD=True)
 (OUT/'verification.json').write_text(json.dumps(facts,indent=2)+'\n');pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 raw=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in OUT.iterdir())+20000;assert raw<8000000
 result=dict(id=OUT.name,all_passed=True,frozen=True,independent_actual_B_layers_reconstruction=True,input_sha256=pins,
  evidence_sha256={p.name:sha(p)for p in OUT.iterdir()},RAW_total_upper=raw,seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(facts,indent=2));print('report_sha256',sha(OUT/'report.json'))
if __name__=='__main__':main()

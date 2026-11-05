#!/usr/bin/env python3
"""Independent fixed-C4 sparse integer longdivision of the original7208 QQ source."""
import hashlib,json,math,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PRIMARY=ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ENTRY_20261105_v63'
OLD=ROOT/'reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63/root239.poly'
OUT=ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C4_ROOT_LONGDIVISION_REVIEW_20261105_v63'
MA=1048576;ME=1024
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));OUT.mkdir(exist_ok=False)
 p=PRIMARY/'report.json';assert sha(p)=='a86f87ecddf58ae08110daf099b3fa61a0d290a018423550c22c9d68cd097131'
 r=json.loads(p.read_text());assert r['all_passed'] and r['result']['quotient_terms']==4500;pins={str(p.relative_to(ROOT)):sha(p)}
 for name,want in r['evidence_sha256'].items():assert sha(PRIMARY/name)==want;pins[str((PRIMARY/name).relative_to(ROOT))]=want
 qp=PRIMARY/'quotient-integer-common-D.txt'
 with qp.open()as f:assert f.readline().strip()=='T0C4QV001';D=int(f.readline())
 assert D.bit_length()==197 and sha(OLD)=='0b4bdc49ed1fc16711329358ac10d37f151d5af3864e00f479e180dcb2bf6755';pins[str(OLD.relative_to(ROOT))]=sha(OLD)
 R={};Q={};l1=0
 with OLD.open()as f:
  assert f.readline().strip()=='4' and f.readline().strip()=='7208'
  for _ in range(7208):
   a,e,d,text=f.readline().split();a,e,d=map(int,[a,e,d]);parts=list(map(int,text.split('/')));assert len(parts)in[1,2]
   num=parts[0];den=parts[1]if len(parts)==2 else 1;assert den>0 and D%den==0
   c=num*(D//den);key=a*MA+e*ME+d;assert key not in R and c
   assert 0<=a<=20 and 0<=e<=40 and 0<=d<=18 and a+e+d<=57
   R[key]=c;l1+=abs(c)
  assert not f.read().strip()
 assert len(R)==7208 and l1<2**326
 # Different construction: enumerate all multinomial choices of the FOUR
 # literal terms of C. No Laurent series, R powers, slab or 100^20 is used.
 C4={}
 for na in range(5):
  for ne2 in range(5-na):
   for ne1 in range(5-na-ne2):
    n0=4-na-ne2-ne1
    c=math.factorial(4)//(math.factorial(na)*math.factorial(ne2)*math.factorial(ne1)*math.factorial(n0))
    c*=100**na*18**ne2*(-9)**ne1*(-2)**n0
    key=na*MA+(2*ne2+ne1)*ME;C4[key]=C4.get(key,0)+c
 C4={k:v for k,v in C4.items()if v};assert len(C4)==25 and C4[4*MA]==100**4
 assert sum(map(abs,C4.values()))<=129**4<2**29
 # Independent finite formula/Gauss proof already bounds truequotient l1
 # by2^303, so any partial integer longdivision remainder is bounded by
 # ||I0||+||C4||*||J0||<2^333. Every mutation is checked before commit.
 assert 2**326+129**4*2**303<2**333
 assert 2*(14592+9690+55)+779+16000<70000
 (OUT/'pre-arithmetic-source-and-budget.json').write_text(json.dumps(dict(
  original_QQ_source7208=True,divisor_multinomial_terms=25,integer_partial_bound333=333,hard_integer518=518,
  whole_owned70000=70000,max_scalarproducts242250=242250,RAW_hard8000000=8000000,
  no_native_finite_series_or_scalar10020=True),indent=2)+'\n')
 peak=0;maximum=0;products=0;cursor={}
 def account(aliases=0):
  nonlocal peak
  owned=2*(len(R)+len(Q)+len(C4))+aliases+16000;peak=max(peak,owned);assert owned<70000
 def guard(value):
  nonlocal maximum
  maximum=max(maximum,abs(value).bit_length())
  if abs(value).bit_length()>333:
   (OUT/'guard-prefix.json').write_text(json.dumps(dict(cursor=cursor,uncommitted=str(value),R={str(k):str(v)for k,v in R.items()},Q={str(k):str(v)for k,v in Q.items()},input_sha256=pins,no_auto_retry=True),separators=(',',':'))+'\n');raise RuntimeError('source333guard')
  return value
 for a in range(20,3,-1):
  aliases=[k for k in R if k//MA==a];assert len(aliases)<=779;account(len(aliases))
  for key in aliases:
   c=R[key];assert c%(100**4)==0;q=c//(100**4);qkey=key-4*MA
   assert q and qkey not in Q;Q[qkey]=guard(q)
   for ck,cc in C4.items():
    cursor={'input_a':a,'quotient_key':qkey,'divisor_key':ck};target=qkey+ck
    aa,rem=divmod(target,MA);ee,dd=divmod(rem,ME);assert aa<=20 and ee<=40 and dd<=18 and aa+ee+dd<=57
    product=guard(q*cc);products+=1;value=guard(R.get(target,0)-product)
    if value:R[target]=value
    else:R.pop(target,None)
    account(len(aliases))
 assert not R and len(Q)==4500 and products==112500
 ql1=sum(map(abs,Q.values()));assert ql1<2**303
 with qp.open()as f:
  assert f.readline().strip()=='T0C4QV001' and int(f.readline())==D and int(f.readline())==4500
  previous=-1
  for _ in range(4500):
   key,c=map(int,f.readline().split());assert previous<key and Q[key]==c;previous=key
  assert not f.read().strip()
 facts=dict(all4500_integer_coefficients_independently_equal=True,full_longdivision_remainder_zero=True,
  original_D197_equal=True,quotient_l1_bits=ql1.bit_length(),method='fixed25-term C4 multinomial divisor, leading-a integer sparse longdivision on original7208 QQ source',
  no_producer_imports_finite_series_R_powers_or_scalar10020=True,products=products,maximum_integer_bits=maximum,owned_upper=peak,
  no_GCD_other_factor_or_tail=True)
 (OUT/'verification.json').write_text(json.dumps(facts,indent=2)+'\n');pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 raw=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in OUT.iterdir())+20000;assert raw<8000000
 result=dict(id=OUT.name,all_passed=True,frozen=True,independent_actual_quotient_reconstruction=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()},RAW_total_upper=raw,seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(facts,indent=2));print('report_sha256',sha(OUT/'report.json'))
if __name__=='__main__':main()

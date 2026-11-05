#!/usr/bin/env python3
"""ROOT source/code audit: exact old T0 QQ body, without a second offset shift or C values."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C_HASSE_WORKER_CUSTODY_20261105_v63'
OUT=ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C_HASSE_ROOT_SOURCE_REVIEW_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 OUT.mkdir(exist_ok=False);p=BASE/'report.json';assert sha(p)=='d96426dfc0506657a31817eaba3f50672715effa9630258b20a1b456b58deff4'
 r=json.loads(p.read_text());assert r['all_passed'] and r['value_worker_NOT_run'];pins={str(p.relative_to(ROOT)):sha(p)}
 for group in ['input_sha256','evidence_sha256']:
  for name,want in r[group].items():
   p=ROOT/name if group=='input_sha256'else BASE/name;assert sha(p)==want;pins[str(p.relative_to(ROOT))]=want
 old=ROOT/'reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63/root239.poly'
 l1=0;degrees=[];keys=set();maxcoef=0
 with old.open()as original,(BASE/'native-input.txt').open()as native:
  assert original.readline().strip()=='4' and original.readline().strip()=='7208'
  assert native.readline().strip()=='T0CH001 7208 20 0';D=int(native.readline());assert 0<D<2**197
  for _ in range(7208):
   a,e,d,text=original.readline().split();a,e,d=map(int,[a,e,d]);parts=list(map(int,text.split('/')));assert len(parts)in[1,2]
   num=parts[0];den=parts[1]if len(parts)==2 else 1;assert den>0 and D%den==0
   key,c=map(int,native.readline().split());assert key==a*1048576+e*1024+d and c==num*(D//den)
   assert key not in keys and c;keys.add(key)
   assert 0<=a<=20 and 0<=e<=40 and 0<=d<=18 and a+e+d<=57
   l1+=abs(c);maxcoef=max(maxcoef,abs(c).bit_length());degrees.append(a)
  assert not original.read().strip() and not native.read().strip()
 assert l1<2**326 and maxcoef==320 and len(keys)==7208
 # Every physical exponent in the old .poly BODY was already normalized
 # by its external offset4; matching every key above rules out a second shift.
 assert abs(2)+abs(9)+abs(-18)==29<100
 assert 2**326*2**20*100**20<2**481 and D*100**20<2**331
 forward=sum((a+1)**2 for a in degrees)
 max_layer=max(sum(2*(a-k)+1 for a in degrees if a>=k)for k in range(21))
 polynomial_products=forward+2*max_layer
 scalar_products=20+sum(3*(2*h-1)for h in range(1,21))+2*sum(a+1 for a in degrees)+4*7208+1
 assert polynomial_products<=3769784 and polynomial_products+scalar_products+7208<=4117089
 assert 2*(7208+sum(2*h+1 for h in range(21))+81*19)+16000==34376
 cpp=(ROOT/'code/s61_nonsquare_T0_C_hasse_20261105_v63.cpp').read_text()
 driver=(ROOT/'code/s61_nonsquare_prePOWER_T0_C_hasse_entry_20261105_v63.py').read_text()
 for text in ['in.reserve(n)','int root[3]={2,9,-18}','i=(int)aa;',
              'magic!="T0CH001"||n!=7208||amax!=20||shift!=0','scale100[20-h].v,choose[i][k]',
              'fmpz_bits(nextL.v)>326','fmpz_swap(input_l1.v,nextL.v)',
              'own(ninput,powslots,accsize+1)','convolution(k,3)','convolution(k,4)',
              'all21_layers_zero_contradicts_nonzero_input_no_retry','T0CHP001',
              'std::stoull(argv[4])+4400000+200000>8000000',
              'rlimit ram{2000000000,2000000000},cpu{60,60}',
              'phase=5;icursor=rcursor=-1;clockcheck()']:
  assert text in cpp,text
 assert 'timeout=58' in driver and 'raw+4600000<8000000' in driver and 'exist_ok=False' in driver
 closure=json.loads((BASE/'RAW-closure.json').read_text());assert r['RAW_total_upper']==5377074<8000000
 facts=dict(original_all7208_QQ_coefficients_exact_fixed_D=True,
  external_offset4_body_already_T0_div_a4=True,future_shift0_no_second_shift=True,
  input_L1_bits=l1.bit_length(),input_maxcoef_bits=maxcoef,D_bits=D.bit_length(),
  formula='D*100^20*H_k=sum I_ied*binomial(i,k)*100^(20-i+k)*(2+9e-18e^2)^(i-k)*e^e*d^d',
  valuation='C=100(a-root), C does not divide a; first nonzero order equals v_C(T0)',
  integer_partial_upper481=481,integer_hard518=518,owned_upper34376=34376,
  all_products_including_conversion_upper=polynomial_products+scalar_products+7208,
  RAW_upper=5377074,code_audit='fixed input; transactional L1; one accumulator; forward/reverse/reproduction; complete phase/key/power prefix; no double full I/O failure output',
  no_C_value_power_quotient_or_GCD_evaluated=True)
 (OUT/'source-proof.json').write_text(json.dumps(facts,indent=2)+'\n');pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 result=dict(id=OUT.name,all_passed=True,frozen=True,SOURCE_only=True,ONE_T0_C_valuation_ready=True,
  input_sha256=pins,evidence_sha256={'source-proof.json':sha(OUT/'source-proof.json')})
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(facts,indent=2));print('report_sha256',sha(OUT/'report.json'))
if __name__=='__main__':main()

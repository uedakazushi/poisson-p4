#!/usr/bin/env python3
"""Independent ROOT input/formula/code/budget audit; never evaluates F-Hasse values."""
import hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_F_HASSE_WORKER_CUSTODY_20261105_v63'
OUT=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_F_HASSE_ROOT_SOURCE_REVIEW_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 OUT.mkdir(exist_ok=False)
 report=json.loads((BASE/'report.json').read_text());assert report['all_passed'] and report['prepared_only'] and report['value_worker_NOT_run']
 pins={str((BASE/'report.json').relative_to(ROOT)):sha(BASE/'report.json')}
 for group in ['input_sha256','evidence_sha256']:
  for name,want in report[group].items():
   p=(ROOT/name)if group=='input_sha256'else(BASE/name);assert sha(p)==want,(name,want);pins[str(p.relative_to(ROOT))]=want
 native=(BASE/'native-input.txt').read_bytes().splitlines(keepends=True)
 original=(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt').read_bytes().splitlines(keepends=True)
 assert native[0]==b'Q1FH001 4811 16 0\n' and original[0]==b'T1C4QN001\n'
 assert original[2]==b'4811\n' and native[1]==original[1] and native[2:]==original[3:]
 D=int(native[1]);assert 0<D<2**197
 l1=0;degrees=[];previous=-1
 for line in native[2:]:
  key,c=map(int,line.split());assert previous<key and c;previous=key
  a,rem=divmod(key,1048576);e,d=divmod(rem,1024)
  assert 0<=a<=16 and 0<=e<=33 and 0<=d<=18 and a+e+d<=50
  l1+=abs(c);degrees.append(a)
 assert len(degrees)==4811 and l1<2**301
 # Absolute l1 majorant bounds every signed intermediate before arithmetic.
 # binomial(i,k)<=2^16; each 500^(16-h)*||R^h||_1<=500^16.
 assert abs(24)+abs(-3)+abs(-90)==117<500
 integer_majorant=2**301*2**16*500**16
 assert integer_majorant<2**463 and D*500**16<2**342
 forward=sum((a+1)**2 for a in degrees)
 max_layer=max(sum(2*(a-k)+1 for a in degrees if a>=k)for k in range(17))
 # At most two additional complete traversals of the first nonzero layer.
 polynomial_products=forward+2*max_layer
 scalar_products=16+sum(3*(2*h-1)for h in range(1,17))+2*sum(a+1 for a in degrees)+4*4811+1
 assert polynomial_products<=1707905 and polynomial_products+scalar_products<=1896419
 assert 2*(4811+sum(2*h+1 for h in range(17))+66*19)+16000==28708
 cpp=(ROOT/'code/s61_nonsquare_Q1_F_hasse_20261105_v63.cpp').read_text()
 py=(ROOT/'code/s61_nonsquare_prePOWER_Q1_F_hasse_entry_20261105_v63.py').read_text()
 for text in ['in.reserve(n)','std::stoull(argv[4])+4000000+200000>8000000','rlimit ram{2000000000,2000000000},cpu{60,60}',
              '500[16-h].v,choose[i][k]','int root[3]={24,-3,-90}','own(ninput,powslots,accsize+1)',
              'fmpz_bits(nextL.v)>301','fmpz_swap(input_l1.v,nextL.v)','convolution(k,3)','convolution(k,4)',
              'all17_layers_zero_contradicts_nonzero_input_no_retry','Q1FHP001','phase=5;icursor=rcursor=-1;clockcheck()']:
  assert text in cpp,text
 assert 'timeout=58' in py and 'raw+4200000<8000000' in py and "exist_ok=False" in py
 closure=json.loads((BASE/'RAW-closure.json').read_text());assert closure['total_upper']==4714969<8000000
 facts=dict(formula='D*500^16*H_k=sum K_ied*binomial(i,k)*500^(16-i+k)*(24-3e-90e^2)^(i-k)*e^e*d^d',
  coefficient_order_valuation='F=500*(a-root);500 is a rational unit, so first nonzero Hasse order equals exact F-adic valuation',
  normalized_input_body_byte_identical=True,input_terms=4811,input_l1_bits=l1.bit_length(),original_D_bits=D.bit_length(),
  all_integer_partials_strict_upper_bits=463,integer_hard_limit=518,denominator_bits_upper=342,
  whole_owned_upper=28708,polynomial_products_upper=polynomial_products,all_products_upper=polynomial_products+scalar_products,
  RAW_upper=4714969,RAW_hard_limit=8000000,
  code_audit='fixed input; transactional L1; exact key box; whole reserve; single accumulator; first nonzero reverse/reproduction; all-phase cursor prefix; no second full output on I/O failure',
  no_F_values_evaluated=True,no_general_factor_or_GCD=True)
 (OUT/'source-proof.json').write_text(json.dumps(facts,indent=2)+'\n')
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 result=dict(id=OUT.name,all_passed=True,frozen=True,SOURCE_only=True,ONE_F_valuation_job_ready=True,
  audit='ROOT independent literal input, formula, signed partial majorant, all products/ownership/RAW and frozen native/driver code',
  input_sha256=pins,evidence_sha256={'source-proof.json':sha(OUT/'source-proof.json')})
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(facts,indent=2))
if __name__=='__main__':main()

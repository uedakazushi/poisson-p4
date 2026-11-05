#!/usr/bin/env python3
"""ROOT literal/input/valuation/code/budget audit for the single named J."""
import hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_J_HASSE_WORKER_CUSTODY_20261105_v63'
OUT=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_J_HASSE_ROOT_SOURCE_REVIEW_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 OUT.mkdir(exist_ok=False);p=BASE/'report.json';assert sha(p)=='64ec22dc6f056fc0d71f4915c257bb53f8673608ad93f87170c30f3f0c8703aa'
 r=json.loads(p.read_text());assert r['all_passed'] and r['value_worker_NOT_run'];pins={str(p.relative_to(ROOT)):sha(p)}
 for group in ['input_sha256','evidence_sha256']:
  for name,want in r[group].items():
   p=ROOT/name if group=='input_sha256'else BASE/name;assert sha(p)==want;pins[str(p.relative_to(ROOT))]=want
 native=(BASE/'native-input.txt').read_bytes().splitlines(keepends=True)
 old=(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt').read_bytes().splitlines(keepends=True)
 assert native[0]==b'Q1JH001 4811 18 0\n' and old[0]==b'T1C4QN001\n' and old[2]==b'4811\n'
 assert native[1]==old[1] and native[2:]==old[3:];D=int(native[1]);l1=0;degrees=[]
 for line in native[2:]:
  key,c=map(int,line.split());a,rem=divmod(key,1048576);e,d=divmod(rem,1024)
  assert 0<=a<=16 and 0<=e<=33 and 0<=d<=18 and a+e+d<=50 and c
  l1+=abs(c);degrees.append(d)
 assert len(degrees)==4811 and l1<2**301 and D.bit_length()==197
 # For every l,k the correlated binomial/root/scalar weight is <=
 # 5^(18-l)*sum_k binomial(l,k)*5^k*10^(l-k)<=15^18.
 assert 2**301*15**18<2**372 and D*5**18<2**239
 for l in range(19):
  assert sum(math.comb(l,k)*5**(18-l+k)*10**(l-k)for k in range(l+1))==5**(18-l)*15**l
 support=[sum(1 for a in range(35)for e in range(52-k)if a+e<=68-k)for k in range(19)]
 assert support[0]==1667 and 2*(4811+190+1667)+16000==29336
 counts=[sum(d-k+1 for d in degrees if d>=k)for k in range(19)]
 eligible=[sum(d>=k for d in degrees)for k in range(19)]
 schedules=[sum(counts[:k+1])+2*counts[k]+2*(sum(eligible[:k+1])+2*eligible[k])+342+18+1 for k in range(19)]
 assert max(schedules)==441936<1200000
 cpp=(ROOT/'code/s61_nonsquare_Q1_J_hasse_20261105_v63.cpp').read_text()
 driver=(ROOT/'code/s61_nonsquare_prePOWER_Q1_J_hasse_entry_20261105_v63.py').read_text()
 for text in ['in.reserve(n)','int root[2]={4,6}','magic!="Q1JH001"||n!=4811||dmax!=18||shift!=0',
              'int h=d-k','scale5[18-h].v,choose[d][k]','(i+18-h)*1024+e+r',
              'degreea>34||e<0||e>51-kcur||degreea+e>68-kcur','fmpz_bits(nextL.v)>301',
              'fmpz_swap(input_l1.v,nextL.v)','convolution(k,3)','convolution(k,4)','Q1JHP001',
              "zz(Dclear)<<' '<<18",'all19_layers_zero_contradicts_nonzero_input_no_retry',
              'std::stoull(argv[4])+4000000+200000>8000000','rlimit ram{2000000000,2000000000},cpu{60,60}']:
  assert text in cpp,text
 assert 'D==D0*5**18 and am==18' in driver and 'timeout=58' in driver
 assert r['RAW_total_upper']==5615693<8000000
 facts=dict(unchanged4811_body_and_D197_exact=True,
  ordinary_valuation_proof='J=5ad-6e-4 is primitive in QQ[a,e][d], coefficients have gcd1; Gauss plus degree-one irreducibility gives prime J. J does not divide a, so localization at named unit a preserves J-adic order. First nonzero Taylor order at d=(6e+4)/(5a) equals v_J.',
  no_unit_ideal_claim_for_J_coefficients=True,
  exact_formula='P_k=sum I_ijl binomial(l,k)5^(18-l+k)a^(i+18-l+k)(4+6e)^(l-k)e^j',
  denominator='numeric D*5^18 times symbolic a^18',integer_partial_bits372=372,hard_bits518=518,
  whole_owned29336=29336,all_schedule_products_upper=441936,RAW_upper=5615693,
  code_audit='one accumulator, correlated exponent boxes beforecommit, transactional L1, first nonzero forward/reverse/reproduction, full symbolic denominator readback, allstate prefix, no double full I/O output',
  no_J_value_quotient_GCD_other_factor_or_tail_evaluated=True)
 (OUT/'source-proof.json').write_text(json.dumps(facts,indent=2)+'\n');pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 result=dict(id=OUT.name,all_passed=True,frozen=True,SOURCE_only=True,ONE_Q1_J_valuation_ready=True,
  input_sha256=pins,evidence_sha256={'source-proof.json':sha(OUT/'source-proof.json')})
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(facts,indent=2));print('report_sha256',sha(OUT/'report.json'))
if __name__=='__main__':main()

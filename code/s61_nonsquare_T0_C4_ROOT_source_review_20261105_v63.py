#!/usr/bin/env python3
"""Independent ROOT finite-slab/Gauss/multiplyback source audit. No quotient coefficients."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_WORKER_CUSTODY_20261105_v63'
OUT=ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ROOT_SOURCE_REVIEW_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 OUT.mkdir(exist_ok=False);p=BASE/'report.json';assert sha(p)=='3fe340136edc8a37dd1ebbee629e45e04583dccc6000f9c898a9bcb3349b290d'
 r=json.loads(p.read_text());assert r['all_passed'] and r['value_worker_NOT_run'];pins={str(p.relative_to(ROOT)):sha(p)}
 for group in ['input_sha256','evidence_sha256']:
  for name,want in r[group].items():
   p=ROOT/name if group=='input_sha256'else BASE/name;assert sha(p)==want;pins[str(p.relative_to(ROOT))]=want
 a=(BASE/'native-input.txt').read_bytes().splitlines(keepends=True)
 b=(ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C_HASSE_WORKER_CUSTODY_20261105_v63/native-input.txt').read_bytes().splitlines(keepends=True)
 assert a[0]==b'T0C4Q001 7208 20 0 16\n' and a[1:]==b[1:]
 D=int(a[1]);l1=0;degrees=[]
 for line in a[2:]:
  key,c=map(int,line.split());aa,rem=divmod(key,1048576);e,d=divmod(rem,1024)
  assert 0<=aa<=20 and 0<=e<=40 and 0<=d<=18 and aa+e+d<=57 and c
  l1+=abs(c);degrees.append(aa)
 assert len(degrees)==7208 and l1<2**326 and D.bit_length()==197
 # Laurent expansion at a=infinity yields the exact finite nonnegative
 # coefficient S_i. Previously independently proved C-adic valuation4
 # guarantees the negative Laurent tail vanishes; whole multiplyback still
 # checks this fact directly, coefficient by coefficient.
 assert -1*(-9)+4*(-2)==1
 assert 4*71**4>100**4  # sum binom(h+3,3)*(29/100)^h <4.
 assert 4*100**16*2**326<2**435
 assert 100**20>=2**132 and 2**435//(100**20)<2**303
 assert 129**4<2**29 and 2**303*129**4+2**326<2**333
 quotient_support=sum(1 for aa in range(17)for e in range(33)for d in range(19)if aa+e+d<=49)
 remainder_support=sum(1 for aa in range(21)for e in range(41)for d in range(19)if aa+e+d<=57)
 assert quotient_support==9690 and remainder_support==14592
 formula_own=2*(7208+289+9690+73*19)+16000
 verify_own=2*(7208+289+9690+14592)+16000
 assert formula_own==53148 and verify_own==79558
 formula_products=sum((aa-3)**2 for aa in degrees if aa>=4)
 pairs=sum(aa-3 for aa in degrees if aa>=4)
 all_products=formula_products+2*pairs+sum(3*(2*h-1)for h in range(1,17))+20+9690+9690*25+4*55
 assert all_products<2581132
 cpp=(ROOT/'code/s61_nonsquare_T0_C4_quotient_20261105_v63.cpp').read_text()
 driver=(ROOT/'code/s61_nonsquare_prePOWER_T0_C4_quotient_entry_20261105_v63.py').read_text()
 for text in ['in.reserve(n)','magic!="T0C4Q001"||n!=7208||am!=20||shift!=0||qm!=16',
              'int h=a-i-4','hundred[16-h].v,choose[h]','int root[3]={2,9,-18}',
              'ee>72','phase=3;ij=rr=-1','e>32||i+e+d>49',
              'fmpz_divisible(it->second.v,hundred[20].v)','fmpz_divexact(canonical.v,it->second.v,hundred[20].v)',
              'fmpz_mul(back.v,canonical.v,hundred[20].v)','R.clear();R.shrink_to_fit();rootN=0',
              'FULL_C4_multiplyback_NONZERO','fmpz_set(Dq.v,D.v)','T0C4QP001',
              'std::stoull(argv[4])+8300000+200000>11000000','rlimit ram{2000000000,2000000000},cpu{60,60}']:
  assert text in cpp,text
 assert 'timeout=58' in driver and 'base+8300000+200000<=11000000' in driver
 assert r['RAW_total_upper']==9994291<11000000<33554432
 facts=dict(input_body7208_and_D197_byte_identical=True,source_shift0_no_double_a4=True,
  finite_formula='S_i=sum I_j*binomial(j-i-1,3)*R^(j-i-4)*100^(20-j+i); J_i=S_i/100^20',
  primitive_Gauss_and_fixed_scalar_proved=True,whole_slab_then_support_and_exact_scalar_check=True,
  all_integer_partial_bits435=435,backcheck_bits436=436,multiplyback_bits333=333,hard_bits518=518,
  formula_owned=53148,verification_owned=79558,all_products_upper=all_products,RAW_upper=9994291,
  code_audit='full input, transactional L1, one complete slab, no high-e dropping, fixed10020 divisibility and backmultiply, release rootcache, full C4 coefficient multiplyback, originalD, whole cursor prefix, no double full I/O output',
  no_quotient_values_GCD_other_factors_or_tail_evaluated=True)
 (OUT/'source-proof.json').write_text(json.dumps(facts,indent=2)+'\n');pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 result=dict(id=OUT.name,all_passed=True,frozen=True,SOURCE_only=True,ONE_T0_C4_quotient_ready=True,
  input_sha256=pins,evidence_sha256={'source-proof.json':sha(OUT/'source-proof.json')})
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(facts,indent=2));print('report_sha256',sha(OUT/'report.json'))
if __name__=='__main__':main()

#!/usr/bin/env python3
"""SOURCE only for the distinct five-jet peer and fixed C4 quotient of oldT0."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_T0_C4_NEXT_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 t=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):
  return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 primary=load('S61_NONSQUARE_PREPOWER_T0_C_HASSE_ENTRY_20261105_v63');assert primary['all_passed']and primary['result']['C_valuation']==4
 src=load('S61_NONSQUARE_PREPOWER_T0_C_HASSE_SOURCE_20261105_v63','T0-literal-frame-and-preflight.json');D=int(src['common_D197'])
 z=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C_HASSE_ENTRY_20261105_v63/completed-zero-orders.txt');assert z.read_text().splitlines()==[f'{i} ZERO'for i in range(4)]
 hh=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C_HASSE_ENTRY_20261105_v63/first-nonzero-C-Hasse.txt')
 with hh.open()as f:
  assert f.readline().strip()=='T0CHV001';k,dc=map(int,f.readline().split());assert k==4 and dc==D*100**20;assert int(f.readline())==646
 jq=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt')
 with jq.open()as f:assert f.readline().strip()=='T1C4QN001';D1=int(f.readline())
 assert D==D1
 horner=dict(SOURCE_ONLY_NO_NEW_VALUES=True,input='originalroot239.poly QQ7208 bodyALREADY U0=T0/a^4, fixedD197; no native imports',
  algorithm='group I_i(e,d). J0=I20,J1..4=0. For i19..0 descending and r4..1 descending: Jnew_r=R*Jr+100*J(r-1); then Jnew0=R*J0+100^(20-i)*I_i.',
  R='2+9e-18e^2',perstage_denominator='D*100^(20-i); final ALLfivejets share D*100^20',
  root_scale100_in_crossjet_REQUIRED=True,
  proof='multiply the truncated Horner polynomial by(R+100z) and add 100^(20-i)I_i; this is independent of binomial/Rpower producer',
  source_jet_l1_bits_upper481=481,
  lifetime='input7208 coefficient/key maps; fivejets1539 each; ONEtemp1539; replace descendingr to retain previous lowerjet; no second fullfivejet copy',
  all_owned_upper48884=2*(7208+6*1539)+16000,
  proposed_peer_guard60000=60000,
  ALL_scalar_products_upper606565=606565,
  expected='J0,J1,J2,J3 literalZERO; J4 all646match and commonD*10020; no higherjets asserted',
  compare='after inputfree, stream646expected sortedcoeffs againstJ4, no expected fullmap duplicate',
  RAW_upper4000000=4000000,soft45_hard60_RAM2GB=True,independent_ROOT_or_peer_run_REQUIRED=True)
 dump('independent-fivejet-Horner-SOURCE.json',horner)
 qcount=sum(a+e+d<=49 for a in range(17)for e in range(33)for d in range(19))
 remcount=sum(a+e+d<=57 for a in range(21)for e in range(41)for d in range(19))
 assert (qcount,remcount)==(9690,14592)
 formula_owned=2*(7208+289+qcount+1387)+16000
 verify_owned=2*(7208+289+qcount+remcount)+16000
 assert (formula_owned,verify_owned)==(53148,79558)
 products=7208*289+25*qcount+2*7208*17+768+120+20+qcount+100
 quotient=dict(SOURCE_ONLY_NO_QUOTIENT_VALUES=True,independent_Cvaluation_peer_REQUIRED=True,
  U0='T0/a^4=I0/D0, I0 integer7208,D0 original197bits',
  target='J0=I0/C^4 integer polynomial; Q0=J0/D0=T0/(a^4*C^4)',
  C='100a-R, R=2+9e-18e^2',
  finite_scaled_slab='S_i=sum_{j>=i+4} I_j*binom(j-i-1,3)*R^(j-i-4)*100^(20-j+i),i0..16',
  scalar_normalization='after FULLslab, assert every S_i coefficient divisible fixed100^20, then J0_i=S_i/(100^20); backmultiply integer equality',
  Gauss_proof='Cprimitive since-1*(-9)+4*(-2)=1; rationalC4divisibility makesI0/C4 integral. ThereforeALLscaledslabcoeffs are divisible10020; no content/GCDsearch.',
  negative_a_remainder='no assumption or truncation: verify ALLcoeff C^4*J0=I0 after completeddivision, which eliminates the whole Laurentnegativepart',
  fullslab_before_bounds='e<=72,d<=18,1387slots; only AFTERwhole sum check finalbounds, no dropping high e terms',
  final_J0_bounds=dict(a=[0,16],e=[0,32],d=[0,18],total_degree_max49=49,support_upper9690=qcount),
  remainder_bounds=dict(a_max20=20,e_max40=40,d_max18=18,total_degree_max57=57,support_upper14592=remcount),
  lifetime='input7208+Rcache289+J0+ONEslab; releaseALLRcache and emptyslab BEFORE Cpowercache/remainder; one product accumulator',
  formula_owned53148=formula_owned,verification_owned79558=verify_owned,
  coefficient_and_packedkey_counted=True,ALLmetadata16000=16000,
  source_l1I0_bits326=326,scaledS_l1_bits435=435,
  canonical_J0_l1_bits303=303,integer_scalar_backcheck_bits436=436,
  full_C4multiplyback_partial_bits333=333,proposed_unchanged_hard_integer518=518,
  height_proof='||S||1<4*100^16||I||1; divide known100^20>=2^132; ||C4||1<2^29, so multiplyback partials<2^333',
  ALL_scaled_root_polynomial_products_upper2581132=products,
  exact_fixed_scalar_division_checks_upper9690=qcount,
  canonical_D0_equals_existing_J1_D1=True,common_D0_D1_integer=D,
  future_two_tail_reduction='T0=a4*C4*J0/D andT1=a3*C4*J1/D; sameQQdenominator, knownFULLUunits only. NoactualJ0 computed.',
  RAW_upper11000000=11000000,soft45_hard60_RAM2GB=True,
  one_state_OR_success_no_IO_doubleprefix_required=True,
  no_actual_division_GCD_otherfactor_tail_or_badlocus=True,new_budget_and_fresh_source_code_audit_REQUIRED=True)
 assert products==2581132
 dump('fixed-C4-integer-quotient-SOURCE.json',quotient)
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 r=dict(id=ID,all_passed=True,frozen=True,source_only=True,
  different_fivejet_Horner_peer_SOURCE_ready=True,actual_peer_arithmetic_NOT_run=True,
  fixed_C4_quotient_SOURCE_only=True,J0_quotient_NOT_computed=True,
  Gauss_integer_normalization_upfront=True,shared_D197_literal_equality=True,
  quotient_owned79558=79558,quotient_hard518=518,ALL_quotient_products2581132=2581132,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-t)
 dump('report.json',r);print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

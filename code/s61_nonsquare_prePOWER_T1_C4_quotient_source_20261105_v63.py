#!/usr/bin/env python3
"""SOURCE ONLY finite known-C^4 quotient operator; no quotient coefficient evaluated."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 t=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 v=load('S61_NONSQUARE_PREPOWER_T1_C_HASSE_ENTRY_20261105_v63');assert v['all_passed']and v['result']['C_valuation']==4
 s=load('S61_NONSQUARE_PREPOWER_T1_C_HASSE_WORKER_CUSTODY_20261105_v63');assert s['all_passed']
 qcount=sum(a+e+d<=50 for a in range(17)for e in range(34)for d in range(19))
 rcount=sum(a+e+d<=58 for a in range(21)for e in range(42)for d in range(19))
 slab=74*19;cache=17**2
 assert (qcount,rcount,slab,cache)==(10013,14991,1406,289)
 own_formula=2*(7595+cache+qcount+slab)+16000
 own_verify=2*(7595+cache+qcount+rcount)+16000
 assert max(own_formula,own_verify)==81776
 formula_pairs=7595*289;verify_pairs=25*qcount
 d=dict(only_SOURCE=True,independent_C4_valuation_peer_REQUIRED_before_any_actual_quotient=True,
  input='U1=T1/a^3=sum I_j(e,d)*a^j/D1; literal source7595, no old149/572',
  C='100*a-R(e)',R='-18*e^2+9*e+2',Q='U1/C^4',
  common_denominator='Dq=D1*100^20',
  exact_finite_operator='Dq*Q_i=sum_{j=i+4..20} I_j*binom(j-i-1,3)*R^(j-i-4)*100^(20-j+i), i0..16',
  coefficient_origin='Laurent polynomial part of (100a-R)^(-4)*U1; exact C^4 divisibility makes the negative remainder zero',
  no_Horner_repeated_rational_division_or_field_inverse=True,
  final_quotient_bounds=dict(a=[0,16],e=[0,33],d=[0,18],total_degree_max=50,support_upper=qcount),
  intermediate_slab_bound=dict(one_a_slice_only=True,e_max=73,d_max=18,support_upper=slab,
   no_cancellation_assumed_until_full_slice_done=True),
  outside_final_bounds='do not discard; a nonzero outside expected quotient support is a STOP before storage/claim',
  product_remainder_bound=dict(a_max=20,e_max=41,d_max=18,total_degree_max=58,support_upper=rcount),
  final_whole_coefficient_identity='C^4*Q_integer-100^20*U1_integer=0 in QQ[a,e,d]',
  constant_C4_support_upper=25,
  all_quotient_coefficients_from_completed_slices=True,
  known_R_power_slots=cache,
  lifetime='input7595 retained; build/stream Q_i one slab at a time, retain only verifiedQ; release all Rpower caches before remultiplication; bounded constant C-power cache<=289 slots',
  owned_formula=own_formula,owned_verification=own_verify,owned_MAX=81776,
  ownership='2 cells per coefficient+packedkey, ALLmetadata/scratch16000; no second full product',
  source_l1_bits363=363,
  quotient_l1_proof='sum_h binom(h+3,3)*(29/100)^h < (100/71)^4 <4 because4*71^4>100^4',
  exact_integer_inequality_4times71pow4=4*71**4,
  integer_Q_l1_bits_upper=472,
  remultiplication_partial_l1_bits_upper=502,
  proof_partial='L1(Qint)<4*100^16*L1(I); L1(C^4)<=129^4<2^29; L1(100^20 I)+L1(C^4 Qint)<2^502',
  proposed_hard_integer_bound518=518,
  all_formula_scalar_pairs_upper=formula_pairs,all_remultiplication_pairs_upper=verify_pairs,
  TOTAL_polynomial_products_upper=formula_pairs+verify_pairs,
  source_partial_product_guard_before_commit_required=True,
  proposed_soft45_hard60_RAM2GB=True,
  RAW_conservative_upper=15000000,
  RAW_scope='literal oneinput+ONE fullstate or final output+all source/cursors/receipts; commonD once, no duplicateQQ JSON',
  no_actual_quotient_GCD_factorization_other_source_or_tail=True,
  new_parent_budget_required=True,
  geometric_effect='same FULLU zero locus after named-unit a^3*C^4 removal; no bad-locus elimination or strict chart enlargement')
 assert d['exact_integer_inequality_4times71pow4']>100**4
 (out/'finite-operator-and-preflight.json').write_text(json.dumps(d,indent=2)+'\n')
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 r=dict(id=ID,all_passed=True,frozen=True,source_only=True,finite_quotient_recipe_complete=True,
  actual_C4_quotient_NOT_computed=True,independent_Cvaluation_peer_pending=True,
  whole_owned_upper81776=81776,integer_bits518=518,polynomial_products_upper=formula_pairs+verify_pairs,
  RAW_upper15000000=15000000,input_sha256=pins,evidence_sha256={'finite-operator-and-preflight.json':sha(out/'finite-operator-and-preflight.json')},seconds=time.monotonic()-t)
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

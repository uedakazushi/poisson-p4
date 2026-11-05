#!/usr/bin/env python3
"""SOURCE ONLY for a different-arithmetic F0 Horner peer; does not evaluate it."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_Q1_F_HORNER_REVIEW_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 t=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):
  return json.loads(pin(ROOT/'reports'/st/name).read_text())
 primary=load('S61_NONSQUARE_PREPOWER_Q1_F_HASSE_ENTRY_20261105_v63');assert primary['all_passed']and primary['result']['F_valuation']==0
 old=load('S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63');assert old['all_passed']and old['normalized_terms']==4811
 inp=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt')
 expected=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_F_HASSE_ENTRY_20261105_v63/first-nonzero-F-Hasse.txt')
 zero=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_F_HASSE_ENTRY_20261105_v63/completed-zero-orders.txt');assert zero.stat().st_size==0
 with inp.open()as f:
  assert f.readline().strip()=='T1C4QN001';D=int(f.readline());assert D.bit_length()==197;assert int(f.readline())==4811
 with expected.open()as f:
  assert f.readline().strip()=='Q1FHV001';k,dd=map(int,f.readline().split());assert k==0 and dd==D*500**16;assert int(f.readline())==665
 own=2*(4811+2*1254)+16000;assert own==30638
 d=dict(SOURCE_ONLY_NO_HORNER_VALUE_RUN=True,
  independent_input='ROOT original-common-D normalized quotient, no producer native or binomial/power import',
  group_input='17 coefficient maps I_i(e,d),4811total; packed e*64+d key; do not retain second raw input vector',
  R_F='24-3*e-90*e^2',cleared_denominator='D0*500^16',
  alternative_algorithm='H=I16; for i15..0: H_new=R_F*H+500^(16-i)*I_i; discard oldH before nextstep',
  exact_prefix_formula='H_i=sum_{j=i..16} I_j*R_F^(j-i)*500^(16-j)',
  proof='H_0 is exactly D0*500^16*Q_C(R_F/500,e,d), the producerk0 layer; no binomial/Taylor derivative machinery',
  source_signed_majorant='||H_i||1<=L1(allI)*500^(16-i), since||R_F||1=117<500',
  integer_partial_bits_upper445=445,hard_integer518_compatible=True,
  all17_input_groups_total4811=4811,H_old1254=1254,H_new1254=1254,
  ALL_metadata16000=16000,whole_owned30638=own,proposed_guard32000=32000,
  coefficient_and_packed_key_cells_counted=True,
  first_layer_expected_comparison='free inputgroups afterH0; stream expected665 sortedcoeffs againstH0; no second expected665map',
  polynomial_scalar_products_upper=16*1254*3+4811+16+1,
  all_zero_orders_to_check='none; completed-zero-orders.txt must be literally empty',
  exact_output_binding='all665 coefficients plus denominator D0*500^16; nonzeroHH proves v_F(Q_C)=0',
  no_qq_crossproducts_polynomial_GCD_point_prime_or_other_factor=True,
  RAW_upper4000000=4000000,soft45_hard60_RAM2GB=True,
  independent_root_or_peer_execution_REQUIRED=True,
  NO_F_division_source='factor not present; no F quotient is to be built')
 (out/'Horner-source-plan.json').write_text(json.dumps(d,indent=2)+'\n')
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 r=dict(id=ID,all_passed=True,frozen=True,source_only=True,Horner_values_NOT_computed=True,
  only_existing_F0_peer_target=True,no_F_division_source_because_valuation0=True,
  whole_owned30638=30638,products_upper65020=65020,integer_bits445=445,RAW_upper4000000=4000000,
  input_sha256=pins,evidence_sha256={'Horner-source-plan.json':sha(out/'Horner-source-plan.json')},seconds=time.monotonic()-t)
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

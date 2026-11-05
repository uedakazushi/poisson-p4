#!/usr/bin/env python3
"""Literal allorder gB tail/order/selector audit; no high coefficient read."""
import json
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_GB_PLANE3_CACHE_SOURCE_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 pins={}
 def read(p,keys=None):pins[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and pins[rp]['sha256']=='5a5d16ebc5b75fc8047e9180d6b5d757650c748a07944c89d05b35c326db9c6d'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/literal-gB-tail-orders-and-existing-source3-source2-selectors.json';d=read(ep);assert pins[ep]['sha256']==r['evidence_sha256'][Path(ep).name]
 chain='reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63';cr=read(chain+'/report.json');assert cr['all_passed'] and cr['actual_chain_identity_UP_QB']
 hist=read('reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW/report.json');assert hist['all_passed'] and hist['input_sha256'][chain+'/report.json']==pins[chain+'/report.json']['sha256']
 op=chain+'/actual-UP-QB-640-typed-operator.json';identity=read(op,['accepted_actual_regular_function_identities'])['accepted_actual_regular_function_identities'];assert hist['input_sha256'][op]==pins[op]['sha256']
 tail=[v for v in identity if v['kind']=='accepted_actual_CB_f63_tail'];assert tail==[{'kind':'accepted_actual_CB_f63_tail','function':28539,'gB':299650,'coordinate_quotient_nodes':[327551,327554,327557,327560,327563,327566,327569]}]
 rel='reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json';q=read(rel,['CB_coordinate_Y9_nodes','CB_f63_node','CB_g_restriction_node','CB_f63_zero7_telescoping_coefficients','all_coordinate_divisions']);z=q['CB_coordinate_Y9_nodes'][:7];h=q['CB_f63_zero7_telescoping_coefficients'];assert z==[27897,27948,28203,28254,28437,28692,28794] and h==tail[0]['coordinate_quotient_nodes'] and q['CB_f63_node']==28539 and q['CB_g_restriction_node']==299650
 # The full-relative accepted peer is a literal input to the successful
 # negative-boundary allorder proof, not a new interpretation of proof text.
 for key in ['relative','relative_peer','C0_peer']:
  p=cr['input_hashes'][key]['path'];assert sha(ROOT/p)==cr['input_hashes'][key]['sha256']
  if key.endswith('peer'):assert read(p)['all_passed']
 div={v['node']:v for v in q['all_coordinate_divisions']};assert all(div[k]['coordinate']==j for j,k in zip(z,h,strict=True))
 first={i:sp(v) for i,v in read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json',['first640'])['first640']['all79_atom_values_full10']};rows=[]
 for zi,hi in zip(z,h,strict=True):
  assert all(first[zi].get(i,Q(0))==0 for i in [0,8,9]) and first[hi].get(0,Q(0))==0
  rows.append({'zB_node':zi,'h_node':hi,'zB_point_r_ell_exact0':True,'h_point_exact0':True,'known_h_r_first':str(first[hi].get(8,Q(0))),'known_h_ell_first':str(first[hi].get(9,Q(0)))})
 assert rows==d['ordered_source_rows']
 source='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63';ids={i for i,v in read(source+'/actual-source105-ell2-raw33-values.json',['all_requested164_root_values'])['all_requested164_root_values']};assert {28539,*z}<=ids and read(source+'_REVIEW/report.json')['all_passed']
 raw=read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',['actual33_basis'])['actual33_basis'];plane=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json',['plane9_monomials','raw33_plane9_map_columns']);assert read('reports/AMBIENT_B18_ELL2_T2_UPPER617_BOUNDARY_PLANE_SOURCE_20261105_v63_REVIEW/report.json')['all_passed']
 selectors=[]
 for ri,pi,m in [(24,3,[2,0]),(25,4,[1,1]),(30,6,[2,1])]:
  assert raw[ri]==[0,*m] and plane['plane9_monomials'][pi]==m and plane['raw33_plane9_map_columns'][ri]==[[pi,'1']]
  assert [i for i,col in enumerate(plane['raw33_plane9_map_columns']) if any(j==pi and Q(x) for j,x in col)]==[ri];selectors.append([m,ri,pi])
 assert selectors==d['source_selector_triplets'] and d['zB3_h0_and_all_h2_h3_terms_zero_by_source_order']
 proof='The accepted exact regular-function identity is fCB-gB=sum(zB_i*h_i) in the actual fullR0 frame. Each zB_i has plane order>=2 and each h_i has plane order>=1. In ordinary coefficient r^2ell, only zB_rr*h_ell+zB_rell*h_r survive. Thus gB_r2ell=fCB_r2ell-sum(those two terms), with no zB3,h2,h3. Literal164 and the certified raw33-plane9 ring map make the selected fCB3/zB2 slots available, but this source audit neither extracts nor multiplies them.'
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ve=out/'independent-negative-boundary-gB-tail-order-and-raw-selector-proof.json';ve.write_text(json.dumps({'accepted_successful_NEGATIVE_BOUNDARY_chain_only':True,'literal_allorder_tail_z7_h7_pairing_exact':True,'all7_zB_point_and_plane_firsts_zero':True,'all7_h_points_zero':True,'only_known_h_firsts_referenced':rows,'literal164_fCB_and_zB7_membership':True,'raw33_plane9_unique_rr_rell_r2ell_selectors':selectors,'proof':proof,'no_second_third_scalar_numbers_extracted_or_multiplied':True},separators=(',',':'))+'\n');rr={'id':STEM,'all_passed':True,'primary_sha256':pins[rp]['sha256'],'source_identity_orders_and_availability_only':True,'no_actual_gB3_or_upper_response':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ve.name:sha(ve)}};(out/'report.json').write_text(json.dumps(rr,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()

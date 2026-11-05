#!/usr/bin/env python3
"""Independent complete first Hilbert layer from peered source values; fixed rank only."""
from fractions import Fraction as Q
from pathlib import Path
import json,hashlib,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,mv
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_FIRST_HILBERT_RANK_20261105_v63';STEM=A+'_REVIEW'
G='AMBIENT_B18_ELL2_T2_FIRST_HILBERT_RANK_INPUT_GATE_20261105_v63'
P={'primary':f'reports/{A}/report.json','value':f'reports/{A}/actual-two-by619-rank-and-complete-first-Hilbert-layer.json',
'gate_report':f'reports/{G}/report.json','gate':f'reports/{G}/actual-two-by619-first-Hilbert-complete-source-input.json',
'old':'reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63/actual104-source-rell-values-and-two-by619-old-knownparts.json',
'old_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63_REVIEW/report.json',
'native':'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json',
'native_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63_REVIEW/report.json',
'adjoint_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63_REVIEW/report.json',
'known':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/report.json',
'known_peer':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63_REVIEW/report.json',
'mini_peer':'reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'FIXED_SOURCE_RANK_AND_COMPLETE_FIRST_LAYER','inputs':inputs});D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()}
 assert all(D[k]['all_passed'] for k in D if k.endswith('peer') or k in ['primary','gate_report','known'])
 for k,stem in [('primary',A),('gate_report',G)]:
  for p,h in D[k]['input_sha256'].items():assert sha(ROOT/p)==h
  for p,h in D[k]['evidence_sha256'].items():assert sha(ROOT/'reports'/stem/p)==h
 E=D['value'];old=D['old']['all619_two_old_knownparts_including_frozen_pivot_first_correction'];native=D['native']['all619_two_actual_native_parts']
 assert old==D['gate']['actual_old_two_by619'] and native==D['gate']['actual_native_two_by619'] and len(old)==len(native)==619
 columns=[{k:Q(a)+Q(b) for k,(a,b) in enumerate(zip(x,y)) if Q(a)+Q(b)} for x,y in zip(old,native)];assert list(map(enc,columns))==E['all619_actual_two_class_columns']
 assert D['old']['all_honest_second_rell_sensitivities_exactly_zero'] and D['native']['L1ell_Omega_r_retained'] and D['native']['negative_intrinsic_pivot_correction_already_in_old_knownparts_not_added_again']
 # No fresh elimination: verify the supplied fixed minor and all619 factors.
 selected=list(map(sp,E['selected_original_columns']));rank=E['rank'];assert len(selected)==rank<=2
 for p,f in zip(columns,E['all619_original_basis_factorizations']):assert mv(selected,sp(f))==p
 if rank:
  minor=E['fixed_original_nonzero_minor'];chosen=minor['columns'];rows=minor['rows'];assert [columns[j] for j in chosen]==selected
  if rank==1:det=columns[chosen[0]][rows[0]]
  else:assert rows==[0,1];u,v=selected;det=u.get(0,Q(0))*v.get(1,Q(0))-u.get(1,Q(0))*v.get(0,Q(0))
  assert det==Q(minor['determinant']) and det
 else:assert all(not p for p in columns) and E['fixed_original_nonzero_minor'] is None
 assert D['known']['known_G3_top_linear_rank']==25 and D['known']['known_bottom_linear_rank']==23 and D['known']['zero_bottom_kernel_count']==619
 hf=[6,54-25-23-rank];assert hf==E['actual_HF_T2_through1']==D['primary']['actual_HF_T2_through1'] and sum(hf)==E['actual_length_T2_mod_m2']==D['primary']['actual_length_T2_mod_m2']
 ep=o/'independent-two-by619-rank-and-first-layer-completeness.json';dump(o,ep.name,{'inputs':inputs,'all619_old_native_additions_and_pivot_correction_applied_exactly_once':True,'all619_fixed_original_basis_factorizations_verified':True,'independent_fixed_nonzero_minor':E['fixed_original_nonzero_minor'],'actual_two_class_rank':rank,'actual_HF_T2_through1':hf,'actual_length_T2_mod_m2':sum(hf),'complete_linear_initial_relation_proof':'The accepted complete minimal6x1060 matrix has zero constants. Hence all degree1 relations are the constant span of its1060 columns, with no higher-coefficient contributions. The418 G3 columns occupy upper rank25. The642 remaining columns have bottom rank23, and exactly619 reviewed kernel combinations. Their remaining two classes are the exact full918 adjoint evaluations checked by old and native source peers. Therefore full first-relation rank is25+23+d.','actual_T2_total_length_H0_length_K_exponents_and_higher_HF_not_computed':True,'new_source_or_implicit_or_gauge_coefficients':0})
 report={'id':STEM,'all_passed':True,'actual_two_class_rank':rank,'actual_HF_T2_through1':hf,'actual_length_T2_mod_m2':sum(hf),'actual_total_T2_length_or_H0_exponents_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'rank':rank,'HF':hf}))
if __name__=='__main__':main()

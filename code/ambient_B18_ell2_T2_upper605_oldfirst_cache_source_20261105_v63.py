#!/usr/bin/env python3
"""Pruned existing OLD/FIRST/native92 cache bindings; no response contractions."""
import argparse,json,resource,signal,time
from pathlib import Path
from fractions import Fraction as F
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER605_OLD_FIRST_CACHE_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);pins={};raw=0
 def read(p):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024;return json.loads((ROOT/p).read_text())
 def field(p,k):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024;return H.selected(p,k)
 # Stage FIRST extraction before loading either large shared polynomial DAG.
 cp='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63/pruned-known-Cfull224-input-and-positive-source-envelopes.json';c=read(cp);piv=c['pivot27_source_indices'];assert len(piv)==27
 old=field(H.PATHS['old640'],'source_columns640_rows56_sparse_dual9');nat=field(H.PATHS['old128'],'source_columns128_rows48_sparse_dual9');recipes=[];u0=set();ue=set()
 assert count([c,old,nat])+15000<100000
 active=[0,1,4,5,7,8]
 for p,orig in enumerate(piv):
  for var in active:
   if var==8:
    kind,arg=c['pivot27_ell_recipes'][p]
    if kind=='literal':r={'kind':'literal56','column':arg}
    else:
     assert kind=='native' and arg==orig-640;r={'kind':'native_ell_SUM','native128_source':arg,'Se':c['native_Se'][arg],'S0':c['native_S0'][arg]};u0.update(k for k,q in r['Se']);ue.update(k for k,q in r['S0'])
   else:
    dual=8 if var==7 else var+1;cols=old[orig] if orig<640 else nat[orig-640];v=[]
    for row,values in cols:
     q=dict(values).get(dual,'0')
     if q!='0':v.append([row,q])
    if var==7:assert v==c['pivot27_r_columns'][p]
    r={'kind':'literal56','column':v,'original_old8_dual_index':dual,'native_bottom8_exact0':orig>=640}
   recipes.append({'pivot27':p,'original768':orig,'variable9':var,'first_column_recipe':r})
 del old,nat
 L0=[[k,c['native_L0'][k]] for k in sorted(u0)];Le=[[k,c['native_Le'][k]] for k in sorted(ue)];Gamma=c['Gamma_two'];del c
 invp='reports/AMBIENT_B18_ELL2_T2_UPPER605_SHARED_SOURCE_20261105_v63/';ir=read(invp+'report.json');assert ir['all_passed'] and pins[invp+'report.json']=='af0728e250dda8979857b60e1d414e19399f311f46bd756cb1382deb4df3cc89';inv=read(invp+'complete605-shared-source-linearity-and-global24-opaque-frontier.json')
 peer=read('reports/AMBIENT_B18_ELL2_T2_UPPER605_SHARED_SOURCE_20261105_v63_REVIEW/report.json');assert peer['all_passed']
 jp='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_INPUT_GATE_20261105_v63/617-direct-old-known-first-formal-Jacobian-positive-input.json';j=read(jp)
 sp='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/combined-OLD2-D-and-four-pivot-source-factor-inventory.json';s=read(sp)
 assert j['Gamma_two']==Gamma and j['original619_source768_combinations']==s['original619_c0_source_combinations']
 typed=dict(s['typed_program']);assert len(typed)==757 and all(typed[i]==n for i,n in j['typed_source_program'])
 aa=field('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json','all_zero_bottom_and_Y0_source_combinations');assert len(aa)==619
 assert all(r['original768_combination']==s['original619_c0_source_combinations'][i] and r['full56_Xgamma1_equals_Xp0a1'] for i,r in enumerate(aa));a1=[r['actual_a1_27_times9'] for r in aa];del aa
 recs={r['free_case619']:r['subtract_actual_basis563609'] for r in s['kernel617_source_recipes']};jf={r['free_case619']:r['constant_case619_source_factors'] for r in j['case_source_recipes617']}
 for i,(a,b) in recs.items():
  expected=[{'case619':i,'scale_operator':['literal','1']}]
  for base,q in [(563,a),(609,b)]:
   if q!='0':expected.append({'case619':base,'scale_operator':['negate_literal',q]})
  assert jf[i]==expected
 nativep='reports/AMBIENT_B18_NATIVE_KERNEL617_TO105_SOURCE_20261105_v63/all617-native-constant-sources-factor-through-fixed105-kernel.json';nm=read(nativep)
 assert all(r['subtract563609']==recs[r['free_case619']] for r in nm['all617_source_maps'])
 del nm
 used={col for v in s['original619_c0_source_combinations'] for col,q in v if col<640};assert used<=set(dict(s['D_old_root_columns'])) and used<=set(dict(j['pruned_original_old_columns_roots']))
 boundary=j['literal_boundary_firsts_point_r_ell_and_membership'];want={r[0] for r in boundary};assert len(want)==96
 f=field('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json','first640');f79={i:dict(v) for i,v in f['all79_atom_values_full10']};del f
 f=field('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json','all_requested164_root_values');physical={i:dict(v) for i,v in f if i in want};del f
 sel=inv['B21_raw33_selectors'];assert sel==[0,7,6,5,4,3,2,1,8,9,22,23,20,21,14,15,12,13,24,25,30]
 leaves=[];physical_count=0;only79_count=0
 for fn,kind,point_r_e in boundary:
  if fn in physical:
   v=[physical[fn].get(k,'0') for k in sel];physical_count+=1
   if fn in f79:assert list(map(F,v[:10]))==[F(f79[fn].get(k,'0')) for k in range(10)]
   mode='actual_physical164_B21'
  else:
   assert fn in f79;v=[f79[fn].get(k,'0') for k in range(10)]+['0']*11;only79_count+=1;mode='FULL9firsts_with_arbitrary_higher_extension_ONLY_after_whole605_certificate'
  assert list(map(F,[v[0],v[8],v[9]]))==list(map(F,point_r_e)) and F(v[0])==0
  leaves.append([fn,mode,v])
 assert physical_count==48 and only79_count==48
 basis=read('reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63/ONE508-source-closure-B21-and-exact-carry-zeros.json')
 assert basis['all660_NF_projection_identities_verified'] and basis['all1089_raw33_products_preserved'] and basis['raw33_to_B21_selector21']==sel
 n92=read('reports/AMBIENT_B18_NATIVE_PIVOT23_SECOND_VALUES_20261105_v63/actual92-Gamma-nativepivot23-plane2-and322-full-product-pieces.json');pr=read('reports/AMBIENT_B18_NATIVE92_ROOT_DIRECT_REVIEW_20261105_v63/report.json');assert pr['all_passed']
 sources92=n92['original768_native23_sources'];assert sources92==piv[4:];targets92=n92['Gamma_Xp2_native23_twochannels_rr_rell'];assert len(targets92)==23
 allrecipes=[[i,recs[i]] for i in inv['cases605']]
 data={'cases605':inv['cases605'],'opaque_cases12_kept_outside_value_scope':inv['opaque_cases12'],'actual_or_formal_B21_boundary96':leaves,'physical48_actual_and_only79_48_extension_labels_explicit':True,'all96_point_r_ell_and_physical_FULL9_overlap_checked':True,'ordinary_B21_basis':basis['ordinary_basis_B21'],'raw33_B21_selector':sel,'shared_literal_C0_program532':j['literal_C0_program'],'shared_full_typed_program757':s['typed_program'],'Gamma_two':Gamma,'D_two_by10_by56':s['D_two_by10_by56'],'OLD_Gamma_source_columns618':j['pruned_original_old_columns_roots'],'OLD_D_source_columns618':s['D_old_root_columns'],'raw619_original768_constant_sources':s['original619_c0_source_combinations'],'raw619_actual_a1_27times9':a1,'corrected605_recipes_basis563609':allrecipes,'point27_source_indices':piv,'individual_FIRST_column_recipes162':recipes,'native_L0_selected_columns':L0,'native_Le_selected_columns':Le,'accepted_native92_original_sources23':sources92,'accepted_native92_Gamma_plane_targets23_twochannels_rr_rell':targets92,'fixed_pair_order10':s['pair_order10'],'source_extension_scope':'The96 shared leaf polynomials include arbitrary higher extensions for48 only79 functions. They are not actual96 jets. Only the COMPLETE corrected605 OLD functional is actual by accepted all617 J1+OLD2 certificates and exactouterOLD4 theorem. Keep all757 FIRSTcross and all29/etc a1 entries.','raw619_partial_response_scope':'Raw OLD/FIRST/native92outer partial results are formal extensions, and must be corrected by subtracting actual source563/609 factors before interpreting605. No raw basisupper is set0. Native corrected605 uses separate true72kernel response map.','shared_FIRST_table_recipe':'T[t,p,v] sums D_pair dot Xp_first_other over pairs containingv; rr once. Then raw619 dot a1 and correct563609. All nativeell individualSUM retained.','OLD4_plane_outer_removed_ONLY_on_corrected605_by_exact_linear_12function_certificate':True,'native_kernel72_cache_SOURCE_external_not_duplicated':True,'all_Q_and_native_FINE_and_whole_owner_height_codec_custody_peer_ROOTGO_still_required':True,'no_source_polynomial_jet_convolution_or_scalar_response_contraction_computed':True}
 ep=out/'pruned605-OLD757-C0532-B21-FIRST162-and-native92-source-input.json';H.dump(ep,data)
 owned=count([j,s,inv,ir,peer,boundary,f79,physical,leaves,a1,recipes,L0,Le,n92,pr,basis,data,pins])+15000;assert owned<150000
 report={'id':STEM,'all_passed':True,'cache_binding_and_recipe_SOURCE_ONLY':True,'boundary96_physical48_actual_only79_48_formal_extension':True,'full757_and532_programs_kept':True,'individual_FIRST_recipes162':True,'all619_full_a1_and_exact563609_recipe_bindings':True,'native92_all23_reference_only_no_pairings':True,'source_owned_numeric_upper':owned,'SOURCE_RAW_read_upper':raw,'value_convolution_or_response_UNRUN':True,'future_positive_whole_source_code_custody_peer_ROOTGO_required':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','source_owned_numeric_upper','SOURCE_RAW_read_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

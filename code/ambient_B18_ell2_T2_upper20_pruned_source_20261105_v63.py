#!/usr/bin/env python3
"""Exact subset/ancestor metadata for20 exceptional605 classes; no values."""
import argparse,json,resource,signal,time
from pathlib import Path
from fractions import Fraction as F
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER20_PRUNED_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);pins={};raw=0
 def read(p):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024;return json.loads((ROOT/p).read_text())
 sp='reports/AMBIENT_B18_ELL2_T2_UPPER605_OLD_FIRST_CACHE_SOURCE_20261105_v63/';gr=read(sp+'report.json');assert gr['all_passed'];s=read(sp+'pruned605-OLD757-C0532-B21-FIRST162-and-native92-source-input.json')
 pp='reports/AMBIENT_B18_ELL2_T2_UPPER605_POSITIVE_WHOLE_20261105_v63/';pr=read(pp+'report.json');assert pr['all_passed'];p=read(pp+'positive-shared605-OLD-FIRST-nativeOuter-nativeFINE-and-Q-common-bounds.json')
 zr=read('reports/AMBIENT_B18_ELL2_T2_UPPER554_VALUES_V2_20261105_v63_REVIEW/report.json');assert zr['all_passed']
 cases=[i for i in p['structural_positive_nonzero_case_labels'] if i!=554];assert len(cases)==20
 zero584=sorted(set(s['cases605'])-set(cases)-{554});assert len(zero584)==584
 rec=dict(s['corrected605_recipes_basis563609']);needed=set(cases)
 for i in cases:
  for base,q in zip([563,609],rec[i]):
   if F(q):needed.add(base)
 assert needed==set(cases)|{563} and rec[617]==['-20','0']
 c0=s['raw619_original768_constant_sources'];a1=s['raw619_actual_a1_27times9'];old=sorted({j for i in needed for j,q in c0[i] if j<640})
 Grows={i for row in s['Gamma_two'] for i,q in row if F(q)};Drows={i for ch in s['D_two_by10_by56'] for row in ch for i,q in enumerate(row) if F(q)}
 gc=dict(s['OLD_Gamma_source_columns618']);dc=dict(s['OLD_D_source_columns618']);gc={j:[[row,tag] for row,tag in gc[j] if row in Grows] for j in old};dc={j:[[row,node] for row,node in dc[j] if row in Drows] for j in old}
 t=dict(s['shared_full_typed_program757']);C=dict(s['shared_literal_C0_program532']);keepT=set();keepC=set();boundary=set();todoT=[node for roots in dc.values() for row,node in roots];todoC=[]
 for roots in gc.values():
  for row,(kind,node) in roots:
   if kind=='typed':todoT.append(node)
   else:assert kind=='C0';todoC.append(node)
 while todoT:
  i=todoT.pop()
  if i in keepT:continue
  keepT.add(i);d=t[i]
  if d[0]=='actual_C0_scalar_function':boundary.add(d[1])
  elif d[0]=='sum':todoT.extend(k for k,a,b in d[1])
  else:assert d[0]=='product';todoT.extend(d[3])
 while todoC:
  i=todoC.pop()
  if i not in C:boundary.add(i);continue
  if i in keepC:continue
  keepC.add(i);d=C[i]
  if d[0]=='n':todoC.append(d[1])
  elif d[0] in ['a','m']:todoC.extend(d[1:3])
  else:assert d[0]=='c'
 leaves={i:(mode,v) for i,mode,v in s['actual_or_formal_B21_boundary96']};assert boundary<=set(leaves)
 pairs=s['fixed_pair_order10'];firstuses=set()
 for i in needed:
  for slot,q in a1[i]:
   pv,var=divmod(slot,9)
   for x,y in pairs:
    if var==x:firstuses.add((pv,y))
    if var==y and x!=y:firstuses.add((pv,x))
 recipes=[r for r in s['individual_FIRST_column_recipes162'] if (r['pivot27'],r['variable9']) in firstuses];assert {(r['pivot27'],r['variable9']) for r in recipes}==firstuses
 u0=set();ue=set()
 for r in recipes:
  d=r['first_column_recipe']
  if d['kind']=='native_ell_SUM':u0.update(i for i,q in d['Se']);ue.update(i for i,q in d['S0'])
 L0={i:v for i,v in s['native_L0_selected_columns']};Le={i:v for i,v in s['native_Le_selected_columns']};assert u0<=set(L0) and ue<=set(Le)
 piv=s['point27_source_indices'];needed92={piv[slot//9] for i in needed for slot,q in a1[i] if slot//9>=4 and slot%9 in [7,8]};old92=s['accepted_native92_original_sources23'];v92=s['accepted_native92_Gamma_plane_targets23_twochannels_rr_rell'];v92={orig:v for orig,v in zip(old92,v92)};assert needed92<=set(v92)
 qp='reports/AMBIENT_B18_ELL2_T2_UPPER617_FULL_COORDINATE_SOURCE_20261105_v63/known-Psi-point-correction-and-full-b1-source-annihilation.json';pins[qp]=H.sha(ROOT/qp);raw+=2*(ROOT/qp).stat().st_size;Q=H.selected(qp,'Q_two6')
 data={'source_metadata_ONLY_no_convolution_or_responses':True,'exceptional_cases20':cases,'required_raw_cases21':sorted(needed),'all605_coverage':{'SOURCE_positive_zero584':zero584,'independently_actual_zero554':554,'remaining_uncomputed20':cases},'584_ZERO_requires_fresh_independent_whole_positive_peer_before_adoption':True,'original_source_allc0_raw21':[[i,c0[i]] for i in sorted(needed)],'original_full_a1_raw21':[[i,a1[i]] for i in sorted(needed)],'corrected20_source_recipes':[[i,rec[i]] for i in cases],'source617_plus20_source563_carry_kept_raw563_NOT_actual_upper0':True,'boundary_cache':[[i,*leaves[i]] for i in sorted(boundary)],'literal_C0_program':[[i,d] for i,d in s['shared_literal_C0_program532'] if i in keepC],'typed_program':[[i,d] for i,d in s['shared_full_typed_program757'] if i in keepT],'Gamma_two':s['Gamma_two'],'D_two_by10_by56':s['D_two_by10_by56'],'OLD_Gamma_columns':[[i,gc[i]] for i in old],'OLD_D_columns':[[i,dc[i]] for i in old],'point27_source_indices':piv,'FIRST_column_recipes':recipes,'native_L0_selected_columns':[[i,L0[i]] for i in sorted(u0)],'native_Le_selected_columns':[[i,Le[i]] for i in sorted(ue)],'accepted_native_Gamma_plane_targets':[[i,v92[i]] for i in sorted(needed92)],'ordinary_B21_basis':s['ordinary_B21_basis'],'pair_order10':pairs,'Q_two6':Q,'positive_global_envelopes':p['positive_common_envelopes'],'inherited_committed_prior354':p['future_value_committed_prior'],'inherited_temporary_prior1063':p['future_value_temporary_prior'],'all_nodes_and_callbacks_exact_subset_of_shared605_proof':True,'monotone_positive_pruning_no_new_coefficient_or_sum_term':True,'native20_pruned_formal_aux_actual_kernel_packet_EXTERNAL_SOURCE_pending':True,'rank_source_prior':'Row0/row1 integerclear byupper-envelope D0,D1 givesabsentries<=U0/U1 (166/167bits);2x2det/numerators<=334bits;factor-checkproducts<=501bits/sums<=502bits TEMP only. Allrank/factorcommits<354,temps<1063.','all_physical_higher_honest_extension_and_outerOLD4_scope_same_shared605_proof':True,'fresh_wholeowner_codec_immutableworker_custody_peer_ROOTGO_REQUIRED':True}
 ep=out/'pruned20-exceptional-OLD-FIRST-outer-source-and605-coverage.json';H.dump(ep,data)
 owned=count([s,p,gr,pr,zr,data,C,t,keepC,keepT,boundary,leaves,recipes,gc,dc,v92,pins])+15000;assert owned<150000
 report={'id':STEM,'all_passed':True,'only_exact_source_ancestor_and_callback_pruning':True,'cases20':len(cases),'rawcases21':len(needed),'OLD_columns':len(old),'C0_nodes':len(keepC),'typed_nodes':len(keepT),'boundary_leaves':len(boundary),'FIRST_recipes':len(recipes),'native92_sources':len(needed92),'pruned_input_numeric_count':count(data),'source_owned_upper':owned,'RAW_read_upper':raw,'all20_values_UNRUN':True,'no_new_physical_gauge_implicit_or_jet_coefficients':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','C0_nodes','typed_nodes','boundary_leaves','FIRST_recipes','native92_sources','pruned_input_numeric_count','source_owned_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

#!/usr/bin/env python3
"""ONE554 six-source/OLD/known-FIRST manifest only; no response dot product."""
import argparse,json,resource,signal,time
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER554_OLD_CARRY_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(exist_ok=False,parents=True);pins={};raw=0
 def read(p):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024
  return json.loads((ROOT/p).read_text())
 def field(p,k):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024
  return H.selected(p,k)
 s=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/combined-OLD2-D-and-four-pivot-source-factor-inventory.json')
 a=field('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json','all_zero_bottom_and_Y0_source_combinations')[554]
 wanted=[[289,'25/4'],[379,'-1'],[643,'29/4'],[645,'25'],[651,'-25/2'],[670,'1']]
 assert a['original768_combination']==s['original619_c0_source_combinations'][554]==wanted
 aa=a['actual_a1_27_times9'];assert len(aa)==29
 plane=[[i,q] for i,q in aa if i%9 in [7,8]];assert plane==[[5*9+7,'-175/3'],[5*9+8,'-377/60']]
 kr=next(v for v in s['kernel617_source_recipes'] if v['free_case619']==554);assert kr['subtract_actual_basis563609']==['0','0']
 h=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_INPUT_GATE_20261105_v63/617-direct-old-known-first-formal-Jacobian-positive-input.json');t=dict(h['typed_source_program']);u=dict(s['typed_program']);assert all(t[i]==n for i,n in u.items() if i in t);t.update(u)
 dcols=dict(s['D_old_root_columns']);gcols=dict(h['pruned_original_old_columns_roots']);old=[[j,q] for j,q in wanted if j<640];nc=[[j-640,q] for j,q in wanted if j>=640]
 keep=set();todo=[i for j,q in old for row,i in dcols[j]];cleaf={0}
 for j,q in old:
  for row,(kind,i) in gcols[j]:
   if kind=='typed':todo.append(i)
   else:assert kind=='C0';cleaf.add(i)
 while todo:
  i=todo.pop()
  if i in keep:continue
  keep.add(i);n=t[i]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
  else:assert n[0]=='actual_C0_scalar_function'
 leaf={t[i][1] for i in keep if t[i][0]=='actual_C0_scalar_function'}|(cleaf-{0});assert leaf=={7199,20117,20138,27897,27948,40895,40916}
 cp='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63/pruned-known-Cfull224-input-and-positive-source-envelopes.json';c=read(cp);piv=c['pivot27_source_indices'];assert piv[5]==644 and c['Gamma_two']==h['Gamma_two']
 factors=[];uses=set();pairs=s['pair_order10'];omitted=[]
 for slot,q in aa:
  p,var=divmod(slot,9);hit=False
  for k,(i,j) in enumerate(pairs):
   if var==i:factors.append([p,k,j,q]);uses.add((p,j));hit=True
   if var==j and i!=j:factors.append([p,k,i,q]);uses.add((p,i));hit=True
  if not hit:assert var in [2,3,6];omitted.append([slot,q])
 source_pre=count([s,a,h,t,c,keep,leaf,factors,uses,pins])+60000;assert source_pre<150000
 oldfirst=field(H.PATHS['old640'],'source_columns640_rows56_sparse_dual9');nativefirst=field(H.PATHS['old128'],'source_columns128_rows48_sparse_dual9');recipes=[]
 for p,var in sorted(uses):
  orig=piv[p]
  if var==8:
   kind,arg=c['pivot27_ell_recipes'][p]
   if kind=='literal':r={'kind':'literal56','column':arg}
   else:
    assert kind=='native' and arg==orig-640;r={'kind':'native_ell_SUM','native128_source':arg,'Se':c['native_Se'][arg],'S0':c['native_S0'][arg],'formula':'L0*Se+Le*S0 (SUM)'}
  else:
   dual=8 if var==7 else var+1;cols=oldfirst[orig] if orig<640 else nativefirst[orig-640];v=[]
   for row,values in cols:
    q=dict(values).get(dual,'0')
    if q!='0':v.append([row,q])
   r={'kind':'literal56','column':v,'old8_dual_index':dual,'native_bottom8_exact0_if_native':orig>=640}
  recipes.append({'pivot27':p,'original768':orig,'variable9':var,'first_column_recipe':r})
 del oldfirst,nativefirst
 S0=field('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json','native_S0_columns')
 b=read('reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63/ONE508-source-closure-B21-and-exact-carry-zeros.json')
 assert read('reports/AMBIENT_B18_NATIVE92_ROOT_DIRECT_REVIEW_20261105_v63/report.json')['all_passed']
 ds=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63/actual617-formal-combined-OLD2-second-symbol-coefficients.json');dd=dict(ds['all617_combined_OLD2_independent_symbol_coefficients'])[554]
 assert all(i<299000 for ch in dd for p,i,q in ch)
 ev={'selected_case619':554,'original768_source_all6':wanted,'OLD_source2':old,'native128_source4':nc,'same617_source_recipe':kr,'actual_full_a1_27times9':aa,'only_plane_a1_pivot5_original644':plane,'OLD_full_typed_program':[[i,t[i]] for i in sorted(keep)],'OLD_C0_constant0_and_physical_direct_roots':sorted(cleaf),'OLD_physical_function_leaf_ids7':sorted(leaf),'OLD_D_source_columns':[[j,q,dcols[j]] for j,q in old],'OLD_Gamma_source_columns':[[j,q,gcols[j]] for j,q in old],'Gamma_two':h['Gamma_two'],'D_two_by10_by56':s['D_two_by10_by56'],'source_pair_order10':pairs,'all_known_FIRSTcarry_factor_recipes':factors,'omitted_a1_slots_only_inactive_Dpair_zero':omitted,'all_needed_Xp_FIRST_column_recipes':recipes,'native_L0':c['native_L0'],'native_Le':c['native_Le'],'individual_Xpell_native_SUM_retained':True,'native_weighted_full_S0_factor_inputs':[[j,q,S0[j]] for j,q in nc],'weighted_native_S0_certificate_owned_by_V62_NATIVE_SOURCE_not_claimed_from_individual_entries':True,'native_pivot644_distinct_mask0_FULLL2S0_retained':True,'native92_plane_reference':{'source_report':'reports/AMBIENT_B18_NATIVE_PIVOT23_SECOND_VALUES_20261105_v63/report.json','independent_direct_peer':'reports/AMBIENT_B18_NATIVE92_ROOT_DIRECT_REVIEW_20261105_v63/report.json','original768_source':644,'native23_position':1,'field':'Gamma_Xp2_native23_twochannels_rr_rell','accepted_target_nonzero_support_only_no_weighted_pairing':True},'outer_plane_factor_recipe':'SUBTRACT GammaXp644_rr*(-377/60)+GammaXp644_rell*(-175/3); BOTHcachedtargetsactual/fullL2S0, notzero like648','ordinary_basis_B21':b['ordinary_basis_B21'],'raw33_to_B21_selector21':b['raw33_to_B21_selector21'],'native_base_name':'J18_with_r2ell; NEVERoldA18_m3','complete_upper_recipe':'Gamma(sum6c0X3-Xp2a1)[r²ell]+D(sum6c0X2-Xp1a1), fullPsi/currentpointchart retained','no_dot_native_first_second_third_or_upper_response_values_computed':True}
 ep=out/'ONE554-six-c0-OLD-and-full29-a1-known-FIRST-source.json';H.dump(ep,ev)
 owned=count([s,a,h,t,c,keep,leaf,factors,uses,recipes,S0,b,ds,ev,pins])+15000;assert owned<150000
 report={'id':STEM,'all_passed':True,'SOURCE_only_selected_nonzero_outer_support_candidate':True,'OLD_typed_nodes':len(keep),'OLD_physical_leaves7':True,'all_c0_terms6_and_a1terms29_retained':True,'known_FIRST_factor_records':len(factors),'point_FIRST_recipes':len(recipes),'source_prearith_bound':source_pre,'source_owned_bound':owned,'whole_SOURCE_RAW_upper':raw,'no_signed_response_values':True,'native_weighted_S0_and_whole_native_positive_external_pending':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','OLD_typed_nodes','known_FIRST_factor_records','point_FIRST_recipes','source_owned_bound','whole_SOURCE_RAW_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

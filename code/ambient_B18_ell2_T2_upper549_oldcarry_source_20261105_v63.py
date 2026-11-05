#!/usr/bin/env python3
"""ONE549 OLD/known-FIRST carry source manifest; no response dot products."""
import argparse,json,resource,signal,time
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER549_OLD_CARRY_SOURCE_20261105_v63'
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
 a=field('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json','all_zero_bottom_and_Y0_source_combinations')[549]
 assert a['original768_combination']==s['original619_c0_source_combinations'][549]==[[298,'1/4'],[649,'1']]
 aa=a['actual_a1_27_times9'];assert len(aa)==32
 plane=[[i,q] for i,q in aa if i%9 in [7,8]];assert plane==[[6*9+7,'-79/24'],[6*9+8,'-11/30']]
 kr=next(v for v in s['kernel617_source_recipes'] if v['free_case619']==549);assert kr['subtract_actual_basis563609']==['0','0']
 h=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_INPUT_GATE_20261105_v63/617-direct-old-known-first-formal-Jacobian-positive-input.json');t=dict(h['typed_source_program']);u=dict(s['typed_program']);assert all(t[i]==n for i,n in u.items() if i in t);t.update(u)
 dr=dict(s['D_old_root_columns'])[298];gr=dict(h['pruned_original_old_columns_roots'])[298];keep=set();todo=[i for row,i in dr]
 while todo:
  i=todo.pop()
  if i in keep:continue
  keep.add(i);n=t[i]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
  else:assert n[0]=='actual_C0_scalar_function'
 leaf={n[1] for i,n in t.items() if i in keep and n[0]=='actual_C0_scalar_function'};assert leaf=={7454,20132,40910}
 assert gr==[[1,['C0',7454]],[55,['C0',0]]]
 cp='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63/pruned-known-Cfull224-input-and-positive-source-envelopes.json';c=read(cp);piv=c['pivot27_source_indices'];assert piv[6]==648
 factors=[];uses=set();pairs=s['pair_order10'];omitted=[]
 for slot,q in aa:
  p,var=divmod(slot,9);hit=False
  for k,(i,j) in enumerate(pairs):
   if var==i:factors.append([p,k,j,q]);uses.add((p,j));hit=True
   if var==j and i!=j:factors.append([p,k,i,q]);uses.add((p,i));hit=True
  if not hit:omitted.append([slot,q]);assert var==6
 old=field(H.PATHS['old640'],'source_columns640_rows56_sparse_dual9');native=field(H.PATHS['old128'],'source_columns128_rows48_sparse_dual9')
 recipes=[]
 for p,var in sorted(uses):
  orig=piv[p]
  if var==8:
   kind,arg=c['pivot27_ell_recipes'][p]
   if kind=='literal':rec={'kind':'literal56','column':arg}
   else:
    assert kind=='native' and arg==orig-640
    rec={'kind':'native_ell_SUM','native128_source':arg,'Se':c['native_Se'][arg],'S0':c['native_S0'][arg],'formula':'L0*Se+Le*S0 (SUM)'}
  else:
   dual=8 if var==7 else var+1;cols=old[orig] if orig<640 else native[orig-640];v=[]
   for row,values in cols:
    q=dict(values).get(dual,'0')
    if q!='0':v.append([row,q])
   rec={'kind':'literal56','column':v,'old8_dual_index':dual,'native_bottom8_exact0_if_native':orig>=640}
  recipes.append({'pivot27':p,'original768':orig,'variable9':var,'first_column_recipe':rec})
 del old,native
 # The full point column must vanish AFTER trace lowering/input mask1.
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';S0=field(mug,'native_S0_columns');assert S0[9]==[]
 lower=read('reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63/ONE508-source-closure-B21-and-exact-carry-zeros.json')
 ev={'selected_case619':549,'original768_source':a['original768_combination'],'same617_source_recipe':kr,'actual_full_a1_27times9':aa,'only_plane_a1_pivot6_source648':plane,'OLD_component_source298_scale':'1/4','OLD_full_typed_program':[[i,t[i]] for i in sorted(keep)],'OLD_physical_function_leaf_ids3':sorted(leaf),'OLD_D_source_roots':dr,'OLD_Gamma_source_roots':gr,'Gamma_two':h['Gamma_two'],'D_two_by10_by56':s['D_two_by10_by56'],'source_pair_order10':pairs,'all_known_FIRSTcarry_factor_recipes':factors,'firstcarry_omitted_slots_only_inactive_z6_D_support_zero':omitted,'all_needed_Xp_first_column_recipes':recipes,'native_L0':c['native_L0'],'native_Le':c['native_Le'],'individual_Xpell_native_SUM_retained':True,'ordinary_basis_B21':lower['ordinary_basis_B21'],'raw33_to_B21_selector21':lower['raw33_to_B21_selector21'],'native_source128_j9':{'original768':649,'S32_column':2,'input_trace_mask':1,'full_intrinsic_S0_column48':S0[9]},'native_pivot648_distinct_inputmask0_S0_not_set0':True,'native92_plane_reference':{'source_report':'reports/AMBIENT_B18_NATIVE_PIVOT23_SECOND_VALUES_20261105_v63/report.json','original768_source':648,'native23_position':2,'field':'Gamma_Xp2_native23_twochannels_rr_rell','reference_only_until_independent_parent_adoption':True,'no_values_read_or_contracted_here':True},'GammaXp2_a1_plane_factor_recipe':'GammaXp648_rr*(-11/30)+GammaXp648_rell*(-79/24), SUBTRACTfromfullupper; keepFULLL2S0','complete_upper_recipe':'Gamma(OLD2983/4+native6493 -Xp2*a1)[r²ell] +D(OLD2982/4+native6492-Xp1*a1). FullPsi/currentM0 alreadyinGamma,D.','native_J18_with_r2ell_source_owned_by_V62':True,'no_knownFIRST_dot_second_third_or_upper_response_values_computed':True}
 ep=out/'ONE549-OLD298-and-full32-a1-known-FIRST-carry-factor-source.json';H.dump(ep,ev)
 owned=count([s,a,aa,h,t,keep,leaf,c,factors,uses,recipes,S0,lower,ev,pins])+15000;assert owned<150000
 report={'id':STEM,'all_passed':True,'SOURCE_metadata_and_existing_first_literal_selection_only':True,'OLD_typed_nodes':len(keep),'known_FIRST_factor_records':len(factors),'unique_point_FIRST_column_recipes':len(recipes),'no_inactive_first_dropped_without_zeroDpair_proof':True,'source_owned_bound':owned,'whole_SOURCE_RAW_upper':raw,'no_new_response_values':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','OLD_typed_nodes','known_FIRST_factor_records','unique_point_FIRST_column_recipes','source_owned_bound','whole_SOURCE_RAW_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

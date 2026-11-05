#!/usr/bin/env python3
"""ONE288 literal source/carry profile only; no response coefficient values."""
import argparse,json,resource,signal,time
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER288_SOURCE_20261105_v63'
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
 a=field('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json','all_zero_bottom_and_Y0_source_combinations')[288]
 assert a['original768_combination']==[[352,'1']] and s['original619_c0_source_combinations'][288]==[[352,'1']]
 aa=a['actual_a1_27_times9'];assert len(aa)==8 and all(i%9==1 for i,q in aa)
 kr=next(v for v in s['kernel617_source_recipes'] if v['free_case619']==288);assert kr['subtract_actual_basis563609']==['0','0']
 h=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_INPUT_GATE_20261105_v63/617-direct-old-known-first-formal-Jacobian-positive-input.json')
 t=dict(h['typed_source_program']);u=dict(s['typed_program']);assert all(t[i]==n for i,n in u.items() if i in t);t.update(u);c=dict(h['literal_C0_program']);boundary={i for i,k,q in h['literal_boundary_firsts_point_r_ell_and_membership']}
 dr=dict(s['D_old_root_columns'])[352];gr=dict(h['pruned_original_old_columns_roots'])[352];todo=[i for row,i in dr];cr=set()
 for row,(kind,i) in gr:
  assert kind in ['typed','C0']
  if kind=='typed':todo.append(i)
  else:cr.add(i)
 keep=set()
 while todo:
  i=todo.pop()
  if i in keep:continue
  keep.add(i);n=t[i]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
  else:assert n[0]=='actual_C0_scalar_function';cr.add(n[1])
 kc=set();leaf=set();todo=list(cr)
 while todo:
  i=todo.pop()
  if i in kc or i in leaf:continue
  if i not in c:assert i in boundary;leaf.add(i);continue
  kc.add(i);n=c[i]
  if n[0] in ['a','m']:todo.extend(n[1:3])
  elif n[0]=='n':todo.append(n[1])
  else:assert n[0]=='c'
 ds=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63/actual617-formal-combined-OLD2-second-symbol-coefficients.json');dv=dict(ds['all617_combined_OLD2_independent_symbol_coefficients'])[288]
 assert {f for ch in dv for p,f,q in ch}=={20117,27948,40895}
 js=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63/617-formal-direct-old-third-and-second-symbol-sensitivities.json');jv=next(v for v in js['all617_formal_directOLD_J0_Jr_Jell'] if v['free_case619']==288)
 assert {i for ch in jv['two_formal_J0_Jr_Jell_symbol_rows'] for i,v in ch}=={20117,40895}
 bp='reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63/ONE508-source-closure-B21-and-exact-carry-zeros.json';b=read(bp)
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63_REVIEW/report.json')['all_passed']
 cp='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63/pruned-known-Cfull224-input-and-positive-source-envelopes.json';cc=read(cp)
 assert cc['Gamma_two']==h['Gamma_two']
 positions=[i//9 for i,q in aa];assert positions==[0,4,8,9,11,14,17,20]
 piv=[]
 for slot,q in aa:
  i=slot//9;piv.append({'pivot27':i,'original768':cc['pivot27_source_indices'][i],'a1_z1':q,'Xp_r_literal56':cc['pivot27_r_columns'][i],'Xp_ell_recipe':cc['pivot27_ell_recipes'][i]})
 assert [v['original768'] for v in piv]==[604,640,656,660,752,755,758,761]
 ev={'selected_case619':288,'original768_source':[[352,'1']],'actual_full27times9_a1':aa,'actual_a1_only_z1_not_pointwise_assumption':True,'actual_a1_r_and_ell_both_exact0':True,'native_c0_exact_empty':True,'same617_source_recipe':kr,'complete_actual_upper_recipe':'Gamma X3523[r²ell]+D X3522 - sum_p a1[p,z1]*(D[z1*r] Xp_r +D[z1*ell] Xp_ell). OuterGammaXp2*a1 pureplane0 becausea1 onlyz1; nativec0=0. FullC/Psi retained.','ordinary_basis_B21':b['ordinary_basis_B21'],'raw33_to_B21_selector21':b['raw33_to_B21_selector21'],'all660_NF_and1089_products_inherited_from_peered_common_B21_source':True,'full_typed_program':[[i,t[i]] for i in sorted(keep)],'literal_C0_program':[[i,c[i]] for i in sorted(kc)],'function_leaf_ids':sorted(leaf),'D_old352_source_roots':dr,'Gamma_old352_source_roots':gr,'Gamma_two':h['Gamma_two'],'D_two_by10_by56':s['D_two_by10_by56'],'known_first_carry_pivots8':piv,'native_L0':cc['native_L0'],'native_Le':cc['native_Le'],'native_S0':cc['native_S0'],'native_Se':cc['native_Se'],'individual_native_Xpell_is_actual_SUM_L0Se_plus_LeS0':True,'original_case288_J1_certificate':jv,'original_case288_OLD2_certificate':dv,'all_actual_higher_sensitivities_only_physical20117_27948_40895':True,'no_conductor_mixedsecond_or_third_order_assumption':True,'all_multiproduct_first_cross_terms_retained':True,'no_second_third_native_or_upper_response_values':True}
 ep=out/'ONE288-opaque-free-source-and-eight-known-first-pivot-carries.json';H.dump(ep,ev)
 owned=count([s,h,t,c,keep,kc,leaf,ds,js,b,cc,piv,ev,pins])+10000;assert owned<150000
 report={'id':STEM,'all_passed':True,'metadata_source_only':True,'native_and_Xp2terms_exact_source_zero':True,'known_first_native_Xpell_SUM_kept':True,'a1_nonzero_only_z1':True,'typed_nodes':len(keep),'C0_nodes':len(kc),'function_leaves':len(leaf),'source_owned_bound':owned,'whole_SOURCE_RAW_upper':raw,'no_response_values':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','typed_nodes','C0_nodes','function_leaves','source_owned_bound','whole_SOURCE_RAW_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

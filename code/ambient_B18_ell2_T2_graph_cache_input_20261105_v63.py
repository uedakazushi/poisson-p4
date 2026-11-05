#!/usr/bin/env python3
"""Final invariant/raw33 + A23 sensitivity/cache input gate, no value job."""
import argparse,json,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,srec
from ambient_B18_native128_first_20261104_v63 import PATHS as OLD_PATHS
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_GRAPH_CACHE_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 factor=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63/actual3960-rho-factor-identities-or-obstructions.json');assert factor['combined_two_rho_factor_raw33']
 raw=read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json')
 known=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json');mini=read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json')
 fold=read('reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63/actual-two-rell-source-folds-and-retained-native-first-operators.json');chain=read(H.PATHS['chain']);ops=chain['operator_scalar_nodes'];cols=list(map(dict,chain['B48_by640_sparse_columns']))
 ad=read('reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_AD_20261105_v63/actual-two-by619-old-formal-rell-second-sensitivities.json');grads={i:H.parse(g) for i,g in ad['typed_node_original_second_gradients']}
 atom=read('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json');a0={i:H.parse(p).get(0,F(0)) for i,p in atom['all79_atom_values']}
 catalog=read('reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json');cg={i:H.parse(g) for i,v,g in catalog['all108_composite_values_and_opaque_second_gradients']}
 actualold=read('reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63/actual104-source-rell-values-and-two-by619-old-knownparts.json');oldlabels={i for i,v,g in actualold['all104_original_C0_values_and_opaque_gradients']}
 hs[H.PATHS['C0']]=H.sha(ROOT/H.PATHS['C0']);P=H.selected(H.PATHS['C0'],'honest_D0_relation_columns');Cbasis=H.selected(H.PATHS['C0'],'actual_C0_basis8_scalar_nodes')
 newroots=set();allrecords=[];direct=set();typed=set()
 for basis,j in enumerate(known['selected_bottom_columns']):
  rawc=mini['all642_zero_constant_combinations'][j]['original768_combination'];c=H.parse(rawc);fs=[]
  for gamma in map(H.parse,fold['two_Gamma_on56_actual_X2']):
   f={}
   for source,a in c.items():
    if source>=640:continue
    for row,b in gamma.items():
     if row<48:
      root=cols[source].get(row)
      if root is None:continue
      n=ops[root]
      if n[0]=='c':continue
      if n[0]=='actual_C0_scalar_function':f=add(f,{n[1]:a*b})
      else:
       typed.add(root)
       if root not in grads:
        newroots.add(root);assert n[0]=='product' and n[1:3]==[-1,1] and len(n[3])==2
        x,y=n[3];assert ops[x][0]==ops[y][0]=='actual_C0_scalar_function';ix,iy=ops[x][1],ops[y][1];assert a0[ix]==a0[iy]==0;grads[root]={}
       f=add(f,grads[root],a*b)
     else:f=add(f,{P[source][row-48+1]:a*b})
   fs.append(f);direct.update(f)
   assert all(i==0 or i in cg for i in f)
   assert all(not cg.get(i,{}) for i in f)
  allrecords.append({'basis_column':basis,'source642':j,'original768':rawc,'two_literal_old_second_functionals_after_typed_AD':[srec(v) for v in fs],'all_honest_second_sensitivities_zero':True})
 assert newroots=={516,785,846};assert direct-oldlabels<=set(Cbasis)
 source=read(OLD_PATHS['source']);rawsrc=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json');roots=dict(rawsrc['all_requested164_root_values']);roots.update(zip(rawsrc['joint105_unknown_order'],rawsrc['all105_unknown_final_values']));nuids=[[source['families'][f]['p_quadratic_relation'][k] for k in ['s','t']] for f in ['CA','CB']];assert nuids==[[11,12],[20537,20538]] and all(i in roots for pair in nuids for i in pair)
 physical=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json');r5=physical['physical12_coefficients'][5];r7=physical['physical12_coefficients'][7];assert r5['physical_row']==67 and r7['physical_row']==63
 native=read('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/report.json');assert native['all_passed'];read('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63_REVIEW/report.json')
 # Existing-cache closures only; coefficients remain referenced, not re-evaluated.
 refs={
  'raw33':str(Path('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json')),
  'physical12':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
  'source105':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',
  'actual_old_C0':'reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63/actual104-source-rell-values-and-two-by619-old-knownparts.json',
  'typed_AD':'reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_AD_20261105_v63/actual-two-by619-old-formal-rell-second-sensitivities.json',
  'native':'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json',
  'Mugate':'reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',
  'old640':H.PATHS['old640'],'old128':H.PATHS['old128'],
  'old640ell':'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json',
  'mini':'reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json'}
 for p in refs.values():hs[p]=H.sha(ROOT/p)
 slots={'raw33_actions_and_multiplication_upper':9*33*33+33*33,'full105finite_source_input_upper':105*33,
  'scalarPhi_input_records_and_bar_rho':10*33+2*3*33,'oldfirst_sparse_and_ell_inputs_upper':4271+7928,
  'native_S_and_L_and_A23_scratch_upper':8000,'selected_cache_fields_and_metadata_upper':12000,
  'output2plus46plus46upper':100}
 total=sum(slots.values());assert total<60000
 ep=out/'actual-combined-invariant-and-A23-cache-source-input.json';H.dump(ep,{'actual_two_three_component_bar_rho':factor['actual_two_three_component_bar_rho33'],
  'raw33_primitive_source_formula':{'Wv':[[23,'-3375/11'],[29,'-50'],[30,'25/2']],
   'Phi0':'(3375/11)*nuB*bB67-50*delta63','Phi1':'-(3375/11)*nuA*bA67+(25/2)*delta63','Phi2':'-(3375/11)*delta67+50*bB63+(25/2)*bA63',
   'physical5_row':67,'physical7_row':63,'actual_monic_trace_s_t_nodes':nuids,'nu_formula':'t+s_trace²/4; trace_s isoriginalsourcefunction, notF9sA','use_full33_nu_with_z_firsts':True},
  'all23_old_source_sensitivity_proofs':allrecords,'three_new_typed_root_zero_second_sensitivity':[[i,ops[i],[[j,ops[j],str(a0[ops[j][1]])] for j in ops[i][3]]] for i in sorted(newroots)],
  'all23_unknown_honest_second_sensitivities_zero':True,
  'new_direct_C0_labels_are_literal_coordinates':sorted(direct-oldlabels),
  'A23_actual_values_not_computed':True,'invariant_dot_values_not_computed':True,
  'A23_native_scope':'eachoriginal23combinationS0cN0, useactualL0S2+L1rOmegaEll+L1ellOmegaR andactualfirstpivotcorrections; noindividualeffectiveXp shortcut',
  'normalized_graph_scope':'canonicalfourmu remainunspecified; coherentunitgauge givesmuprime=(alpha,0,beta,0), alpha=mu0r,beta=mu0ell-mu1r, so normalizedg0=A0+ce*alpha+cr*beta, normalizedg1=A1. Noallorderhigherrelationclaim.',
  'existing_cache_paths':refs,'whole_cache_slot_bound':total,'slot_breakdown':slots,
  'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'QQ_coefficients':60000,'bits':2048},
  'no_new_physical_implicit_gauge_or_jet_coefficients':True})
 report={'id':STEM,'all_passed':True,'all23_honest_second_sensitivities_zero':True,'raw33_two_invariant_dot_inputs_complete':True,
  'actual_A23_or_invariant_values_computed':False,'whole_cache_slot_bound':total,'input_sha256':hs,
  'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','all23_honest_second_sensitivities_zero','whole_cache_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

#!/usr/bin/env python3
"""Positive denominator/absolute-numerator bounds only; no Phi/A sums."""
import argparse,json,math,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_GRAPH_CACHE_HEIGHT_BOUND_20261105_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);hs={};rows=[]
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def emit(name,D,U,kind,parents=(),terms=1):
  assert D>0 and U>=0;rows.append({'name':name,'common_denominator':str(D),'abs_cleared_numerator_upper':str(U),'denominator_bits':D.bit_length(),'numerator_bits':U.bit_length(),'kind':kind,'parents':list(parents),'terms':terms});return(D,U,name)
 def inputs(name,vs):
  vs=list(map(F,vs));D=1
  for v in vs:D=math.lcm(D,v.denominator)
  U=max([abs(v.numerator)*(D//v.denominator) for v in vs]+[0]);return emit(name,D,U,'known_input_LCM_max_abs_cleared_numerator',terms=len(vs))
 def prod(name,*es,terms=1):return emit(name,math.prod(e[0] for e in es),terms*math.prod(e[1] for e in es),'positive_product_sum_bound',[e[2] for e in es],terms)
 def union(name,*es):
  D=math.lcm(*(e[0] for e in es));return emit(name,D,max(e[1]*(D//e[0]) for e in es),'union_max',[e[2] for e in es])
 def summ(name,*es):
  D=math.lcm(*(e[0] for e in es));return emit(name,D,sum(e[1]*(D//e[0]) for e in es),'positive_common_denominator_sum',[e[2] for e in es])
 gate=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_INPUT_GATE_20261105_v63/actual-combined-invariant-and-A23-cache-source-input.json');refs=gate['existing_cache_paths']
 NF=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json')
 table=inputs('raw33_product_coefficient_upper',[c for v in NF['all660_actual_raw33_normal_forms'] for i,c in v]);src=read(refs['source105']);root=dict(src['all_requested164_root_values']);root.update(zip(src['joint105_unknown_order'],src['all105_unknown_final_values']))
 traces=gate['raw33_primitive_source_formula']['actual_monic_trace_s_t_nodes'];s=inputs('trace_s_source33',[c for pair in traces for i,c in root[pair[0]]]);t=inputs('trace_t_source33',[c for pair in traces for i,c in root[pair[1]]]);quarter=inputs('quarter',[F(1,4)]);one=inputs('one',[1])
 # ν=t+s²/4: each output hasatmost33² productsandone tablecoefficient.
 nu=summ('nu_envelope',t,prod('s_squared_over4',s,s,quarter,table,terms=33*33))
 physical=read(refs['physical12']);r5=physical['physical12_coefficients'][5];r7=physical['physical12_coefficients'][7]
 b=inputs('all_primitive_b_source33',[c for r in [r5,r7] for key in ['bA','bB'] for i,c in r[key]]);delta=inputs('primitive_delta_source33',[c for r in [r5,r7] for i,c in r['delta']]);cw=inputs('Wv_constants',[F(-3375,11),F(-50),F(25,2)])
 nub=prod('nu_times_b',nu,b,table,terms=33*33);phi=prod('Phi3vector_raw33',cw,union('primitive_product_or_delta',nub,delta,b),terms=3)
 bar=inputs('both_combined_bar_rho',[c for inv in gate['actual_two_three_component_bar_rho'] for component in inv for i,c in component]);invariants=prod('both_invariants',bar,phi,terms=3*33)
 ad=read(refs['typed_AD']);q=inputs('typed_known_rell',[v[3] for i,v in ad['typed_node_known4_coefficients']]);grad=inputs('typed_second_gradients',[c for i,v in ad['typed_node_original_second_gradients'] for j,c in v]);cv=read(refs['actual_old_C0']);oldvalues=inputs('known_original_C0_rell',[v[3] for i,v,g in cv['all104_original_C0_values_and_opaque_gradients']])
 first=read('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json');firstR=inputs('all79_atom_r_firsts',[dict(v).get(8,'0') for i,v in first['all79_atom_values']])
 # New3 products'ell firsts lie in copied original79 first-source data. Bound
 # original restoration84 arithmetic withknown input derivative denominators
 # by their already completed actualell scalar outputs from the µ compiler.
 mu=read(refs['Mugate']);leaf={r['original_node']:H.parse(r['origin_and_ds']) for r in mu['all_exact_cached_leaf_bindings_origin_and_ds']}
 # Only first-jet arithmetic here: these functionfirsts were already used by
 # completedµjob; no second/higher/Φ/A coefficient is evaluated.
 for i,n in mu['old_exact84_original_arithmetic_nodes']:
  if n[0]=='c':v={0:F(n[1],n[2])} if n[1] else {}
  elif n[0]=='a':v=H.parse([[j,str(c)] for j,c in __import__('ambient_B18_ell2_raw33_quotient_20261105_v63').add(leaf[n[1]],leaf[n[2]]).items()])
  elif n[0]=='n':v={j:-c for j,c in leaf[n[1]].items()}
  else:v=H.pairmul(leaf[n[1]],leaf[n[2]])
  leaf[i]=v
 firstE=inputs('all79_atom_ell_firsts',[leaf[i].get(1,F(0))/5 for i in mu['old79_atom_ids']]);newproduct=prod('new3_zero_origin_products',firstR,firstE,terms=2)
 mini=read(refs['mini']);comb=inputs('selected23_source_coefficients',[c for r in gate['all23_old_source_sensitivity_proofs'] for i,c in r['original768']]);Gamma=inputs('both_Gamma',[c for v in read('reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63/actual-two-rell-source-folds-and-retained-native-first-operators.json')['two_Gamma_on56_actual_X2'] for i,c in v]);Nc=max(len(r['original768']) for r in gate['all23_old_source_sensitivity_proofs'])
 old_direct=prod('old_direct_source_sum',comb,Gamma,oldvalues,terms=Nc*9);old_typed=prod('old_typed_known_source_sum',comb,Gamma,union('q_or_newproduct',q,newproduct),terms=Nc*9)
 old_grad=prod('old_typed_consolidated_known_leaf_sum',comb,Gamma,grad,oldvalues,terms=Nc*9*104)
 oldpart=summ('old_A23_part',old_direct,old_typed,old_grad)
 native=read(refs['native']);S0=inputs('S0',[c for v in mu['native_S0_columns'] for i,c in v]);Se=inputs('S1ell',[c for v in mu['native_S1ell_columns'] for i,c in v]);S2=inputs('S2rell',[c for v in native['intrinsic_S2rell_columns'] for i,c in v]);L0=inputs('L0',[c for v in mu['L0_Res4_Ntranspose_columns'] for i,c in v]);Lr=inputs('L1r',[c for v in native['actual_L1r_columns'] for i,c in v]);Le=inputs('L1ell',[c for v in native['actual_L1ell_columns'] for i,c in v])
 gs=read('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json')['actualG8'];Sr=inputs('S1r',[c for v in gs for c in v[4:8]])
 nx2=summ('native_X2_one_column',prod('L0S2',L0,S2,terms=48),prod('LrSe',Lr,Se,terms=48),prod('LeSr',Le,Sr,terms=48))
 nativepart=prod('native_A23_part',Gamma,comb,nx2,terms=48*Nc)
 oxr=H.selected(H.PATHS['old640'],'source_columns640_rows56_sparse_dual9')+H.selected(H.PATHS['old128'],'source_columns128_rows48_sparse_dual9');oldr=inputs('existing_all_source_rfirsts',[dict(p).get(8,'0') for col in oxr for row,p in col]);olde=inputs('existing_old640_ellfirsts',[c for col in read(refs['old640ell'])['actual_first640_ell_columns'] for i,c in col]);Ne=summ('actual_individual_native_ellfirst',prod('L0Se',L0,Se,terms=48),prod('LeS0',Le,S0,terms=48));X1=union('all_source_firsts_r_ell',oldr,olde,Ne)
 inv=inputs('fixed_B30_inverse',[c for col in mini['fixed_bottom30_inverse_columns'] for i,c in col]);a1=prod('all_a1_r_ell',inv,comb,X1,terms=30*Nc);corr=prod('all_pivot_first_corrections',Gamma,X1,a1,terms=2*27*56)
 A23=summ('all_actual_A23',oldpart,nativepart,corr)
 coefficients=inputs('normalized_graph_cr_ce',[r['four_unknown_mu_coefficients_per_channel'][k] for r in read('reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json')['all23_linear_graph_source_records'] for k in ['M_r_H0v','M_ell_H0v']]);graph=summ('all_normalized_graph',A23,prod('alpha_ce',invariants,coefficients),prod('beta_cr',invariants,coefficients))
 maximum=max(max(r['denominator_bits'],r['numerator_bits']) for r in rows);assert maximum<=2048,('source_derived_height_guard',maximum)
 report={'id':STEM,'all_passed':True,'maximum_derived_intermediate_or_output_height_bits':maximum,'authorized_cap_bits':2048,
  'proof':'For eachfield allcoefficients havecommonD andabs(D*c)<=U. Products multiplyD/U andtermcounts; sums useLCM andpositiveclearednumeratorbound. No signedPhi,A23,invariantorgraphsumwas evaluated. Provenunknownsenszero permitsdroppingonlyidenticallyzero unknown terms.',
  'bounds':rows,'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'new_second_higher_physical_or_value_job_started':False,
  'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','maximum_derived_intermediate_or_output_height_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

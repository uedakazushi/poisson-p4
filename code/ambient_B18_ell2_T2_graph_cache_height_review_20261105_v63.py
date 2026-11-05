#!/usr/bin/env python3
"""Independent LCM/positive-envelope postcheck; no signed Phi/A/invariant values."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,math,time
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,add,scale
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_GRAPH_CACHE_HEIGHT_BOUND_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','cache':'reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_INPUT_GATE_20261105_v63/actual-combined-invariant-and-A23-cache-source-input.json','cache_peer':'reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_INPUT_GATE_20261105_v63_REVIEW/report.json',
'NF':'reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json','atom':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json','fold':'reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63/actual-two-rell-source-folds-and-retained-native-first-operators.json','nativeN':'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json','graph':'reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def bound(values):
 v=list(map(Q,values));D=math.lcm(*(q.denominator for q in v));U=max([abs(q.numerator)*(D//q.denominator) for q in v]+[0]);return D,U,len(v)
def prod(*es,terms=1):return math.prod(e[0] for e in es),terms*math.prod(e[1] for e in es)
def summ(*es):
 D=math.lcm(*(e[0] for e in es));return D,sum(e[1]*(D//e[0]) for e in es)
def union(*es):
 D=math.lcm(*(e[0] for e in es));return D,max(e[1]*(D//e[0]) for e in es)
def height(e):return max(e[0].bit_length(),e[1].bit_length())
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};R=D['primary'];assert R['all_passed'] and D['cache_peer']['all_passed']
 for p,h in R['input_sha256'].items():assert sha(ROOT/p)==h
 C=D['cache'];refs=C['existing_cache_paths'];E={}
 for k in ['source105','physical12','typed_AD','actual_old_C0','Mugate','mini','native','old640ell']:
  p=refs[k];inputs[k]={'path':p,'sha256':sha(ROOT/p)};E[k]=json.loads((ROOT/p).read_text())
 qsrc=E['source105'];roots=dict(qsrc['all_requested164_root_values']);roots.update(zip(qsrc['joint105_unknown_order'],qsrc['all105_unknown_final_values']));traces=C['raw33_primitive_source_formula']['actual_monic_trace_s_t_nodes'];p5,p7=E['physical12']['physical12_coefficients'][5],E['physical12']['physical12_coefficients'][7];AD=E['typed_AD'];mu=E['Mugate'];native=E['native']
 leaf={r['original_node']:sp(r['origin_and_ds']) for r in mu['all_exact_cached_leaf_bindings_origin_and_ds']}
 for i,n in mu['old_exact84_original_arithmetic_nodes']:
  if n[0]=='c':v={0:Q(n[1],n[2])} if n[1] else {}
  elif n[0]=='a':v=add(leaf[n[1]],leaf[n[2]])
  elif n[0]=='n':v=scale(leaf[n[1]],-1)
  else:
   x,y=leaf[n[1]],leaf[n[2]];v={0:x.get(0,Q(0))*y.get(0,Q(0)),1:x.get(0,Q(0))*y.get(1,Q(0))+x.get(1,Q(0))*y.get(0,Q(0))};v={j:c for j,c in v.items() if c}
  leaf[i]=v
 oldr=[]
 for key,field in [('old640','source_columns640_rows56_sparse_dual9'),('old128','source_columns128_rows48_sparse_dual9')]:
  p=refs[key];inputs[key]={'path':p,'sha256':sha(ROOT/p)}
  oldr.extend(sp(v).get(8,Q(0)) for col in fields(ROOT/p,[field])[field] for row,v in col)
 values={
 'raw33_product_coefficient_upper':[c for p in D['NF']['all660_actual_raw33_normal_forms'] for i,c in p],
 'trace_s_source33':[c for pair in traces for i,c in roots[pair[0]]],
 'trace_t_source33':[c for pair in traces for i,c in roots[pair[1]]],
 'quarter':[Q(1,4)],'one':[Q(1)],
 'all_primitive_b_source33':[c for p in [p5,p7] for k in ['bA','bB'] for i,c in p[k]],
 'primitive_delta_source33':[c for p in [p5,p7] for i,c in p['delta']],
 'Wv_constants':[Q(-3375,11),Q(-50),Q(25,2)],
 'both_combined_bar_rho':[c for q in C['actual_two_three_component_bar_rho'] for p in q for i,c in p],
 'typed_known_rell':[v[3] for i,v in AD['typed_node_known4_coefficients']],
 'typed_second_gradients':[c for i,p in AD['typed_node_original_second_gradients'] for j,c in p],
 'known_original_C0_rell':[v[3] for i,v,g in E['actual_old_C0']['all104_original_C0_values_and_opaque_gradients']],
 'all79_atom_r_firsts':[sp(v).get(8,Q(0)) for i,v in D['atom']['all79_atom_values']],
 'all79_atom_ell_firsts':[leaf[i].get(1,Q(0))/5 for i in mu['old79_atom_ids']],
 'selected23_source_coefficients':[c for p in C['all23_old_source_sensitivity_proofs'] for i,c in p['original768']],
 'both_Gamma':[c for p in D['fold']['two_Gamma_on56_actual_X2'] for i,c in p],
 'S0':[c for p in mu['native_S0_columns'] for i,c in p],
 'S1ell':[c for p in mu['native_S1ell_columns'] for i,c in p],
 'S2rell':[c for p in native['intrinsic_S2rell_columns'] for i,c in p],
 'L0':[c for p in mu['L0_Res4_Ntranspose_columns'] for i,c in p],
 'L1r':[c for p in native['actual_L1r_columns'] for i,c in p],
 'L1ell':[c for p in native['actual_L1ell_columns'] for i,c in p],
 'S1r':[c for p in D['nativeN']['actualG8'] for c in p[4:8]],
 'existing_all_source_rfirsts':oldr,
 'existing_old640_ellfirsts':[c for p in E['old640ell']['actual_first640_ell_columns'] for i,c in p],
 'fixed_B30_inverse':[c for p in E['mini']['fixed_bottom30_inverse_columns'] for i,c in p],
 'normalized_graph_cr_ce':[p['four_unknown_mu_coefficients_per_channel'][k] for p in D['graph']['all23_linear_graph_source_records'] for k in ['M_r_H0v','M_ell_H0v']]}
 envelopes={};records=[]
 for r in R['bounds']:
  name=r['name'];parents=[envelopes[k] for k in r['parents']];kind=r['kind'];terms=r['terms']
  if kind=='known_input_LCM_max_abs_cleared_numerator':
   d,u,n=bound(values[name]);assert n==terms;e=d,u
  elif kind=='positive_product_sum_bound':e=prod(*parents,terms=terms)
  elif kind=='positive_common_denominator_sum':e=summ(*parents)
  else:assert kind=='union_max';e=union(*parents)
  assert e==(int(r['common_denominator']),int(r['abs_cleared_numerator_upper'])) and e[0].bit_length()==r['denominator_bits'] and e[1].bit_length()==r['numerator_bits'];envelopes[name]=e;records.append((name,height(e)))
 assert max(h for n,h in records)==538==R['maximum_derived_intermediate_or_output_height_bits']
 Nc=max(len(r['original768']) for r in C['all23_old_source_sensitivity_proofs']);assert Nc==27
 # Source term-count safety, independently checked rather than inferred from
 # the final inequalities. Sign normalization never increases these bounds.
 required={'s_squared_over4':1089,'nu_times_b':1089,'Phi3vector_raw33':3,'both_invariants':99,'new3_zero_origin_products':2,'old_direct_source_sum':Nc*9,'old_typed_known_source_sum':Nc*9,'old_typed_consolidated_known_leaf_sum':Nc*9*104,'L0S2':48,'LrSe':48,'LeSr':48,'native_A23_part':48*Nc,'L0Se':48,'LeS0':48,'all_a1_r_ell':30*Nc,'all_pivot_first_corrections':2*27*56}
 for r in R['bounds']:
  if r['name'] in required:assert r['terms']>=required[r['name']]
 # Operationalfactoring may materialize Omega and the S1pivot correction;
 # add conservative positive envelopes for these temporary vectors too.
 B=envelopes;ar=B['all_a1_r_ell'];cn=B['selected23_source_coefficients']
 omegaR=summ(prod(B['S1r'],cn,terms=Nc),prod(B['S0'],ar,terms=27));omegaE=summ(prod(B['S1ell'],cn,terms=Nc),prod(B['S0'],ar,terms=27));corrS=summ(prod(B['S1r'],ar,terms=27),prod(B['S1ell'],ar,terms=27))
 extra=[omegaR,omegaE,corrS,prod(B['both_Gamma'],B['L1r'],omegaE,terms=48*56),prod(B['both_Gamma'],B['L1ell'],omegaR,terms=48*56),prod(B['both_Gamma'],B['L0'],corrS,terms=48*56)]
 maximum=max([538]+list(map(height,extra)));assert maximum<=2048
 ep=o/'independent-source-LCM-positive-bounds-and-operational-scratch.json';dump(o,ep.name,{'inputs':inputs,'all_input_class_denominator_LCMs_and_max_abs_cleared_numerators_verified':True,'all_source_derived_positive_envelopes_replayed_without_signed_values':True,'producer_maximum_derived_envelope_bits':538,'independent_operational_scratch_bits':list(map(height,extra)),'independent_conservative_height_bits':maximum,'authorized_cap_bits':2048,
'semantic_coverage':'Raw33tablecoefficientsarecontainedinNFofmonomials<=4; higherbasisproductszerobym5. All23nativeS0cNzero permitsdroppingL2ONLY. ActualindividualnativeellfirstsincludeL0Se+LeS0. B30inverseactsonthe30literalcomplementrows; a1 usesatmost30*Nc products. Pivotfirstcorrectionscoverbothr/ellandall27pivots. Newlytypedproductsuseexactzerooriginandknownfirstcrosses. All sumsareboundedbysumofabsoluteintegerclearednumerators, so partialaccumulationsandallreduceddenominatorsarecovered.',
'no_signed_Phi_A23_invariant_or_graph_value_evaluated':True})
 report={'id':STEM,'all_passed':True,'prior_source_derived_height_bound_accepted':True,'maximum_conservative_height_bits':maximum,'authorized_cap_bits':2048,'actual_invariant_A23_or_graph_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'max_bits':maximum,'seconds':report['seconds']}))
if __name__=='__main__':main()

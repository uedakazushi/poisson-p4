#!/usr/bin/env python3
"""Read-only source/index/support preflight for the native128 first derivative."""
import collections,hashlib,json,gzip
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63'
PATHS={
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'source105':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-source105-Artin20-values.json',
 'source105_peer':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW/report.json',
 'interface':'reports/AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63/actual-source-physical-firstjet-interface.json',
 'Q_operator':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
 'proper':'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63/literal-proper-source14-preflight.json',
 'equations':'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_20261104_v63/actual-proper-source14-equations17.json',
 'first5':'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63/actual12-flag-divided-polynomials.json',
 'second5':'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json',
 'membership':'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/membership-and-source-cofactors.json',
 'core_source':'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json',
 'source_group_order':'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz',
 'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/phase2-native-adapter.json',
 'X0_report':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 d={k:json.loads(gzip.decompress((ROOT/p).read_bytes()).decode() if p.endswith('.gz') else (ROOT/p).read_text()) for k,p in PATHS.items()}
 assert d['source105_peer']['all_passed'] and d['X0_report']['all_passed']
 S,U,I,O=d['source'],d['source105'],d['interface'],d['Q_operator']
 unknown={i:dict((k,Q(c)) for k,c in a) for i,a in zip(U['joint105_unknown_order'],U['all105_unknown_final_values'])}
 first={i:[str(p.get(j,Q(0))) for j in range(9)] for i,p in unknown.items()}
 assert all(p[0]=='0' for p in first.values())
 params=[];roots=[]
 for fam,f in S['families'].items():
  pairs=[]
  for name,(a,b) in f['coefficient_unknown_pairs'].items():
   pairs.append({'name':name,'node_pair':[a,b],'even_first9':first[a],'odd_first9':first[b]})
  params.append({'family':fam,'pairs23':pairs,'s_node':f['p_quadratic_relation']['s'],
   's_first9':first[f['p_quadratic_relation']['s']],
   't_node':f['p_quadratic_relation']['t'],'t_first9':first[f['p_quadratic_relation']['t']]})
  roots.extend(i for pair in f['complete_normalized_physical_w_coefficients'][0] for i in pair)
 nodes=S['DAG_scalar_nodes'];closure=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in closure:continue
  closure.add(i);n=nodes[i]
  if n[0] in ['a','m','n','i']:todo.extend(n[1:])
 order=sorted(closure);uses=collections.Counter()
 for i in order:
  n=nodes[i]
  if n[0] in ['a','m','n','i']:uses.update(n[1:])
 uses.update(roots);live=set();peak=0
 for i in order:
  live.add(i);peak=max(peak,len(live));n=nodes[i]
  if n[0] in ['a','m','n','i']:
   for j in n[1:]:
    uses[j]-=1
    if uses[j]==0:live.remove(j)
 group=O['group30_to120'];odd=list(map(Q,group['group30_odd_origin']))
 assert odd==[Q(row[7]) for row in I['group15_first9']]
 assert group['original_group30_det_origin']=='4/9'
 packs=[('F17',d['equations']['equations14variable_sparse'],5),
  ('H12',[r['divided_primitive_polynomial'] for r in d['first5']['all12_actual_remaining_and_divided_records']],5),
  ('firstroots',[r['numerator'] for r in d['first5']['five_explicit_rational_roots'].values()],4),
  ('secondroots',[r['numerator'] for r in d['second5']['five_actual_rational_roots'].values()],4),
  ('membershipC',[p for row in d['membership']['cofactor3_by4'] for p in row],4)]
 polycounts={name:{'polynomials':len(ps),'raw_terms':sum(map(len,ps)),
   'needed_terms':sum(sum(sum(t[0])<=n for t in p) for p in ps),
   'degree_cutoff':n,'max_terms':max(sum(sum(t[0])<=n for t in p) for p in ps)} for name,ps,n in packs}
 E={
  'actual_Rbar_coordinates':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],
  'dual_index_order':'0 constant;1..7 zA0..6;8 rA; all epsilon products0',
  'all105_existing_unknown_first9':[[i,first[i]] for i in U['joint105_unknown_order']],
  'both48_parameter_firsts':params,
  'physical_R0_rows9':U['actual_physical_R0_rows'],
  'physical_R0_readback_first9':[[str(dict((i,Q(c)) for i,c in p).get(j,Q(0))) for j in range(9)] for p in U['actual_physical_R0_readback']],
  'nu_A_B_first9':[[str(dict((i,Q(c)) for i,c in p).get(j,Q(0))) for j in range(9)] for p in U['nu_CA_CB_values']],
  'all_physical150_pair_roots_CA_CB':{fam:f['complete_normalized_physical_w_coefficients'][0] for fam,f in S['families'].items()},
  'dual9_source_evaluation':{'closure_order':order,'source_ops':dict(collections.Counter(nodes[i][0] for i in order)),
   'peak_live_vectors':peak,'peak_dual9_coefficients':peak*9,
   'meaning':'evaluate original rational scalar coefficient DAG with these saved105firsts; no IFT solve or physical source jet'},
  'frame_input':'H_i=G_gamma_i*core_group_i(p_i), p_i=tau_i+s_i/2; frozen core_group_i is quadratic in p_i; all gamma pairfirsts are above',
  'GL18_source_parameter_positions':d['source_group_order']['GL18_generator_indices'],
  'frozen_core_group_matrices':{f:r['core_group'] for f,r in d['core_source']['families'].items()},
  'nil_filtered_polynomial_support':polycounts,
  'finite_algebras':{
   'gauge_J_AD':'QQ[eta,delta]/(eta2,delta2), length4,15separate groupdirections',
   'gauge_RHS':'dual9[eta]/(eta2-nuA1), length18, maximalideal fourthpower0',
   'full_native':'dual9[tauA,tauB]/(tauA2-nuA1,tauB2-nuB1), length36, maximalideal fifthpower0',
   'all_mixed_scalar_indices':'epsilon_l*tauA^a*tauB^b, l=1..8,a,b in{0,1}; internal tau powers2/3 reduced using nu1; never frozen to tau2=0 at first order'},
  'group30_derivative':{
   'unknown_order':'15 even then15 odd;30x8=240 QQ first coefficients',
   'saved_center_even15':[0]*15,'saved_center_odd15':list(map(str,odd)),
   'A0_stored15_Jacobian':I['group15_unit_Jacobian'],
   'A1_needed225_indices':[[constraint,parameter,'eta*delta'] for constraint in range(15) for parameter in range(15)],
   'A1_evaluation':'at saved eta-nil center, vary one G15 coefficient by delta, normalize own flag and compute the SAME15 constraints; extract eta*delta',
   'constraints':'five own flag plus ten fixed source-fiber coefficients, exact same source callback as origin',
   'J30':'[[A0,0],[A1,A0]], det4/9; do not assume A1zero',
   'rhs_new240_indices':[[constraint,eta,l] for eta in [0,1] for constraint in range(15) for l in range(8)],
   'rhs':'evaluate saved-center gauge with actualH_A/P_A firsts in length18; solve once with fixed30Jac; re-evaluate all30componentconstraints',
   'source_scalar_derivative_job_executed':False},
  'off_intersection_identities':{
   'F17':'F17(actual14)=R17*physicalF12, not0 in general',
   'first5':'full ordered monomial telescopes actualpivot->rationalroot5; derivative shortcut invalid awayfromR0origin',
   'R12':'d*H12=R12*physicalF12; d never inverted',
   'T':'T*H12=(G8,0)',
   'N':'N*physicalF12=d*(G8,0), not0 in general',
   'Qsyzygy':'Res4(N^t*S12by32), same native source signs, verify all actualN0coeffs match frozenX0'},
  'source_support_ceiling':{'new_J30_mixed_QQ':225,'new_RHS_QQ':240,'new_solution_QQ':240,
   'largest_output_scalar_slots':48*128*9,
   'all_diagnostics_conservative_QQ_slots':120000,'byte_guard':32<<20,
   'whole_sourcejet_or_105solve_or_GB_or_point_or_prime':False},
  'nil_center_reproduction':{'all144_nativeN_entries4components':'must match saved X0 exactly',
   'all8_G_C0':'must match saved X0, including nonzero tauA*tauB coefficients',
   'G8_tauA_tauB_coefficients':[p[3] for p in d['X0']['G8_C0']]},
  'requested_output_schema':{'key':'source_columns128_rows48_sparse_dual9',
   'columns':128,'row_range':[0,48],'entry':'[row,[[dualindex,QQstring],...]]',
   'bottom8rows':'identically0 in fullX','data_sha_must_be_pinned_in_report':True},
  'all_actual_coefficient_inputs_present':True,
  'actual_native128_origin_first_array_computed':False,
  'input_hashes':{k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}}
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=True)
 ep=out/'actual-native128-first-input-and-mixed-index-profile.json';ep.write_text(json.dumps(E,separators=(',',':'))+'\n')
 rep={'all_passed':True,'source_index_profile_only':True,'all_inputs_existing105firsts_and_literal_source_operators':True,
  'new_fixed30_mixed_Jacobian_coefficients':225,'new_first_order_rhs_and_unknown_coefficients':240,
  'native_algebra_length':36,'native128_derivative_evaluation_executed':False,
  'off_intersection_identity_corrections_explicit':True,
  'source_dual9_closure':len(order),'source_peak_coefficients':peak*9,
  'input_hashes':E['input_hashes'],'evidence':{'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)},'code_sha256':sha(Path(__file__))}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':main()

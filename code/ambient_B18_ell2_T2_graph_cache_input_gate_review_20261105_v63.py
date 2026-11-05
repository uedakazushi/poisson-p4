#!/usr/bin/env python3
"""Independent cache-only invariant/A23 availability and sensitivity gate."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_GRAPH_CACHE_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','gate':f'reports/{A}/actual-combined-invariant-and-A23-cache-source-input.json',
'factor_peer':'reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63_REVIEW/report.json',
'known':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json',
'mini':'reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json',
'fold':'reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63/actual-two-rell-source-folds-and-retained-native-first-operators.json',
'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
'AD':'reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_AD_20261105_v63/actual-two-by619-old-formal-rell-second-sensitivities.json',
'AD_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63_REVIEW/report.json',
'atom':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
'catalog':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
'old':'reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63/actual104-source-rell-values-and-two-by619-old-knownparts.json',
'old_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63_REVIEW/report.json',
'mu_gate':'reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',
'native':'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json',
'native_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63_REVIEW/report.json',
'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
'phi':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',
'physical':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
'physical_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};D={k:json.loads((ROOT/p).read_text()) for k,p in P.items() if k!='C0'};assert D['primary']['all_passed'] and all(D[k]['all_passed'] for k in D if k.endswith('peer'))
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 E=D['gate'];ops=D['chain']['operator_scalar_nodes'];B=list(map(dict,D['chain']['B48_by640_sparse_columns']));C0=fields(ROOT/P['C0'],['honest_D0_relation_columns','actual_C0_basis8_scalar_nodes']);Pcols=C0['honest_D0_relation_columns']
 grads={i:sp(g) for i,g in D['AD']['typed_node_original_second_gradients']};orig={i:sp(p).get(0,Q(0)) for i,p in D['atom']['all79_atom_values']};catalog={i:sp(g) for i,v,g in D['catalog']['all108_composite_values_and_opaque_second_gradients']};oldlabels={i for i,v,g in D['old']['all104_original_C0_values_and_opaque_gradients']}
 gamma=list(map(sp,D['fold']['two_Gamma_on56_actual_X2']));S0=list(map(sp,D['mu_gate']['native_S0_columns']));newroots=set();newdirect=set();records=[]
 for index,j in enumerate(D['known']['selected_bottom_columns']):
  c=sp(D['mini']['all642_zero_constant_combinations'][j]['original768_combination']);stored=E['all23_old_source_sensitivity_proofs'][index];assert stored['basis_column']==index and stored['source642']==j and sp(stored['original768'])==c
  cn={k-640:a for k,a in c.items() if k>=640};assert not mv(S0,cn)
  for channel,g in enumerate(gamma):
   f={}
   for source,a in c.items():
    if source>=640:continue
    for row,b in g.items():
     if row>=48:f=add(f,{Pcols[source][row-47]:a*b});continue
     if row not in B[source]:continue
     root=B[source][row];n=ops[root]
     if n[0]=='c':continue
     if n[0]=='actual_C0_scalar_function':f=add(f,{n[1]:a*b});continue
     if root not in grads:
      assert n[0]=='product' and n[1:3]==[-1,1] and len(n[3])==2
      assert all(ops[x][0]=='actual_C0_scalar_function' and orig[ops[x][1]]==0 for x in n[3]);grads[root]={};newroots.add(root)
     f=add(f,scale(grads[root],a*b))
   assert enc(f)==stored['two_literal_old_second_functionals_after_typed_AD'][channel]
   assert all(i==0 or i in catalog and not catalog[i] for i in f)
   newdirect.update(set(f)-oldlabels)
  records.append({'source642':j,'native_S0_cN_zero':True,'both_honest_second_sensitivities_zero':True})
 assert newroots=={516,785,846} and sorted(newdirect)==E['new_direct_C0_labels_are_literal_coordinates']==[7148,7688,7943]
 assert newdirect<=set(C0['actual_C0_basis8_scalar_nodes'][:7])
 assert [[i,ops[i],[[j,ops[j],str(orig[ops[j][1]])] for j in ops[i][3]]] for i in sorted(newroots)]==E['three_new_typed_root_zero_second_sensitivity']
 # Original trace IDs, keptdistinctfromphysical F9(sA,rA,ell).
 nuids=[[D['source']['families'][f]['p_quadratic_relation'][k] for k in ['s','t']] for f in ['CA','CB']];assert nuids==[[11,12],[20537,20538]]==E['raw33_primitive_source_formula']['actual_monic_trace_s_t_nodes']
 roots=dict(D['phi']['all_requested164_root_values'])
 for i,p in zip(D['phi']['joint105_unknown_order'],D['phi']['all105_unknown_final_values']):
  if i in roots:assert sp(roots[i])==sp(p)
  roots[i]=p
 assert all(i in roots for p in nuids for i in p) and D['phi']['all105_equations_exactly_zero']
 physical=D['physical']['physical12_coefficients'];assert physical[5]['physical_row']==67 and physical[7]['physical_row']==63
 assert E['raw33_primitive_source_formula']['Wv']==[[23,'-3375/11'],[29,'-50'],[30,'25/2']]
 # Exactfree4 multiplication: coeffindexmask0/1/2/3, tauA²=nuA,tauB²=nuB.
 # Symbolic scalar labels avoid querying any sourcevalues.
 def term(*names):return tuple(sorted(names))
 def times_mask(coeff,inp):
  out=[{} for _ in range(4)]
  for mask,p in enumerate(coeff):
   carry=mask&inp;result=mask^inp
   for m,c in p.items():
    names=list(m)+(['nuA'] if carry&1 else [])+(['nuB'] if carry&2 else []);out[result]=add(out[result],{tuple(sorted(names)):c})
  return out
 F67=[{term('d67'):Q(1)},{term('a67'):Q(1)},{term('b67'):Q(-1)},{}];F63=[{term('d63'):Q(1)},{term('a63'):Q(1)},{term('b63'):Q(-1)},{}]
 terms=[times_mask(F67,3),times_mask(F63,1),times_mask(F63,2)];weights=[Q(-3375,11),Q(-50),Q(25,2)];odd=[]
 for mask in [1,2,3]:odd.append(add(*(scale(v[mask],c) for v,c in zip(terms,weights))))
 expected=[{term('nuB','b67'):Q(3375,11),term('d63'):Q(-50)},{term('nuA','a67'):Q(-3375,11),term('d63'):Q(25,2)},{term('d67'):Q(-3375,11),term('b63'):Q(50),term('a63'):Q(25,2)}];assert odd==expected
 limits=E['limits'];assert E['whole_cache_slot_bound']==47182==sum(E['slot_breakdown'].values())<60000 and limits['bits']==2048
 ep=o/'independent-A23-sensitivity-and-primitive-cache-dot-source-gate.json';dump(o,ep.name,{'inputs':inputs,'all46_actual_selected23_old_callbacks_refolded_and_opaque_second_sensitivities_zero':True,'all23_native_S0_cN_zero_explicitly_verified':records,'new_typed_three_products_zero_second_sensitivity':True,'new_direct_labels_three_literal_z_coordinates':[7148,7688,7943],
'primitive_Vw_exact_mask_law_verified':'(-50d63+(3375/11)nuB*bB67, (25/2)d63-(3375/11)nuA*bA67, -(3375/11)d67+50bB63+(25/2)bA63)',
'original_monic_trace_ids':nuids,'trace_free_monic_proof':'tau=p-s_trace/2 andp²=s_trace*p+t implytau²=t+s_trace²/4. Thesephysicalsource105tracenodesarenotF9sA.',
'available_native_individual_pivot_ell_firsts':'L0*S1ell+L1ell*S0 fromtheacceptednativecache; neverreplaceindividualcolumnsbyeffectiveL0S1ellonly. WholezeroorigincombinationsdropL2butretainL1ellOmegaR.',
'whole_cache_slot_bound':47182,'declared_coefficient_bit_guard':2048,'resource_scope':'Thedeclared2048bitguardisreportedexplicitly; itisnotclaimedtofollowfromtheolder384bitgates. ParentmustuseitsexactdeclaredscopeforvalueGO.',
'no_A23_Phi_invariant_graph_or_new_physical_gauge_values_computed':True})
 report={'id':STEM,'all_passed':True,'all46_A23_and_actual_primitive_dot_inputs_source_bound':True,'all23_native_S0_eligibility_verified':True,'actual_A23_or_invariant_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()

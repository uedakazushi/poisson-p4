#!/usr/bin/env python3
"""Literal original five-head connecting operator, with formal RegDivell."""
import hashlib,json,re,signal,time
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_FIVE_FORMAL_SOURCE_20261105_v63'
PATHS={
 'q7':'reports/AMBIENT_B18_C_GENUINE_PROJECTION7_SOURCE_20261105_v63/sparse-actual-Dlinear-projector7-source-and-five-head-order-frontier.json',
 'q7_report':'reports/AMBIENT_B18_C_GENUINE_PROJECTION7_SOURCE_20261105_v63/report.json',
 'seven':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63/honest-seven-source-quotient-four-fixed-kernel-classes-and-edge-constraints.json',
 'seven_peer':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63_REVIEW/report.json',
 'unit':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'unit_peer':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63_REVIEW/report.json',
 'natural':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/actual-natural-Told-to-Tnil-source-map-and-origin.json',
 'natural_peer':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW/report.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'chain_peer':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW/report.json',
 'clones':'reports/AMBIENT_B18_NIL_CONNECTING16_SOURCE_OPERATOR_20261105_v63/actual-formal-connecting16-and-ten-source-obstructions.json',
 'clones_peer':'reports/AMBIENT_B18_NIL_CONNECTING16_SOURCE_20261105_v63_REVIEW/report.json',
 'clones_alias':'reports/AMBIENT_B18_NIL_CONNECTING16_NAMESPACE_BINDING_20261105_v63/report.json',
 'R0':'reports/AMBIENT_B18_C_GENUINE_R0_TRANSPORT_SOURCE_20261105_v63/genuine-four-R0-action-transport-and-cached-HB-operator-identities.json',
 'R0_peer':'reports/AMBIENT_B18_C_GENUINE_R0_TRANSPORT_SOURCE_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def scalar_prefix(path,last,receipts):
 b=path.read_bytes();s=b.decode();key='operator_scalar_nodes'
 hits=list(re.finditer(r'"'+key+r'"\s*:',s));assert len(hits)==1
 pos=hits[0].end()
 while s[pos].isspace():pos+=1
 assert s[pos]=='[';start=pos;pos+=1;dec=json.JSONDecoder();rows=[]
 for i in range(last+1):
  while s[pos].isspace()or s[pos]==',':pos+=1
  row,pos=dec.raw_decode(s,pos);rows.append(row)
 lo=len(s[:start].encode());hi=len(s[:pos].encode())
 receipts.append({'path':str(path.relative_to(ROOT)),'field':key,'decoded_prefix_rows':len(rows),
   'byte_slice':[lo,hi],'slice_sha256':hashlib.sha256(b[lo:hi]).hexdigest(),
   'numeric_prefix_count':S.count(rows),'unused_suffix_and_B640_never_decoded':True})
 return rows
def main():
 started=time.monotonic();signal.alarm(60);receipts=[]
 pins={p:sha(ROOT/p)for p in PATHS.values()};pins[str(Path(S.__file__).relative_to(ROOT))]=sha(Path(S.__file__))
 fields={
  'unit':['actual_variable_order','V56_target_order','V56_column_order','typed_regular_function_operators',
    'actual_V56_function_operator_refs','typed_m24_by56','typed_n56_by24','typed_unit24_inverse_node',
    'fixed_selected24_columns','fixed24_minor_determinant','typed_map_proof','actual_Xfull_source'],
  'natural':['old_minimal_generator_free56_rows','old_minimal_honest116_lifts','matrix_operator_nodes','matrix_operator_shapes','source_artifacts','linearity_and_descent'],
  'chain':['exact_existing_Q48_C4','exact_existing_U9_C4','scalar_source_artifacts','actual_scalar_reference_artifact',
    'old_multiplier_last_is_f63_not_new_gA','typed_operator_nodes'],
  'clones':['two_distinct_original105_source_clones','source_namespace_bindings','formal_coordinate_division_proof'],
  'R0':['basis_change_inverse4','actual_r_matrix4','actual_z0_through_z6_and_ell_actions'],
  'q7':['projector7_by660_sparse_columns','all220_monomials_degree_at_most3','same_actual_R0_variables','kernel_projection'],
 }
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 assert d['q7_report']['all_passed']and not d['q7_report']['edge_values_or_new_affine_responses_computed']
 c=d['chain'];u=d['unit'];n=d['natural']
 assert c['actual_scalar_reference_artifact']=='C0'and c['old_multiplier_last_is_f63_not_new_gA']
 Q48=c['exact_existing_Q48_C4'];U9=c['exact_existing_U9_C4']
 assert len(Q48)==48 and all(len(x)==4 for x in Q48)and len(U9)==9
 assert U9[0]==[1,0,0,0]and all(U9[j]==[0,j+1,0,0]for j in range(1,9))
 # Full V56 has the three nonconstant C4 components, with the ACTUAL
 # negative U sign. Reference pairs are [original_scalar_id, +/-1].
 full=[[[i,1]for i in row[1:]]for row in Q48]
 full +=[[[i,-1 if i else 1]for i in row[1:]]for row in U9[1:]]
 assert len(full)==56
 ops=u['typed_regular_function_operators'];oldrefs=u['actual_V56_function_operator_refs']
 for col in range(56):
  for target in range(3):
   op=ops[oldrefs[col][target]]
   got=[0,1]if op==['zero']else[op[1],-1 if op[0]=='negative_actual_chain_scalar_restrict_ell0'else 1]
   assert op==['zero']or op[0]in['actual_chain_scalar_restrict_ell0','negative_actual_chain_scalar_restrict_ell0']
   assert got==full[col][target],(col,target,got,full[col][target])
 roots=sorted({i for col in full for i,sign in col});assert max(roots)==2160
 prefix=scalar_prefix(ROOT/PATHS['chain'],max(roots),receipts);closure=set()
 def visit(i):
  if i in closure:return
  closure.add(i);op=prefix[i]
  if op[0]=='c':assert len(op)==3 and op[2];return
  if op[0]=='actual_C0_scalar_function':assert len(op)==2;return
  if op[0]=='sum':
   assert len(op)==2
   for j,num,den in op[1]:assert den;assert j<i;visit(j)
  elif op[0]=='product':
   assert len(op)==4 and op[2]
   for j in op[3]:assert j<i;visit(j)
  else:raise AssertionError(op)
 for i in roots:visit(i)
 selected=[[i,prefix[i]]for i in sorted(closure)];prefix_peak=S.count(prefix)
 del prefix
 heads=n['old_minimal_generator_free56_rows'];assert heads==[22,23,29,30,55]
 nn=n['matrix_operator_nodes'];ns=n['matrix_operator_shapes']
 assert nn[0]==['actual_external_matrix','old_unit','typed_m24_by56']
 assert nn[1]==['actual_external_matrix','old_unit','typed_n56_by24']
 assert ns[:8]==[[24,56],[56,24],[56,5],[56,56],[56,56],[56,56],[56,56],[56,5]]
 assert nn[2][0]=='QQ_matrix_rows'
 E=nn[2][1];assert E==[[str(int(i==h))for h in heads]for i in range(56)]
 assert n['old_minimal_honest116_lifts']=='w84=0; free32 literal tail vectors; both old Schur codomain right inverses preserve these lifts exactly'
 assert all(0<=i<len(ops)for row in u['typed_m24_by56']+u['typed_n56_by24']for i in row)
 assert len(u['typed_m24_by56'])==24 and all(len(r)==56 for r in u['typed_m24_by56'])
 assert len(u['typed_n56_by24'])==56 and all(len(r)==24 for r in u['typed_n56_by24'])
 clones=d['clones']['two_distinct_original105_source_clones'];assert len(clones)==2
 assert[cl['namespace']for cl in clones]==['FULL_R0_F9','ELL0_S_F9']
 assert clones[0]['target9'][-1]==['sum',['r'],['scale','1/5',['ell']]]
 assert clones[1]['target9'][-1]==['r']
 assert clones[0]['equation_target_indices']==list(range(48,57))==clones[1]['equation_target_indices']
 assert clones[0]['origin_Jacobian_determinant']==clones[1]['origin_Jacobian_determinant']=='1202315964973056/390625'
 assert d['clones_alias']['all_passed']
 # Every operation has an explicit shape and foreign namespace. No expanded
 # scalar coefficient, 56x768 matrix, or higher jet is evaluated.
 nodes=[];shapes=[]
 def put(op,shape):nodes.append(op);shapes.append(shape);return len(nodes)-1
 m=put(['actual_external_ELL0_matrix','old_unit','typed_m24_by56','ELL0_S_F9'],[24,56])
 nv=put(['actual_external_ELL0_matrix','old_unit','typed_n56_by24','ELL0_S_F9'],[56,24])
 e=put(['QQ_matrix_rows',E],[56,5]);identity=put(['identity_matrix',56],[56,56])
 nm=put(['matrix_product',nv,m],[56,56]);complement=put(['matrix_difference',identity,nm],[56,56])
 sigma=put(['matrix_product',complement,e],[56,5])
 lift=put(['ell_independent_ring_section_Rbar_into_R0',sigma,'DO_NOT_REPLACE_BY_FULL_SOURCE_m_n'],[56,5])
 vf=put(['actual_full_V3_by56_from_original_chain_Q48_negative_U8',full,'original_chain','FULL_R0_F9'],[3,56])
 numerator=put(['matrix_product',vf,lift],[3,5])
 eta=put(['entrywise_formal_RegDiv_ell_before_finite_quotient',numerator,'ell','Vbar_sigma=0'],[3,5])
 coefficients=put(['raw_coefficient_extraction_through_degree3',eta,'actual_z7_r_ell','ordinary_coefficients'],[660,5])
 projected=put(['actual_external_Dlinear_q7_sparse','q7','projector7_by660_sparse_columns',coefficients],[7,5])
 coords=put(['select_rows_of_genuine_ker_E7',projected,[1,2,5,6],'eta_is_genuine'],[4,5])
 inv=put(['QQ_matrix_rows',d['R0']['basis_change_inverse4']],[4,4])
 jordan=put(['matrix_product',inv,coords],[4,5])
 constraints=put(['actual_three_chi_constraints_from_source7',jordan,'chi=(0,0,c/4,c,d)'],[3,1])
 required=[]
 for col in d['q7']['projector7_by660_sparse_columns']:
  take=[[i,a]for i,a in col['target7_affine']if i in[1,2,5,6]]
  if not take:continue
  exp=col['monomial'][:];exp[8]+=1
  required.append({'raw_target':col['raw_target'],'numerator_monomial':exp,'q7_free_rows_affine':take})
 highest=[r for r in required if sum(r['numerator_monomial'])==4]
 assert len(highest)==4 and{r['raw_target']for r in highest}=={0,1}
 assert{tuple(r['numerator_monomial'][7:])for r in highest}=={(3,1),(2,2)}
 current=S.count(d)+S.count(selected)+S.count(receipts)+S.count(nodes)+S.count(shapes)+S.count(required)+S.count(full)+S.count(closure)+12000
 owned=max(current,current+prefix_peak);assert owned<=150000
 data={
  'scope':'SOURCE-only actual original connecting delta:Told11→genuine C[ell]. No edge scalar values; not a new formal lift of every nil generator.',
  'source_namespace_bindings':{k:{'path':p,'sha256':pins[p]}for k,p in PATHS.items()},
  'original_chain_scalar_source_artifacts':c['scalar_source_artifacts'],
  'full_V56_reference_pairs_columns56_targets3':full,
  'full_to_ELL0_all168_literal_source_sign_bindings_verified':True,
  'selected_original_chain_scalar_nodes':selected,
  'all_selected_scalar_nodes_are_original_immutable_records':True,
  'old_eighth_generator':'literal old f63; NOT replace by new curve gA without the accepted allorder z-column operation',
  'same_original105_two_separate_source_clones':clones,
  'old_five_heads_free56_rows':heads,
  'actual_old_minimal_to_cycle_lift':'sigma=(I56−nbar*mbar)E56×5. Accepted old honest116 lifts have w84=0 and literal free32 tails, preserved by both Schur sections. This is the literal source used in the accepted natural map, not an arbitrary choice of jets.',
  'matrix_operator_nodes':nodes,'matrix_operator_shapes':shapes,
  'matrix_roots':{'five_formal_cycles_ELL0':sigma,'ell_independent_lifts':lift,'full_actual_V3_by56':vf,
    'full_numerators':numerator,'formal_RegDivell_eta':eta,'genuine_q7_coordinates':coords,'J2_plus_two_k_five_head_edge':jordan,'three_chi_constraints':constraints},
  'all15_regular_division_memberships':'The unit construction gives mbar*nbar=I24 and Vbar=d1*mbar. Hence Vbar*(I−nbar*mbar)E=0 entrywise exactly, so every numerator entry belongs to ell*R0. RegDivell is in the formal regular ring BEFORE any truncation/normal form. No division by ell is performed in C4 or another Artin ring.',
  'genuine_target_proof':'ell*eta=Vfull*lift(sigma) lies in imVfull; hence [eta]∈C[ell]. Actual D-linear q takes it into kerE7. Free rows1,2,5,6 give its genuine four coordinates; the accepted inverse basis change gives J2+k² coordinates. q is not asserted R0-linear on arbitrary C.',
  'full_boundary_descent':'Full comparison X has640 columns(B48,P8) with V X=QB−UP=0, plus128 complete native syzygies(S128,0). These are the actual same-R0 free complex, with constant target1 eliminated by its exact identity. For a mod-ell boundary, lift the FULL original X and its coefficient vector: the ell-independent representative differs by ell*t, so eta differs by Vfull*t, zero in C. Thus all768 original relation paths, and their Schur/Euler images825, descend. No uncontrolled homology tensor or missing higherTor assumption is used.',
  'Rbar_linearity':'The ell-independent lift is the literal ring section QQ[[z7,r]]→QQ[[z7,r,ell]]. Formal RegDivell and the boundary argument give the usual Rbar-linear connecting map; actual target transport, not naive q action, supplies its J2+k² structure.',
  'honest_kernel_reduction':'All values factor through the already peered honest M7=Told/(z7,r²)Told. ITold length4 is an ACTUAL kernel subspace. Twenty head coordinates obey three independent chi constraints, leaving seventeen unknown parameters; none is assigned here.',
  'required_raw_numerator_coefficient_functionals':required,
  'highest_four_raw_numerator_monomials':highest,
  'highest_order_status':'V4*sigma0, V3*sigma1, V2*sigma2, V1*sigma3 are all retained. A later literal plane-order gate may cancel selected terms; this construction does NOT assume it.',
  'no_new_scalar_response_or_physical_gauge_IFT_or_jet':True,
  'selective_source_extraction_receipts':receipts,
  'whole_source_owned_numeric_upper':owned,
  'caps':{'owned_numeric':150000,'seconds':60,'memory_bytes':2147483648,'RAW_bytes':33554432},
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'actual-five-head-formal-RegDivell-operator-and-source-binding.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_FIVE_ACTUAL_HEAD_CYCLES_FORMAL_REGDIVELL_GENUINE4_NO_VALUES',
  'all168_full_to_ell0_reference_bindings_verified':True,'regular_formal_division_memberships':15,
  'actual_original_connecting_values_computed':False,'highest_fourth_response_cancellation_proved':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'whole_numeric_upper':owned,'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)
if __name__=='__main__':main()

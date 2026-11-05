#!/usr/bin/env python3
"""Exact SOURCE three genuine bar-pair operators; all higher jet data formal."""
import hashlib,json,signal,time
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KSQUARE_THREE_PAIR_OPERATOR_SOURCE_20261105_v63'
PATHS={
 'natural':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/actual-natural-Told-to-Tnil-source-map-and-origin.json',
 'natural_peer':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW/report.json',
 'zero3':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_MAP_20261105_v63/actual-three-first-columns-and-fixed-quotient-proof.json',
 'zero3_peer':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_SOURCE_20261105_v63_REVIEW/report.json',
 'seven':'reports/AMBIENT_B18_C_KBAR_ACTUAL_SEVEN_BASIS_SOURCE_20261105_v63/actual-seven-barbasis-typed56-cycle-operator-and-full11-unit-complement.json',
 'seven_peer':'reports/AMBIENT_B18_C_KBAR_ACTUAL_SEVEN_BASIS_SOURCE_20261105_v63_REVIEW/report.json',
 'three':'reports/AMBIENT_B18_C_KBAR_THREE_GENERATOR_RANK_SOURCE_20261105_v63/report.json',
 'three_peer':'reports/AMBIENT_B18_K_THREE_ROOT_QUOTIENT_REVIEW_20261105_v63/report.json',
 'finite_nil_peer':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json',
 'old_exact':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json',
 'q2':'reports/AMBIENT_B18_C_KBAR_Q2_BAR_ANNIHILATOR_SOURCE_20261105_v63/report.json',
 'H_order':'reports/AMBIENT_B18_C_KBAR_Q2_BAR_ANNIHILATOR_SOURCE_20261105_v63/specific-H-template-order-and-raw-domain-companion.json',
 'point':'reports/AMBIENT_B18_NBAR14_KERNEL_SOURCE_PREFLIGHT_20261105_v63/five-cache-functionals-and-genuine-fourteen-kernel-source.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'minimal_peer':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63_REVIEW/report.json',
 'honest':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63/actual-honest-Tnil128-972-formal-operator.json',
 'honest_peer':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63_REVIEW/report.json',
 'product':'reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63/actual-nil-product44-and-cube-DAG.json',
 'product_peer':'reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63_REVIEW/report.json',
 'upper':'reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63/actual-trace-first-gate-and-entire-conductor-proof.json',
 'upper_peer':'reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={
  'natural':['matrices','actual_origin_map_columns5'],
  'zero3':['source_frames','actual_three_minimal_first_columns40','actual_old_free56_labels','all576_actual_first_Koszul_checks'],
  'seven':['actual_correction_coefficients'],
  'point':['minimal_head_inclusion_rows128','literal_direct_section_path','head42_origin_cycle44_sparse','all24_m0_e43_entries_zero'],
  'minimal':['matrices','operator_shapes'],
  'honest':['matrices','matrix_operator_shapes','shared_homotopy_template_v1','ell0_pullback'],
  'product':['all_new_R0_scalar_nodes','literal_bilinear44_output_roots','actual_beta8_C0_node_ids','cycle_representative_even_b0','cycle_representative_odd_scalar_B'],
  'upper':['J_annihilates_K_squared','actual_joint_ideal','m_R0_square_subset_J','Nil_D0_times_K_squared_zero'],
 }
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 prior=S.count(d)+20000;assert prior<150000
 assert d['zero3']['actual_old_free56_labels']==[23,29,30]and d['zero3']['all576_actual_first_Koszul_checks']
 assert d['natural']['matrices']['actual_natural_map_minimal5_by5']==15
 assert d['minimal']['matrices']['direct_unit_section_n']==8 and d['minimal']['matrices']['full_target_projection5_by128']==59
 assert d['honest']['matrices']['new_m24_by44']==4 and d['honest']['matrices']['d2_new']==23
 assert d['point']['minimal_head_inclusion_rows128']==[42,43,50,68,70]
 assert d['product']['actual_beta8_C0_node_ids'][-1]==505938
 inputs={x[1]:[]for x in d['product']['all_new_R0_scalar_nodes']if x[0]=='arbitrary_cycle_input'}
 for i,x in enumerate(d['product']['all_new_R0_scalar_nodes']):
  if x[0]=='arbitrary_cycle_input':inputs[x[1]].append([x[2],i])
 assert set(inputs)=={'z','zp'}and all(sorted(x for x,i in v)==list(range(44))for v in inputs.values())
 assert len(d['product']['literal_bilinear44_output_roots'])==44
 assert d['upper']['J_annihilates_K_squared']and d['upper']['Nil_D0_times_K_squared_zero']
 # Actual old-five generator coordinates of U,S,q2. The entries are exact
 # ordinary S-polynomials; no natural-map coefficient is evaluated.
 coeff=d['seven']['actual_correction_coefficients'];z=[0]*8;r=z.copy();r[7]=1;rr=z.copy();rr[7]=2
 T=[[[]for _ in range(3)]for _ in range(5)]
 # Use structured signs, not a malformed doubleminus rational string.
 from fractions import Fraction as Q
 T[1][0]=[[z,str(-Q(coeff['a3=b3/b1']))]];T[1][1]=[[z,str(-Q(coeff['aS=(b2-b3/4)/b1']))]]
 T[2][1]=[[z,'1']];T[3][0]=[[z,'1']];T[3][1]=[[z,'-1/4']]
 T[4][0]=[[r,str(-Q(coeff['c3=(p3-a3*p1)/d']))]];T[4][1]=[[r,str(-Q(coeff['cS=(p2-p3/4-aS*p1)/d']))]];T[4][2]=[[rr,'1']]
 E=[[int(i==j)for j in[42,43,50,68,70]]for i in range(128)]
 nodes=[];shapes=[]
 def put(op,shape):nodes.append(op);shapes.append(shape);return len(nodes)-1
 gamma=put(['external_actual_Rbar_matrix','natural','matrix_operator_nodes',15],[5,5])
 kgen=put(['actual_oldfive_polynomial_coordinates',T,'original8_z0_to6_r'],[5,3])
 Gamma=put(['matrix_product',gamma,kgen],[5,3])
 m=put(['external_actual_Rbar_matrix','honest','matrix_operator_nodes',4],[24,44])
 n=put(['external_actual_Rbar_matrix','minimal','operator_nodes',8],[44,24])
 d2=put(['external_actual_Rbar_matrix','honest','matrix_operator_nodes',23],[24,84])
 inc=put(['QQ_selector128_by5',E],[128,5]);ef=put(['row_slice',inc,0,44],[44,5]);ew=put(['row_slice',inc,44,84],[84,5])
 ident=put(['identity_matrix',44],[44,44]);nm=put(['matrix_product',n,m],[44,44]);complement=put(['matrix_difference',ident,nm],[44,44])
 free=put(['matrix_product',complement,ef],[44,5]);wed=put(['matrix_product',n,put(['matrix_product',d2,ew],[24,5])],[44,5])
 section=put(['matrix_sum',free,wed],[44,5]);cycles=put(['matrix_product',section,Gamma],[44,3])
 Pi=put(['external_actual_Rbar_matrix','minimal','operator_nodes',59],[5,128])
 outputs=[];allpairs=[]
 for i,j in[(0,0),(0,1),(1,1)]:
  left=put(['column',cycles,i],[44,1]);right=put(['column',cycles,j],[44,1])
  prod=put(['apply_exact_external_Product44_Rbar','product','all_new_R0_scalar_nodes','literal_bilinear44_output_roots',
   {'z':left,'zp':right},'ALLactual_R0_functions_functorially_restrict_via_honest_ell0_pullback'],[44,1])
  mc=put(['matrix_product',m,prod],[24,1]);lift=put(['apply_NON_Rbar_LINEAR_SHARED_H_TO_WHOLE_FUNCTION_VECTOR',mc,[1,1,2],
   'honest.shared_homotopy_template_v1','each_exact_cycle','NOcoefficientpullout_or_finiteNF'],[84,1])
  honest=put(['vertical_concatenate',[prod,lift]],[128,1]);minimal=put(['matrix_product',Pi,honest],[5,1])
  rmin=put(['Rbar_scalar_multiply','rA',minimal],[5,1]);outputs.extend([minimal,rmin]);allpairs.append({'pair':[i,j],'free44_cycle':prod,'honest128':honest,'minimal5_representative':minimal,'r_multiple':rmin})
 data={
  'scope':'Exact full-source canonical THREE pair operator plusrmoduleclosure, NO product coefficients or rank evaluated. Coverage bythreepairs uses q2-ann orderlemma onlyafterindependentadoption.',
  'source_namespace_bindings':{k:{'path':p,'sha256':pins[p]}for k,p in PATHS.items()},
  'three_actual_K_generators':['U','S','q2'],'actual_oldfive_coordinates5_by3':T,
  'operator_nodes':nodes,'operator_shapes':shapes,'exports':{'canonical_K_free44_cycles3':cycles,'pair_and_r_minimal5_representatives6':outputs},
  'all_three_pairs_and_full_source_routes':allpairs,
  'Product44_input_node_labels_original_indices':inputs,
  'canonical_resection_proof':'ExactthreeunitSchur isomorphism usesliteralretainedraw128selector. ApplyingPhi_direct givesfree44cycle(canonicalminclass). It agreeswiththefullactualoldnaturalcycleclass via acceptedsectionchange/boundaries. All3KclassesareinTRUEkerconnecting, soactualKliftsEXIST; applyingProduct44 modellcomputesactualNbarproducts. NoTnilalgebra asserted.',
  'all_normalizer_carry_retained':'The complete492-nodeactualProduct44 withQrad_even36,beta8,nuA isappliedafterTHE SAMEactualell0sourcepullback. No unit/gauge/kappa/C4coefficient is dropped; no oldf63forliteralcurve gA505938 substitution.',
  'homotopy_rule':'Eachwholem*c is anexact24cycle; applysharedH_n sectors1,1,2 BEFOREfiniteNF. Never S-linearly extendh780or pullactualcoefficientsoutofH. raw128freepartisc exactly;Pi59 isfullallorderprojection.',
  'complete_module_closure_conditional':'If q2*Nbar0 sourceisindependentlyadopted, Kbar² isRbar-generatedbyU²,US,S². AcceptedJK²0 givesz7*classes0 andr²*classes0 underell0. ThusitsWHOLEQQimageequals spanoftheSIXlistedclasses(pair,rpair). Noothermonomialmultiplesaremissing.',
  'known_first_input_cache':{
   'exact_corrected_ZERO3_firsts40_by3':d['zero3']['actual_three_minimal_first_columns40'],
   'natural_origin_columns5':d['natural']['actual_origin_map_columns5'],
   'Gamma1_U':'ZERO3(h3)-a3*ZERO3(h1)-c3*r*M0(h4)',
   'Gamma1_S':'ZERO3(h2)-ZERO3(h3)/4-aS*ZERO3(h1)-cS*r*M0(h4)',
   'Gamma2_q2_lead':'r²*M0(h4)',
   'canonical_free44_firsts':'section0*Gamma1;section0freeheads are(I-n0m0)e42,e43 andthreewedgeorigins0. Thisisacanonicalrepresentative, NOTzeroorigin=>mNbar.',
   'full_Γ2_Γ3_uncomputed':True},
  'conservative_jet_frontier_to_minimal_degree4':{
   'canonical_input_cycles':'c_U,c_S throughdegree3;section throughdegree2 andGamma throughdegree3',
   'Product44_functions':'Qrad_even36 throughdegree3;CAfirst7fi throughdegree3;literal gA throughdegree3;nuA throughdegree2',
   'honest_output':'productfree44 degrees3..5;regularm degrees0..2;Happliedwholeandmayloseone',
   'full_minimal_projection':'Pi0,Pi1,Pi2 sincehonestpair startsdegree2, unlikecube/quartic-onlycase',
   'finite_class_correction':'Knownfiltered16basis alone isNOTanactualfullNF/actiontable. Needboundedcoherentadjoints/sourceclearingsoflowerrelations BEFORErank ontrueclasses; nofull2475x~140kArtinmatrix.'},
  'mandatory_UNKNOWNs_retained':['Gamma2','Gamma3','section1/2 inclm1/m2andn1/n2','Qeven2/3','fi2/3','nuA2','m1/m2','Pi1/Pi2','lowerrelationcorrectionstoactualfilteredNF'],
  'not_cached_yet_by_this_packet':'OnlyGamma1/origininputsareboundtocorrectedactualcache here. Nohigherdataavailabilityor sensitivitycancellationisclaimed. NextSOURCEmustfoldcompleteoutputadjoints andannihilateorbindeveryunavailablecallback before anyvalueworker.',
  'H0_bridge_positivecriterion_ONLY':'AllsixclassesareactualproductsinH0bar becauseJK²0. If theiractualQQrank is3 inTRUE Nbar/Tnil injection, mu_DH0=3 andellNakayama implyH0=D*K²; acceptedNK²0 thenNH00,N4zero,ell²H00. NOrank,H0typeor nilimprovementisclaimedhere.',
  'no_pair_class_or_higher_jet_value_computed':True,'no_signed_product_GO_from_sourcepacket':True,
  'selective_receipts':receipts,
 }
 owned=S.count(d)+S.count(data)+S.count(nodes)+S.count(shapes)+20000;assert owned<150000;data['whole_SOURCE_owned_upper']=owned
 data['next_positive_numeric_worker_gate_OPEN']=True
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'exact-three-pair-six-closure-class-source-operator-and-uncomputed-jet-frontier.json';p.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_EXACT_THREE_PAIR_OPERATOR_NO_PRODUCT_VALUES_HIGHER_CACHE_GATE_OPEN',
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{p.name:sha(p)},'whole_SOURCE_owned_upper':owned,'operator_nodes':len(nodes),'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()

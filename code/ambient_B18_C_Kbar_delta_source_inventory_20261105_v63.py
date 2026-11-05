#!/usr/bin/env python3
"""Read-only SOURCE inventory for the remaining coherent r-Jordan scalar."""
import hashlib,json,signal,time
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KBAR_DELTA_SOURCE_INVENTORY_20261105_v63'
PATHS={
 'three':'reports/AMBIENT_B18_C_KBAR_THREE_GENERATOR_RANK_SOURCE_20261105_v63/formal-deep-action-unit-three-generator-minor-and-exact-rank-four.json',
 'operator':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'operator_peer':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW/report.json',
 'linear_rank':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
 'linear_source':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
 'zero620':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/zero-linear620-source-combinations.json',
 'quadratic_rank':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63/actual-three-by620-quadratic-image-certificate.json',
 'old39':'reports/AMBIENT_B18_T_KNOWN_INITIAL_MODULE_LENGTH_UPPER24_20261104_v63/known-initial-module-and-length-bound.json',
 'quadratic464':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase2-certified-quadratic-input.json',
 'initial':'reports/AMBIENT_B18_C_KBAR_R2_INITIAL_SOURCE_20261105_v63/two-cached-alpha-reductions-and-projected-r2-seven-source-certificates.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={
  'operator':['actual_minimal_generator_presentation_shape','actual_relation_matrix_node','actual_relation_order','source_manifest'],
  'linear_rank':['complete_linear_rank_certificate'],
  'linear_source':['actual_selected37_original768_constant_combinations'],
  'zero620':['all620','normalized_Xp27_original768_combinations'],
  'quadratic_rank':['complete_fixed_rank_certificate'],
  'old39':['extra_quadratic_source_membership'],
  'quadratic464':['source_columns464'],
  'initial':['actual_alpha1_alpha2'],
  'three':['known_initial_scalars','remaining_one_scalar_obstruction'],
 }
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 prior=S.count(d)+20000;assert prior<150000
 assert d['operator']['actual_minimal_generator_presentation_shape']==[5,825]
 assert d['operator']['actual_relation_order']==['168 literal d3 images','657 remaining actual source cycles']
 l=d['linear_rank']['complete_linear_rank_certificate'];assert l['rank']==37
 chosen=l['selected_original_columns'];assert len(chosen)==37 and all(0<=j<657 for j in chosen)
 lifts=d['linear_source']['actual_selected37_original768_constant_combinations'];assert len(lifts)==37
 rank=d['quadratic_rank']['complete_fixed_rank_certificate'];assert rank['rank']==1 and rank['selected_original_columns']==[563]
 assert rank['basis_normalized_columns']==[[[0,'1'],[1,'11767/3652']]]
 assert d['old39']['extra_quadratic_source_membership'][1]['source_quad_index']==177
 zero=d['zero620']['all620'];assert len(zero)==620;selected=zero[563]
 assert selected['combination_in_actual768']==[[343,'100/3'],[541,'-16/3'],[643,'60'],[666,'-25/4'],[690,'1']]
 labels=d['quadratic464']['source_columns464'];assert len(labels)==464
 assert sum(x['kind']=='variable_times_actual_linear_relation'for x in labels)==296
 assert sum(x['kind']=='certified_actual_d3_quadratic'for x in labels)==168
 data={
  'scope':'SOURCE/CACHE INVENTORY only for one coherent scalar. No oldfiltered gamma/H, matrixM2/M3, new primitive jet or r-Jordan choice evaluated.',
  'actual_coefficient_names':['gamma1=q2 coefficient of rh1-alpha1*v','gamma2=q2 coefficient of rh2-rh3/4-alpha2*v','H=g coefficient of r²h3+c*q1'],
  'one_coherent_functional':'Delta=lambdaU*gamma2+(lambdaS*a3-lambdaU*aS)*gamma1-lambdaS*H-lambdaU*cS+lambdaS*c3',
  'same_original_minimal_operator':'OLD Told actual5x825 root29, first168 literal d3 pluslast657 sourcecycles. Not a348-dimensional scalar matrix, and notthe5x849 Tnil operator.',
  'actual_minimal_matrix_foreign_root':{'namespace':'OLD_TOLD_5_BY825','path':PATHS['operator'],'root':d['operator']['actual_relation_matrix_node']},
  'actual_selected37_linear_lift_bindings':[[j,168+j,combo]for j,combo in zip(chosen,lifts)],
  'binding_scope':'The saved37 original768 combinations are constant-source roots with full Schur/Euler transport in the accepted allorder operator. They are NOT independently expanded allorder coefficient tables.',
  'linear_lift_recipe':'Use fixed rank37 elimination receipts to lift L1=rh1-alpha1*v and L2=rh2-rh3/4-alpha2*v to constant825 coefficientvectors u01,u02 supportedon168+selected37. gamma_i=-pi_q2(M2*u0i). DoNOTdrop sourceunit/section corrections fornonzero-linear combinations.',
  'actual_old620_quadratic_pivot_namespace':'OLD_TOLD_ZERO_LINEAR620 index563 fromthisEXACTpath; do notidentify bynumericalindexalonewithotherT/T2 kernels.',
  'actual_old620_pivot_source_binding':selected,
  'pivot_projected_minor':rank['selected_minor_determinant'],
  'canonical_initial_G38':'h3*r²+(11767/3652)*h4*z5*r',
  'known177_lift_inventory':{'full_source_labels464':labels,'constant_rank177_basis':'quadratic_certificate in actual-linear-quadratic-T-certificate.json','source_count':{'linear37_times8':296,'literal_d3':168}},
  'H_lift_recipe':'Normalize the actual620 pivot byits saved-250/3 minor. ItsM2 differsfromcanonicalG38 bytheKNOWN177quadspace; retainandliftALL464-coefficientcorrections (native168 pluslinear37*vars). If resultingallorderliftis M(v0+v1+...), H=-pi_g(M3*v0+M2*v1), aftercorrectdegree2sourceclearings. HomogeneousG38 alone DOESNOTsetH0.',
  'source_availability_conclusion':'CompleteallorderOLDmatrix,37sourceconstantroots,620pivotrootandknown464source labels AREserialized. A separatelyexpandedfullallorderR37/177cofactorTaylor table isNOTknown. Thecoherentthree-adjointM2/M3liftmuststillbeconstructed andprimitivehighest/lowercache sensitivitiesmustbeproved beforeanycoefficientquery.',
  'unresolved_operational_gate':['exactlinear_lift37cofactors','canonicalG38allorder177correctionwithoutdroppingrawM2unknowns','onecombinedDeltaadjoint','primitiveorderandcacheavailability','positiveheight/owned/RAW/code/custody/freshpeer/ROOTONEGO'],
  'deferred_priority':'ParentnowprioritizesSIX_K_squared_product_SOURCE fordirectoriginalH0bridge. ThisDelta inventoryisfrozenwithoutanyvaluejob.',
  'all_three_coefficients_and_Delta_still_FORMAL':True,'no_originalH0_ell_torsion_nil_claim':True,
  'input_sha256':pins,'code_sha256':sha(__file__ if isinstance(__file__,Path)else Path(__file__)),'selective_receipts':receipts,
  'prior_source_owned_upper':prior,'all_passed':True}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'one-delta-three-old-filtered-coefficients-source-lift-and-cache-frontier.json';p.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_DELTA_THREE_COEFFICIENT_LIFT_INVENTORY_NO_RESPONSES',
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{p.name:sha(p)},'prior_source_owned_upper':prior,'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()

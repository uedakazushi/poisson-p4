#!/usr/bin/env python3
"""Finite source/operator availability; no comparison matrix evaluation."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CBAR_COMPARISON_DIVISION_PREFLIGHT_20261104_v63'
PATHS={
 'Cbar_resolution':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/actual-Cbar-minimal-resolution-and-three-targets.json',
 'Artin105':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-source105-Artin20-values.json',
 'Artin_peer':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW/report.json',
 'chainlift':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'chainlift_peer':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW/report.json',
 'complete_Q_operator':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
 'T_preflight':'reports/AMBIENT_B18_SPECIALIZED_COMPARISON_T_PREFLIGHT_20261104_v63/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 d={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert d['Artin_peer']['all_passed'] and d['chainlift_peer']['all_passed'] and d['T_preflight']['all_passed']
 assert d['Cbar_resolution']['ranks']==[3,24,84,168,210,168,84,24,3]
 assert d['Artin105']['all105_equations_exactly_zero']
 assert d['chainlift']['all640_UP_equals_QB_verified']
 E={
  'actual_parameter_frame':'Rbar=QQ[[zA0,...,zA6,rA]], ell=5(sA-rA)',
  'actual_targets3':['tauA','tauB','tauA*tauB'],
  'sector_annihilators':[['z7','r3'],['z7+c*r2','r3'],['z7','r2']],
  'actual_c_coefficients':d['Cbar_resolution']['CB_shift_coefficients7'],
  'canonical_m_operator':{
   'input':'each actual coefficient f_s of V3x40 lies in its sector regular ideal',
   'B_sector_coordinate_change':'u_i=z_i+c_i*r2; other sectors u_i=z_i',
   'ordered_hadamard':'h_i=(f(u0,...,ui,0after,r)-f(u0,...,u_(i-1),0after,r))/u_i',
   'last':'h_7=f(0,r)/r^n, n=3,3,2; exact because f is zero in the actual cyclic quotient',
   'regularity':'coordinate Hadamard and monomial Weierstrass division of the actual regular function; no fresh evaluation of an old r/g quotient',
   'shape':[24,40],'regular_coefficient_slots':960,
   'actual_m_coefficients_computed':False},
  'constant_unit_minor_source':{
   'linear_z_slots':'coefficient of z_j in the source A20 value',
   'last_slots':'coefficient r^n after z=0 or z=-c*r2, n=3,3,2',
   'curve_restriction_ring':'H20 -> QQ[r]/r5',
   'why_well_defined':'on either curve IA*IB has order at least5; exact normalized pair is(z7,r3),(z7+c*r2,r3)',
   'known_full_V56':'three nonconstant rows of the source Q48 and negative eight C0 odd columns',
   'Q48_values':'construct using existing Artin105 delta,bA,bB,nuA,nuB and the fixed16 tensor products',
   'U8_values_needed':'z7 exactly; gA_plane starts(125/24)r3, so its last division constant and zero linear slots are already cached',
   'full_m0_shape':[24,56],'full_m0_scalar_slots':1344,
   'proper_V40_extra_interface':'apply the accepted Res4 normalizer source-unit adapter; its scalar coefficient center/restriction must actually be bound before claiming the40column numerical matrix',
   'm0_numeric_matrix_or_pivot_selection_computed':False},
  'surjective_m_and_exact_cancellation':{
   'reason':'imV=imd1 and d1 has24minimal generators imply m mod maximalideal is surjective',
   'selection':'deterministically select24columns with nonsingular QQ center matrix',
   'M':'selected24x24 block of actual m; determinant is an actual formal unit',
   'n':'selected-column embedding times M^-1',
   'identities':['d1*m=V','V*n=d1','m*n=I24'],
   'h':'zero exactly',
   'column_change':'after source permutation m=[M,N], use W=[[M^-1,-M^-1*N],[0,I16]], giving V*W=[d1,0]',
   'kernel_decomposition':'kerV is source-isomorphic to im(d2) direct_sum Rbar16',
   'T_complete_conditional_shape':[100,920],
   'relations':'d3 supplies168; the752 actual X cycles supply752',
   'actual_unit_minor_or_inverse_or_column_change_built':False},
  'canonical_Koszul_cycle_contraction':{
   'input':'a sector8-vector a with sum g_i*a_i=0',
   'last_coordinate_step':'modulo seven coordinate generators, r^n*a_7(0,r)=0 in QQ[[r]], hence a_7(0,r)=0; split a_7 by ordered coordinate Hadamards',
   'boundary_step':'use the seven (i,7) wedges to remove a_7; remaining relation lies in the ordinary seven-coordinate Koszul complex',
   'recursion':'the standard ordered-coordinate Koszul contraction gives its explicit wedge2 preimage',
   'input_cycles':'m*X_j, since d1*m*X_j=V*X_j=0',
   'actual_X_transformed':'(m*X_j, X_unselected_rows_j)',
   'readback':'(KoszulLift(m*X_j), X_unselected_rows_j) in84+16 generators',
   'no_Gaussian_or_SyzP_needed':True,
   'all752_actual_cycle_coefficient_lifts_computed':False},
  'fallback_full56_source_route':{
   'before_proper_redundancy_removal':'Vfull3x56 and Xfull56x768 use fullQ48/U9 and completeQ128',
   'after24annunitcancellations':'conditionalT116/936',
   'then':'the16 known free-zero source syzygy units from the actual Res4 adapter yield100/920',
   'not_a_new_source_solve':True},
  'remaining_exact_gates':[
   'serialize full or proper V source functions in the actual three targets',
   'compute source-bound m0 and its fixed24column unit minor; no generic-rank extrapolation',
   'serialize canonical regular Hadamard/monomial-division operators and the formal unit inverse for M',
   'apply the source-bound Koszul contractions to the actual752 or768 cycles',
   'only then restrict the complete honest T presentation to the actual348ring',
   'old r/g are original functions under regular composition; never divide their specialized numerator by their specialized denominator'],
  'no_m_n_T_K_or_depth_computation':True,
  'input_hashes':{k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}}
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=True)
 ep=out/'actual-comparison-Hadamard-and-unit-cancellation-preflight.json';ep.write_text(json.dumps(E,indent=2)+'\n')
 R={'all_passed':True,'source_and_operator_preflight_only':True,
  'conditional_complete_T_shape':[100,920],'h_can_be_zero_after_actual_m_unit_selection':True,
  'full56_m0_existing_source_availability':True,'proper40_adapter_evaluation_gate_retained':True,
  'actual_m0_or_m_n_T_kernel_torsion_computed':False,'input_hashes':E['input_hashes'],
  'evidence':{'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)},'code_sha256':sha(Path(__file__))}
 (out/'report.json').write_text(json.dumps(R,indent=2)+'\n');print(json.dumps(R))
if __name__=='__main__':main()

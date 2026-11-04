#!/usr/bin/env python3
"""Source/theorem preflight only; no T matrix, kernel or torsion calculation."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_SPECIALIZED_COMPARISON_T_PREFLIGHT_20261104_v63'
GATES={
 'actual_chainlift':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/report.json',
 'chainlift_peer':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW/report.json',
 'complete_Q_peer':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63_REVIEW/report.json',
 'literal_Cbar_resolution':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/report.json',
 'actual_Cbar_peer':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW/report.json',
 'generic_C_rank_peer':'reports/AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW/report.json',
 'conductor_H20_peer':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 g={k:json.loads((ROOT/p).read_text()) for k,p in GATES.items()}
 assert all(d['all_passed'] for d in g.values())
 assert g['literal_Cbar_resolution']['ranks']==[3,24,84,168,210,168,84,24,3]
 assert g['generic_C_rank_peer']['actual_comparison_C_generic_DVR_rank']==4
 assert g['generic_C_rank_peer']['actual_mu_D_C']==8
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=True)
 evidence={
 'rings':{'D':'QQ[[ell]]','R':'actual regular R0','Rbar':'R/(ell)=QQ[[zA0,...,zA6,rA]]','I':'image(IA*IB) in Rbar'},
 'source_low_complex':{'V_shape':[3,40],'X_shape':[40,752],
  'construction':'apply accepted actual C-valued12->regular8 plus zero4 unit source adapter to Q48 and B48; discard16 zero source columns and use complete112 Koszul syzygies; form (Q_reduced,-U), then cancel its actual unit1 source/target pair',
  'sign':'V before target cancellation=(Q_reduced,-U); X includes (B_reduced,P) and (SyzQ_reduced,0)',
  'source_P':[9,640],'source_B':[48,640],'Q_reduced':[4,32],'SyzQ_reduced':[32,112],
  'actual_X752_and_V40_scalar_matrices_serialized_here':False},
 'definitions':{'J':'im(V)','L':'ker(V)','K':'L/im(X)','C':'coker(V)',
  'T':'ker(V mod ell)/im(X mod ell)'},
 'exact_specialization':{
  'J_ell_torsionfree':'J is a submodule of the free R^3, so multiplication by ell is injective',
  'first':'0 -> L/ellL -> ker(Vbar) -> C[ell] -> 0',
  'second':'0 -> K/ellK -> T -> C[ell] -> 0',
  'quotient_proof':'im(Xbar) is the image of the actual specialized columns inside L/ellL; right-exact cokernel specialization gives (L/ellL)/im(Xbar)=K/ellK',
  'no_naive_specialized_kernel_equals_Kbar':True},
 'actual_Cbar_presentation':{
  'targets':['tauA','tauB','tauA*tauB'],
  'Cbar_basis_lifts':[['e0','r*e0','r2*e0'],['e1','r*e1','r2*e1'],['e2','r*e2']],
  'actual_minimal_resolution_ranks':[3,24,84,168,210,168,84,24,3],
  'differentials':'three actual triangular regular-sequence Koszul complexes; literal source-bound3072 polynomial entries saved in the companion',
  'generic_DVR_resolution_needed':False},
 'complete_kernel_adapter':{
  'm_shape':[24,40],'n_shape':[40,24],'h_shape':[84,24],
  'required_actual_identities':['d1*m=Vbar','Vbar*n=d1','d2*h=I24-m*n'],
  'Lspecial_generators':'S=[n*d2, I40-n*m]',
  'S_shape':[40,124],
  'complete_relation_matrix':'[(d3,0),(-h,n)]',
  'relation_shape':[124,192],
  'completeness_proof':'if S(w,x)=0, y=m*x-d2*w satisfies x=n*y and d2*(w+h*y)=0; use ker(d2)=im(d3)',
  'lift_each_known_X_cycle':'for v=Xbar_j choose w with d2*w=m*v, then S(w,v)=v',
  'T_cokernel_shape':[124,944],
  'complete_SyzP_required':False,
  'actual_m_n_h_and752_cycle_lifts_constructed':False},
 'finite_support':{
  'premise':'the actual normalization tensor comparison has (IA*IB)K=(IA*IB)C=0',
  'I_annihilates_Kbar_and_C_ell':True,'I_squared_annihilates_T':True,
  'm5_subset_I':'accepted exactH20 has maximalideal fifthpower0',
  'm10_T_zero':True,'I_annihilates_T_claimed':False,
  'H20_naive_kernel_replacement_permitted':False,
  'correct_finite_base':'Rbar/I^2, not Rbar/I',
  'ring348_optional_source':'AMBIENT_B18_CONDUCTOR_SQUARE_CLOSED348_20261104_v63; its separate source/peer gate must be retained when using its concrete348actions'},
 'DVR_count':{
  'mu_D_C':8,'generic_DVR_rank_C':4,'dimension_QQ_C_ell':4,
  'formula':'mu_D(K)=length_QQ(T)-4',
  'why':'a finite DVR module of minimal generator count8 and rank4 has four torsion summands, and each contributes one dimension to C[ell]',
  'C_torsion_total_length_computed':False,
  'K_torsion_boolean':'nonzero iff mu_D(K)>rank_D(K); actual generic rankK still required',
  'T_length_K_mu_K_rank_K_torsion_depth_computed':False},
 'remaining_finite_completion_gates':[
  'serialize actual Vbar and Xbar through the already source-bound Q/B adapters',
  'construct and verify actual comparison m,n,h and752 cycle lifts using the complete literalCbar resolution',
  'restrict the honest complete T cokernel to actual Rbar/I^2; right-exactness is legal only after that presentation is complete',
  'evaluate all coefficient function references by their actual regular restriction operators; never redivide old r/g after specialization',
  'use adapted348 basis/cancellations and sparse operators before any large rational matrix',
  'compare the resulting mu_D(K) with independently established rank_D(K) for a torsion/depth conclusion'],
 'conservative_raw_348_counts':{'QQ_generators':124*348,'QQ_relation_columns':944*348,
  'dense_matrix_authorized':False,'full348_T_matrix_built':False},
 'source_hashes':{k:{'path':p,'sha256':sha(ROOT/p)} for k,p in GATES.items()}}
 ep=out/'actual-specialized-T-exact-sequence-and-comparison-preflight.json'
 ep.write_text(json.dumps(evidence,indent=2)+'\n')
 report={'all_passed':True,'proof_and_source_preflight_only':True,
  'exact_sequence':'0->K/ellK->T->C[ell]->0','C_ell_QQ_dimension':4,
  'actual_Cbar_literal_minimal_resolution_bound':True,'complete_T_conditional_cokernel_shape':[124,944],
  'I_squared_T_zero_and_m10_T_zero':True,'actual_T_or_K_length_rank_torsion_computed':False,
  'source_hashes':evidence['source_hashes'],'code_sha256':sha(Path(__file__)),
  'evidence':{'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)}}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()

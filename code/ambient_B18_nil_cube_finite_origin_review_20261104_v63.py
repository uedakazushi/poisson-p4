#!/usr/bin/env python3
"""Independent dependency/scope and pure-mathematical audit; no CAS job."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63_REVIEW'
PRIMARY='reports/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63/report.json'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for x in iter(lambda:f.read(1<<20),b''):h.update(x)
 return h.hexdigest()
def main():
 p=json.loads((ROOT/PRIMARY).read_text());assert p['all_passed']
 checked={}
 for path,want in p['input_sha256'].items():
  assert sha(ROOT/path)==want,path
  checked[path]=want
  if path.endswith('/report.json'):
   assert json.loads((ROOT/path).read_text())['all_passed'],path
 rank=json.loads((ROOT/'reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json').read_text())
 assert rank['rank_D_K']==4 and rank['rank_D_C']==4
 assert rank['generic_kernel_branch_lengths']=={'CA_critical_V':1,'CA_selfpair_Pa':1,'CB_critical_V':0,'CB_selfpair_Qb':2}
 core=json.loads((ROOT/'reports/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW/report.json').read_text())
 assert core['core_and_normalized_E4_exact_nilradical_index3_accepted']
 tail=json.loads((ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json').read_text())
 assert tail['tail_nil_principal_square_zero_and_cyclic_same_R0_regular_curve']
 assert tail['Dalg_polynomial_Laurent_scope_correction_required_and_bound']
 scope=json.loads((ROOT/'reports/AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63.json').read_text())
 assert scope['all_passed'] and scope['frozen_w0_completed_ring_not_dehomogenized']
 assert p['all_positive_homogeneous_coefficient_maps_A_n_to_M_injective']
 assert not p['global_positive_comparison_kernel_vanishing_assumed']
 assert not p['nil_cube_zero_decided'] and not p['actual_length_T_or_origin_depth_decided']
 assert not p['nil_cube_zero_implies_depth20_claimed']
 proof={
  'homogeneous_nil_projection':{
   'argument':'the nilradical is a closed ideal stable under the continuous positive torus action; finite Vandermonde weight projections modulo positive-degreeN isolate each component up to a tail in maximalidealN; closedness then puts every component in Nil(A)',
   'completion':'Noetherian completed ring; nilradical finitely generated and closed; weight0 coefficients already complete; no w=1 substitution in the completed w0 ring'},
  'mixed_products':{
   'argument':'for nil homogeneous a_i,b_j, both images in the actual stabletail are nil; tailNil squared0; injectivity in degree i+j>=1 forces product0',
   'includes_degree0_positive_mixed_products':True,
   'conclusion':'Nil(A)^2=Nil(A0)^2 and Nil(A)^3=Nil(A0)^3',
   'positive_annihilator':'for positive c, b0*c is positive homogeneous nil, so a0*(b0*c)=0; pass coefficientwise to completed series'},
  'kernel_and_finiteness':{
   'actual_kernel':'K=ker(D0->E4) for the same physical source, with integral Spec-surjective comparison',
   'cube_in_kernel':'Nil(E4)^3=0 sends the degree0 nilcube to0',
   'finite_D':'actual IA*IB annK, Noetherell regular on the four actual conductor branches',
   'generic_zero':'use whole actual local D rings on each generic branch, including their nilthickenings; not merely reduced point values',
   'conclusion':'Nil(A)^3 is finite length and lies in tors_DK=H0(D0)'},
  'conditional_index':{
   'exact_K_rank_gate':4,'mu_DK_formula':'length(T)-4',
   'lengthT8':'mu_DK=rank_DK=4 implies DVRtorsionfree, so cube0',
   'square_nonzero':'the actual Qb family specializes toP0 and has whole nilsquare nonzero generically; coherent support closedness and faithful completion preserve nonzero square atP0',
   'conclusion':'lengthT8 implies exactnilindex3; cubeNonzero implies depth19; cube0 does not imply depth20'},
  'group_restoration':{
   'module':'QQ[[G19]] completed_tensor_QQ finite_length_relative_cube',
   'finite_free_rank':'relative cube length, uncomputed',
   'support_bound_cone_projective':[19,18],
   'no_global_associated_closure_or_component_claim':True},
  'subject':'actual CA intersect CB, not whole Jacobi scheme',
  'no_new_scalar_source_jet_point_rank_or_GB':True,
  'primary_report':{'path':PRIMARY,'sha256':sha(ROOT/PRIMARY)},'checked_input_hashes':checked}
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=True)
 ep=out/'independent-nil-cube-proof-and-gates.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
 rep={'all_passed':True,'independent_pureproof_and_all_dependency_hashes':True,
  'full_completed_nil_cube_equals_degree0_cube':True,'mixed_positive_products_zero':True,
  'nil_cube_finite_length_relative_origin':True,'conditional_lengthT8_exactnilindex3':True,
  'rank_K4_distinguished_from_rank_C4':True,
  'actual_nilcube_zero_Tlength_depth_or_globalclosure_decided':False,
  'primary_report_sha256':sha(ROOT/PRIMARY),'input_hashes':checked,
  'evidence':{'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)},
  'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md'))}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':main()

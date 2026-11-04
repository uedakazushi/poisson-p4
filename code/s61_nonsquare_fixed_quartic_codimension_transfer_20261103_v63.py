#!/usr/bin/env python3
"""Static premise/hash audit for the conceptual fixed-cut codimension transfer.

No coefficients, matrices, points, fields or Groebner bases are recomputed.
"""
from pathlib import Path
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_FIXED_QUARTIC_CODIMENSION_TRANSFER_20261103_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();names=['S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_20261102','S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_REVIEW_20261102','S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102','S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102','S61_BINARY_HYPERPLANE_CUBIC_CODIMENSION_20261102','S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102'];paths=[ROOT/'reports'/(n+'.json')for n in names];data=[json.loads(p.read_text())for p in paths]
 for p,d in zip(paths,data):
  assert d['all_passed']
  for rel,want in d.get('input_sha256',{}).items():assert sha(ROOT/rel)==want,rel
 tr,rv,ch,q,c,g=data;assert tr['exceptional_core_dimension_upper']==1 and tr['coordinate_free_existence_of_transverse']and tr['relative_module_rank_on_zero_nuisance_section']==65
 assert rv['proof_audit_accepted']and rv['completed_regular_base_dimension']==11 and not rv['H11_in_reviewed_scope']
 assert ch['family_parameter_minor']=='-(100*a + 18*e**2 - 9*e - 2)/25'and ch['family_parameters']==[1,2,5]and ch['normal_weights34']==[1,1,1,1,1,2,1,1]
 assert ch['selected_right_inverse_checked']and ch['full_linear_kernel_checked']and ch['implicit_solution_exact_degree_bound']==4
 assert q['bad_locus_dimension_upper_on_three_core_open']==1 and c['cubic_rank_drop_dimension_upper']==1 and c['all_core_cubic_rank_upper']==57
 assert 1+7+(29-6)+(91-57)==65 and ch['weighted_column_counts']['4']==239
 note=ROOT/'research/notes'/(ID+'.md');paths.extend([Path(__file__),note,ROOT/'research/notes/S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_20261102.md',ROOT/'research/notes/S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_REVIEW_20261102.md',ROOT/'research/notes/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.md',ROOT/'code/s61_binary_hyperplane_weighted_chart_20261102.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py'])
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);report={'id':ID,'all_static_premise_hash_checks_passed':True,'new_mathematical_proof_in_note':True,'parameterchange_uses_nil_thickening_CM_and_intrinsic_generic_length65':True,'fixed_full_quartic_defect_codimension_at_least2_on_aFBC':True,'alternate_C_cut_not_identified_with_fixed_row_space':True,'C_zero_outside_stated_base':True,'G6_and_true135_retained':True,'no_new_coefficient_matrix_point_field_or_GB':True,'full239_all_points_NOT_proved':True,'eighty_pointwise_independent_source_kernel_directions_still_missing':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'seconds':time.monotonic()-start,'frozen':True,'independent_conceptual_review_pending':True};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2))
if __name__=='__main__':main()

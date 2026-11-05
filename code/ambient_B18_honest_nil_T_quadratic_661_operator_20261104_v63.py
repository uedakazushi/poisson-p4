#!/usr/bin/env python3
"""Complete6x661 residual quadratic operator, with no new coefficient query."""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import mv,enc,sparse,nqq
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_QUADRATIC_661_OPERATOR_20261104_v63'
PATHS={
 'first':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'first_peer':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63_REVIEW/report.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'quadratic':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json',
 'd3_source':'reports/AMBIENT_B18_HONEST_NIL_T_D3_QUADRATIC_SOURCE_PREFLIGHT_20261104_v63/actual153-d3-quadratic-adjoints-and12-callback-sensitivities.json',
 'd3_values':'reports/AMBIENT_B18_HONEST_NIL_T_D3_QUADRATIC_VALUES_20261104_v63/actual153-d3-quadratic-quotient-values-and-known-rank.json',
 'd3_report':'reports/AMBIENT_B18_HONEST_NIL_T_D3_QUADRATIC_VALUES_20261104_v63/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 D={k:read(ROOT/p)for k,p in PATHS.items()}
 assert D['first_peer']['all_passed']and D['d3_report']['all_passed']and D['d3_report']['actual_d3_quadratic_projection_rank']==0
 A=D['first'];lp=A['rank_minor_columns'];lr=A['rank_minor_rows'];Li=[[Q(q)for q in v]for v in A['rank_minor_inverse']]
 assert len(lp)==35 and all(j<681 for j in lp)
 cols=[]
 for sp in A['all849_actual_minimal_first_columns40_sparse']:
  v=[Q(0)]*40
  for i,c in sp:v[i]=Q(c)
  cols.append(v)
 def lift(j):
  a=mv(Li,[cols[j][i]for i in lr]);v={j:Q(1)}
  for k,c in enumerate(a):
   if c:v[lp[k]]=v.get(lp[k],Q(0))-c
  assert all(sum((cols[k][i]*c for k,c in v.items()),Q(0))==0 for i in range(40))
  return [[i,c]for i,c in sorted(v.items())if c]
 pure=[j for j in range(681)if j not in lp];assert len(pure)==646
 dp=D['minimal']['d3_only_linear']['pivot_columns'];assert len(dp)==15
 mixed=[681+j for j in dp];source=[lift(j)for j in pure+mixed]
 assert len(source)==661
 other=[[[681+i,Q(c)]for i,c in v]for v in D['d3_source']['all153_constant_d3_source_combinations']]
 assert len(other)==153
 # The block identity pivots prove all814 source vectors independent.
 for j,v in zip(pure+mixed,source):
  assert dict(v)[j]==1 and all(i in lp or i==j for i,c in v)
 dpfree=[i for i in range(168)if i not in dp]
 for j,v in zip(dpfree,other):
  assert dict(v)[681+j]==1 and all(i in [681+k for k in dp]or i==681+j for i,c in v)
  assert all(sum((cols[k][i]*c for k,c in v),Q(0))==0 for i in range(40))
 assert 646+15+153==849-35
 # Matrix objects are finite coefficient-extraction operators, not153k scalar
 # series slots; every input entry is the already pinned actual T_min function.
 operators=[
  {'operation':'actual_external_matrix','namespace':'minimal','field':'minimal_T5_by849','shape':[5,849]},
  {'operation':'entrywise_formal_homogeneous_degree_coefficient','input':0,'degree':2,'variables':A['variable_order'],
   'output_row_order':'target-major5 then saved36 unordered quadratic monomials','shape':[180,849]},
  {'operation':'literal_sparse_QQ_rows','rows':D['quadratic']['exact_quotient_duals_rows_sparse'],'shape':[6,180]},
  {'operation':'literal_sparse_QQ_columns','columns':source,'shape':[849,661]},
  {'operation':'matrix_product','inputs':[2,1,3],'shape':[6,661]}]
 data={'actual_variables':A['variable_order'],'actual_targets':A['target_order'],
  'known_quadratic_span_rank':174,'actual_d3_only153_quotient_values_zero':True,
  'all646_Xonly_zero_linear_source_labels':pure,'all15_mixed_source_labels_original849':mixed,
  'all661_source_columns_in_original849':source,'all153_complementary_d3_only_columns_in_original849':other,
  'all814_source_columns_first_relations_zero_checked':True,
  'complete_constant_kernel_proof':'The646 X-only vectors have unique free-X identity coordinates, the15 mixed vectors have unique d3-pivot identity coordinates, and the153 d3-only kernel vectors have unique free-d3 identity coordinates. All lie in ker M1, are independent by these triangular blocks, and total849−rank35=814.',
  'actual_operator_nodes':operators,'actual_residual_operator_shape':[6,661],
  'exact_formula_for_actual_HF2':'6−rank of the actual evaluated operator6x661 over QQ',
  'source_images_not_numerically_evaluated':True,'no_full180_by849_scalar_coefficients_enumerated':True,
  'frozen_input_interpretation':'minimal matrix operator node is read from the actual all-order5x849 source-bound Schur representation; formal Taylor coefficients taken before any finite normal form or quotient',
  'smallest_current_constant_source_split_only_not_rank_minimal_claim':True,
  'remaining_source_order_gates':{
   'top2_quadratic_adjoints':'actual X2 plus first unit coefficients; mixed z/r directions selected by fixed six duals',
   'three_wedge_quadratic_readbacks':'sector0 row7 of (mX)3, requiring X3 row20 and m1row7 X2 + m2row7 X1 before any source-specific cancellation',
   'highest_formal_monomials_before_cancellation':'only z0*r²,z4*r²,z5*r² in row20 X3 for the wedge r² duals',
   'known_m1':'true C44 row7 data from FIRST_LINEAR, retained',
   'm2':'pure-r³ coefficient of Vgamma_A before dividing by r, not present in B11',
   'source_A0_common105_and_own57_high_orders':'existing caches are plane/univariate; they do not by themselves supply all z_i*r² source directions',
   'unknown_unit_pivot_sensitivities':'the81 pivot first/quadratic frame callbacks must be transported together; final681 first sensitivity proof does not establish their higher cancellation',
   'mandatory_before_any_new_values':'prove cancellation or bind every required callback to an existing cache and derive full saved+live operator budget',
   'new_X3_m2_source_gauge_jet_authorized_or_run':False},
  'no_actual_HF2_total_length_finite_GB_annihilator_or_nil_cube_answer':True}
 assert nqq(data)<150000
 p=out/'actual-complete-six-by661-quadratic-source-operator.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;p.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'PROVED_complete6x661_actual_quadratic_functional_with_scalar_values_open',
  'known_quadratic_rank':174,'residual_actual_operator_shape':[6,661],
  'exact_source_split':[646,15,153],'all814_sources_complete':True,'actual_HF2_proved':None,
  'new_X3_m2_source_gauge_jet_or_scalar_residual_values':False,'visible_QQ_support':nqq(data),'bytes':p.stat().st_size,
  'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()

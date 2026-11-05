#!/usr/bin/env python3
"""Complete actual cubic3x917 SOURCE operator. No M2/M3 value evaluation."""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_CUBIC917_SOURCE_OPERATOR_20261105_v63'
PATHS={
 'formal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'first':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'quad':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json',
 'C661':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_661_OPERATOR_20261104_v63/actual-complete-six-by661-quadratic-source-operator.json',
 'C153':'reports/AMBIENT_B18_HONEST_NIL_T_D3_QUADRATIC_SOURCE_PREFLIGHT_20261104_v63/actual153-d3-quadratic-adjoints-and12-callback-sensitivities.json',
 'quad_values':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_VALUES_20261105_v63/phase3-all661-actual-six-quadratic-values-and-rank.json',
 'quad_peer':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_REMAINING_20261105_v63_REVIEW/report.json',
 'cubic_input':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_PREFLIGHT_20261105_v63/literal-known-cubic-input-and-exact24-row-reduction.json',
 'cubic_rank':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_20261105_v63/fixed-known-cubic-rank-and-length-completion.json',
 'cubic_peer':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_20261105_v63_REVIEW/report.json',
 'joint_report':'reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63/report.json',
 'joint_rhs':'reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63/phase2-known-third-RHS-and-proper-response.json',
 'joint_values':'reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63/phase4-selected-joint135-third-adjoint-values.json',
 'joint_peer':'reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def put(d,k,c):
 if c:
  d[k]=d.get(k,Q(0))+c
  if not d[k]:del d[k]
def add(d,v,c=Q(1)):
 d=dict(d)
 for k,a in v.items():put(d,k,c*a)
 return d
def rec(sp):return {i:Q(c)for i,c in sp if Q(c)}
def sp(d):return [[i,c]for i,c in sorted(d.items())if c]
def count(x):
 if isinstance(x,Q):return 1
 if isinstance(x,dict):return sum(count(v)for v in x.values())
 if isinstance(x,(list,tuple)):return sum(count(v)for v in x)
 return 0
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};assert all(D[k]['all_passed']for k in ['quad_peer','cubic_peer','joint_report','joint_peer'])
 assert D['formal']['matrices']['minimal_T5_by849']==51
 first=[rec(v)for v in D['first']['all849_actual_minimal_first_columns40_sparse']];R35=D['first']['rank_minor_columns'];assert len(R35)==35
 C=[rec(v)for v in D['C661']['all661_source_columns_in_original849']]+[{681+i:Q(c)for i,c in v}for v in D['C153']['all153_constant_d3_source_combinations']]
 assert len(C)==814
 for col in C:
  z={}
  for j,c in col.items():z=add(z,first[j],c)
  assert not z
 # Constant coefficients surviving the full quadratic constraint.
 V=D['quad_values'];cp=V['fixed_rank_columns'];assert len(cp)==3
 kernel811=[];labels=[]
 for j in range(814):
  if j in cp:continue
  f=[Q(c)for c in V['all661_rank_factorizations'][j]]if j<661 else[Q(0)]*3
  v={j:Q(1)}
  for p,c in zip(cp,f):put(v,p,-c)
  assert all(sum((Q(V['all661_actual_quadratic_six_values'][k][i])*c for k,c in v.items()if k<661),Q(0))==0 for i in range(6))
  kernel811.append(v);labels.append(j)
 assert len(kernel811)==811 and sum(j<661 for j in labels)==658
 # All106 syzygies of L=(M1R35)*S1, including the mixed-variable ones.
 A=D['quad'];L=[rec(v)for v in A['all280_known_actual_variable_multiple_columns180_sparse']];lp=A['fixed_rank_pivot_columns'];lr=A['fixed_rank_pivot_rows'];iv=[[Q(c)for c in row]for row in A['fixed_minor_inverse']];assert len(lp)==174
 ivsp=[[(j,c)for j,c in enumerate(row)if c]for row in iv];syzygies=[];slabels=[]
 for j,col in enumerate(L):
  if j in lp:continue
  rhs={i:col[r]for i,r in enumerate(lr)if r in col};f=[sum((c*rhs.get(i,Q(0))for i,c in row),Q(0))for row in ivsp];v={j:Q(1)}
  for p,c in zip(lp,f):put(v,p,-c)
  check={}
  for k,c in v.items():check=add(check,L[k],c)
  assert not check
  syzygies.append(v);slabels.append(j)
 assert len(syzygies)==106
 # Three actual cubic NF quotient duals, preserving integer row scaling.
 B,K=D['cubic_input'],D['cubic_rank'];p=K['fixed_pivot_columns'];r=K['fixed_pivot_rows'];free=[i for i in range(24)if i not in r];assert free==[7,15,23]
 ints=[rec(v)for v in B['equivalent_integer24_by1308_sparse_columns']];inv=[[Q(c)for c in row]for row in K['fixed_minor_inverse']];den=list(map(Q,B['row_denominator_units24']));qcols=[rec(v)for v in B['canonical_q24_by600_columns_sparse']]
 pr=[]
 for f in free:
  v={f:Q(1)}
  for a,rr in enumerate(r):put(v,rr,-sum((ints[p[j]].get(f,Q(0))*inv[j][a]for j in range(21)),Q(0)))
  pr.append({j:c*den[j]/den[f]for j,c in v.items()})
 delta=[]
 for row in pr:delta.append({j:sum((c*col.get(i,Q(0))for i,c in row.items()),Q(0))for j,col in enumerate(qcols)if sum((c*col.get(i,Q(0))for i,c in row.items()),Q(0))})
 for raw in B['all1284_known_cubic_columns600_sparse']:
  col=rec(raw);assert all(sum((c*col.get(i,Q(0))for i,c in row.items()),Q(0))==0 for row in delta)
 NF=[119,479,599];assert [[row.get(j,Q(0))for j in NF]for row in delta]==[[Q(i==j)for j in range(3)]for i in range(3)]
 # Finite source operators, factored; Taylor callbacks are exact named functions.
 nodes=[]
 def node(op,shape):i=len(nodes);nodes.append({'op':op,'shape':shape});return i
 M=node(['actual_all_order_formal_matrix',PATHS['formal'],'matrices.minimal_T5_by849',51],[5,849])
 M1=node(['homogeneous_coefficients_original_S',M,1,[5,8]],[40,849])
 M2=node(['homogeneous_coefficients_original_S',M,2,[5,36]],[180,849])
 M3=node(['homogeneous_coefficients_original_S',M,3,[5,120]],[600,849])
 RC=node(['constant_column_inclusion',R35,849],[849,35]);CC=node(['literal_sparse_columns',list(map(sp,C))],[849,814]);DD=node(['literal_sparse_columns',list(map(sp,kernel811))],[814,811])
 HH=node(['fixed_right_inverse_on_known174','embed rows',lp,'inverse external',PATHS['quad'],'fixed_minor_inverse','extract rows',lr],[280,180])
 KK=node(['literal_sparse_columns',list(map(sp,syzygies))],[280,106]);QQ=node(['literal_sparse_rows',list(map(sp,delta))],[3,600])
 CD=node(['matrix_product',CC,DD],[849,811]);Q2=node(['matrix_product',M2,CD],[180,811]);U1=node(['matrix_product',HH,Q2],[280,811]);MR2=node(['matrix_product',M2,RC],[180,35]);C3=node(['matrix_product',M3,CD],[600,811])
 correction=node(['multiply_degree2_relation_columns_by_degree1_coefficients',MR2,U1,'domain order:35 original selected relations then8 variables','ordinary monomial convolution, no divided powers'],[600,811])
 const=node(['matrix_product',QQ,['matrix_difference',C3,correction]],[3,811]);linear=node(['matrix_product',QQ,['multiply_degree2_relation_columns_by_degree1_coefficients',MR2,KK]],[3,106]);allop=node(['horizontal_blocks',const,linear],[3,917])
 jointmons=D['joint_rhs']['monomials5'];assert jointmons==[[7,7,7],[0,7,7],[1,7,7],[4,7,7],[5,7,7]]
 data={'actual_S_variables':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],'actual_minimal_relation_order849':D['formal']['retained_original_relation_columns849'],
  'complete814_constant_first_kernel_columns_in_original849':list(map(sp,C)),
  'all811_constant_quadratic_kernel_columns_in_C814':list(map(sp,kernel811)),'all811_source_labels_in_C814':labels,
  'constant_split':[658,153],'all814_first_images_and811_quadratic_NF_images_zero':True,
  'selected35_actual_first_relations':R35,'linear_coefficient_domain280_labels':A['all280_source_labels_relation_and_variable'],
  'all106_first_linear_syzygies280_sparse_columns':list(map(sp,syzygies)),'all106_free_source_labels':slabels,'all106_L_images_zero_verified':True,
  'canonical_H174_definition':'embed at174 pivot columns; multiply cached174x174 inverse; extract174 pivot rows. It is a right inverse on im L, not allF180.',
  'cubic_NF_dual3_by600_sparse':list(map(sp,delta)),'cubic_NF_representatives600':[119,479,599],'all1284_known_cubic_images_killed_and_r3_embedding_I3':True,
  'typed_operator_nodes':nodes,'exported_complete_cubic_operator':allop,'actual_complete_cubic_shape':[3,917],
  'all_order_Taylor_callback_definition':'exact ordinary homogeneous coefficients in the SAME physical CA-F9 ring S of the accepted all-order5x849 operator51, retaining all unit inverses/regular r,g quotients and mixed d3 paths; not scalar values filled by a recipe',
  'two_stage_correction':'constant D811 first kills the rank3 quadratic NF. For Q2=M2*C*D in im L, actual source degree1 correction is -R35*H174(Q2).',
  'completeness_proof':'For any formal source a=a0+a1+a2..., degree1zero gives a0=Cc. Degree2zero gives c=Dd and a1=-R H174(M2CD)d + R K106 e + C*v1. The C*v1 contribution atdegree3 lies known cubic closure because M2C modulolinear products spans exactly3 NF quadrics; M1*a2 also lies that closure. Thus the917 columns give every additional actual cubic initial relation.',
  'exact_actual_HF3_formula':'HF3(Tnil)=3-rank(actual3x917)',
  'exact_completion_target':{'rank0_HF':[5,5,3,3],'rank1_HF':[5,5,3,2,1],'rank2_HF':[5,5,3,1,1,1],'rank3_impossible':'length16 and nonzero remaining3'},
  'highest_source_sensitivity_operator':{'delta_M3':'Delta3 * deltaM3 * C * D','delta_M2_constant_block':'-Delta3*((deltaM2*R)*H174*(M2*C*D)+(M2*R)*H174*(deltaM2*C*D)) with degree2-times-degree1 convolution','delta_M2_linear_block':'Delta3*(deltaM2*R)*K106','source_warning':'deltaM2 anddeltaM3 must be pulled back COHERENTLY through the same physical/formal source; they are not arbitrary independent cached arrays'},
  'primitive_highest_chain_gate_open':True,
  'joint135_cache_scope':{'monomial_tuples':jointmons,'available':'152 selected adjoint values (including83 gauge), only requested implicit -rhoRHS contributions for old619 functionals; not a solved135-vector and not full native M3','base_name_caution':'joint B25/cover100 means base20 plus5 THIRD triples; gauge120 base25/cover100 means base20 plusfirst z6 plus4 SECOND pairs. Match tuples, never indices ordimension100.',
   'known_values_must_not_zero_fill_missing':'No direct binding from all M2/M3 Taylor callbacks above to the selected joint functionals has yet been proved.'},
  'source_order_frontier':{'naive_top_cubic':'raw X degree3 and all unit-frame carries','naive_wedge_cubic':'canonical H dividesonezi, so raw row20 X degree4 appears; no new coefficient is authorized','retained':'m3/m2 and last-sector source contractions plus higher unit corrections remain in exactTaylor callbacks until complete projected sensitivity proves cancellation','known_cancellations':'earlier raw B/P highest quotient-gradient identities are order-independent constant maps, but their lower mixed products and the106 linear-coefficient paths still require a separate coherent gate'},
  'actual_M2_M3_or_new_cubic_values_computed':False,'new_source_gauge_jet':False}
 assert count(data)<150000
 pth=out/'complete-cubic917-typed-source-and-three-NF-adjoints.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;pth.write_text(blob)
 r={'id':STEM,'all_passed':True,'status':'COMPLETE_SOURCE_OPERATOR_values_and_primitive_cache_gate_open','complete_actual_cubic_shape':[3,917],'constant_kernel_split':[658,153],'first_linear_syzygy_count':106,'all15_mixed_and153_d3_sources_preserved':True,'all_first_and_second_NF_kernel_checks':True,'actual_cubic_scalar_values_or_HF3_decided':False,'primitive_highest_source_sensitivities_evaluated':False,'new_source_gauge_jet':False,'typed_QQ_slots':count(data),'bytes':len(blob.encode()),'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{pth.name:sha(pth)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':main()

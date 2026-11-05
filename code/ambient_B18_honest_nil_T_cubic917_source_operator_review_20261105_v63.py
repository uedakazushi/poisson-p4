#!/usr/bin/env python3
"""Independent complete cubic-source/type theorem, values remain OPEN."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_HONEST_NIL_T_CUBIC917_SOURCE_OPERATOR_20261105_v63';ID=BASE+'_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),indent=2)+'\n')
 rep=load(BASE);assert pins[f'reports/{BASE}/report.json']=='d2c9b13ef33b9cd18d549b9a0fe7e467fb4d3bdc19323ae8239848f6a141bc8b'
 for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h
 for n,h in rep['evidence_sha256'].items():assert sha(ROOT/'reports'/BASE/n)==h
 D=load(BASE,'complete-cubic917-typed-source-and-three-NF-adjoints.json');Min=load('AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63','actual-minimal-Tnil5-849-and-d3-linear-operator.json')
 first=load('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase3-all849-minimal-first-relations-and-rank.json');B=load('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63','known-quadratic-closure-and-complete814-source-projection.json')
 C661=load('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_661_OPERATOR_20261104_v63','actual-complete-six-by661-quadratic-source-operator.json');C153=load('AMBIENT_B18_HONEST_NIL_T_D3_QUADRATIC_SOURCE_PREFLIGHT_20261104_v63','actual153-d3-quadratic-adjoints-and12-callback-sensitivities.json')
 assert load('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_REMAINING_20261105_v63_REVIEW')['actual_complete6x661_rank']==3
 V=load('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_VALUES_20261105_v63','phase3-all661-actual-six-quadratic-values-and-rank.json')
 assert load('AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_20261105_v63_REVIEW')['actual_m3_equals_rA_m2']
 cub=load('AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_20261105_v63_REVIEW','independent-homogeneous600-closure-and-r-reduction-proof.json')
 length=load('AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63');assert length['actual_Tnil_length']==16
 assert Min['matrices']['minimal_T5_by849']==51 and D['actual_minimal_relation_order849']==Min['retained_original_relation_columns849']
 save('phase.json',dict(phase='all actual all-order/first/quadratic/known-cubic peers pinned; source-only'))
 M1=[Q.pairs(v)for v in first['all849_actual_minimal_first_columns40_sparse']];R35=D['selected35_actual_first_relations'];assert R35==first['rank_minor_columns']and max(R35)<681
 C=[Q.pairs(v)for v in D['complete814_constant_first_kernel_columns_in_original849']]
 assert C==[Q.pairs(v)for v in C661['all661_source_columns_in_original849']]+[{681+i:F(c)for i,c in v}for v in C153['all153_constant_d3_source_combinations']]
 assert all(not Q.plus(*(Q.scaled(M1[j],a)for j,a in col.items()))for col in C)
 # Explicit unit triangular free-coordinate witness, not a column count.
 fx=[j for j in range(681)if j not in R35];assert len(fx)==646
 d3free=C153['all153_zero_linear_d3_source_labels'];d3pivot=[j for j in range(168)if j not in d3free];assert len(d3pivot)==15
 xleaders=[]
 for col in C[:646]:
  sp={j:a for j,a in col.items()if j in fx};assert len(sp)==1 and next(iter(sp.values()))==1;xleaders.append(next(iter(sp)))
 assert set(xleaders)==set(fx)and len(set(xleaders))==646
 mixed=[]
 for col in C[646:661]:
  sp={j-681:a for j,a in col.items()if j>=681};assert len(sp)==1 and next(iter(sp.values()))==1;mixed.append(next(iter(sp)))
 assert set(mixed)==set(d3pivot)and len(set(mixed))==15
 for label,col in zip(d3free,C[661:]):assert {j-681:a for j,a in col.items()if j-681 in d3free}=={label:F(1)}and all(j>=681 for j in col)
 # In order fx646, pivot-d3 15, free-d3 153, C has block-upper-triangular
 # unit diagonal. Thus its814 independent columns are WHOLE ker(M1).
 K811=[Q.pairs(v)for v in D['all811_constant_quadratic_kernel_columns_in_C814']];labels=D['all811_source_labels_in_C814'];cp=V['fixed_rank_columns'];assert labels==[j for j in range(814)if j not in cp]
 qvals=[Q.sparse(v)for v in V['all661_actual_quadratic_six_values']]+[{}for _ in range(153)]
 for j,col in zip(labels,K811):
  assert {i:a for i,a in col.items()if i not in cp}=={j:F(1)}
  assert not Q.plus(*(Q.scaled(qvals[i],a)for i,a in col.items()))
 assert [sum(j<661 for j in labels),sum(j>=661 for j in labels)]==[658,153]
 assert all(j in labels for j in range(646,661)) # ALL15 mixed paths
 L=[Q.pairs(v)for v in B['all280_known_actual_variable_multiple_columns180_sparse']];lp=B['fixed_rank_pivot_columns'];lr=B['fixed_rank_pivot_rows'];inv=[Q.sparse(v)for v in B['fixed_minor_inverse']]
 assert D['linear_coefficient_domain280_labels']==B['all280_source_labels_relation_and_variable']
 expectedlabels=[[j,v]for j in R35 for v in range(8)];assert B['all280_source_labels_relation_and_variable']==expectedlabels
 H=[]
 for col in L:
  a=Q.mvec(inv,{i:col[row]for i,row in enumerate(lr)if col.get(row,F(0))});pre={lp[i]:x for i,x in a.items()};assert Q.plus(*(Q.scaled(L[j],x)for j,x in pre.items()))==col;H.append(pre)
 K106=[Q.pairs(v)for v in D['all106_first_linear_syzygies280_sparse_columns']];free280=[j for j in range(280)if j not in lp];assert D['all106_free_source_labels']==free280
 for j,col in zip(free280,K106):
  assert {i:a for i,a in col.items()if i not in lp}=={j:F(1)}
  assert not Q.plus(*(Q.scaled(L[i],a)for i,a in col.items()))
 assert len(K106)==106
 Delta=[Q.pairs(v)for v in D['cubic_NF_dual3_by600_sparse']];assert Delta==[Q.pairs(row.items())for row in cub['known_cubic_quotient3_by600']]
 assert D['cubic_NF_representatives600']==[119,479,599]and [[row.get(j,F(0))for j in[119,479,599]]for row in Delta]==[[F(i==j)for j in range(3)]for i in range(3)]
 # Verify each finite typed node against the actual matrix source and all
 # literal scalar blocks; no Taylor coefficient callback is evaluated.
 N=D['typed_operator_nodes'];assert len(N)==19
 expectedshapes=[[5,849],[40,849],[180,849],[600,849],[849,35],[849,814],[814,811],[280,180],[280,106],[3,600],[849,811],[180,811],[280,811],[180,35],[600,811],[600,811],[3,811],[3,106],[3,917]]
 assert [v['shape']for v in N]==expectedshapes
 ops=[v['op']for v in N];assert ops[0][0]=='actual_all_order_formal_matrix'and ops[0][3]==51
 for j,d in[(1,1),(2,2),(3,3)]:assert ops[j]==['homogeneous_coefficients_original_S',0,d,[5,[8,36,120][d-1]]]
 assert ops[4]==['constant_column_inclusion',R35,849]
 assert ops[5]==['literal_sparse_columns',D['complete814_constant_first_kernel_columns_in_original849']]
 assert ops[6]==['literal_sparse_columns',D['all811_constant_quadratic_kernel_columns_in_C814']]
 assert ops[7][0]=='fixed_right_inverse_on_known174'and ops[7][2]==lp and ops[7][-1]==lr
 assert ops[8]==['literal_sparse_columns',D['all106_first_linear_syzygies280_sparse_columns']]
 assert ops[9]==['literal_sparse_rows',D['cubic_NF_dual3_by600_sparse']]
 for j,op in[(10,['matrix_product',5,6]),(11,['matrix_product',2,10]),(12,['matrix_product',7,11]),(13,['matrix_product',2,4]),(14,['matrix_product',3,10])]:assert ops[j]==op
 assert ops[15][0]=='multiply_degree2_relation_columns_by_degree1_coefficients'and ops[15][1:3]==[13,12]
 assert 'ordinary monomial' in ops[15][4]
 assert ops[16]==['matrix_product',9,['matrix_difference',14,15]]
 assert ops[17]==['matrix_product',9,['multiply_degree2_relation_columns_by_degree1_coefficients',13,8]]
 assert ops[18]==['horizontal_blocks',16,17]and D['exported_complete_cubic_operator']==18
 # Fixed-rank data C,D,H,Delta are already actual. Frechet variations here
 # hold them fixed; future primitive variations must also preserve all
 # previously evaluated lower-cache/NF constraints and be coherent.
 sens=D['highest_source_sensitivity_operator'];assert 'deltaM3' in sens['delta_M3']and 'deltaM2*R' in sens['delta_M2_constant_block']and 'deltaM2*C*D' in sens['delta_M2_constant_block']
 joint=load('AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63');jrhs=load('AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63','phase2-known-third-RHS-and-proper-response.json')
 assert joint['all_passed']and not joint['full_general_source105_gauge30_cubic_jet_constructed']
 assert jrhs['monomials5']==[[7,7,7],[0,7,7],[1,7,7],[4,7,7],[5,7,7]]
 gauge=load('AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63','phase1-common105-and-physical-source-base25.json');gmon=[tuple(m)for m in gauge['base25_monomials']]
 assert all(tuple(m)not in gmon for m in jrhs['monomials5'])
 assert D['primitive_highest_chain_gate_open']and not D['actual_M2_M3_or_new_cubic_values_computed']and not D['new_source_gauge_jet']
 proof=dict(whole814_M1_kernel_unit_triangular_witness=dict(free_X_identity646=True,mixed_pivot_d3_identity15=True,free_d3_identity153=True),
  all811_constant_quadratic_kernel_independent_and_complete=True,constant_split=[658,153],all15_mixed_paths_retained=True,
  all106_linear_syzygies_independent_and_complete=True,H174_L_equals_identity_on_imL_verified_all280=True,H174_NOT_right_inverse_on_F180=True,
  actual_Delta3_source_bound=True,all19_nodes_types_and_literal_blocks_verified=True,
  completeness_proof='Expand a=a0+a1+a2. M1a0=0 gives a0=Cc. Projecting M2a0+M1a1=0 by Pi6 gives c=Dd. Since M2CD lies in im L, a1=-R H(M2CD)d + R K e + C v1. Conversely these satisfy both low-order equations. M1a2 lies in the35-linear quadratic-multiple closure; M2C v1 lies in the known full quadratic-multiple closure. Delta3 therefore yields exactly the811+106 additional cubic source images, with no missing coefficient-linear paths.',
  exact_operator='Delta3 [M3 C D - convolution(M2 R,H174 M2 C D) | convolution(M2 R,K106)]',
  exact_formula='HF3(Tnil)=3-rank(actual3x917)',
  actual_total_length=16,actual_HF_through2=[5,5,3],rank_completion={0:[5,5,3,3],1:[5,5,3,2,1],2:[5,5,3,1,1,1]},rank3_impossible_by_Nakayama_and_remaining_length3=True,
  fixed_lower_data_Frechet_signs_checked=True,primitive_variations_must_preserve_known_lower_NF_and_cache=True,
  primitive_highest_sensitivities_and_all917_values_OPEN=True,selected_joint135_not_full_source_or_M3=True,
  two_different_B25_cover100_namespaces_distinguished_by_tuples=True,
  no_new_source_jet_M2_M3_actual_cubic_value_or_tail_choice=True)
 p=OUT/'independent-complete917-source-and-homogeneous-theorem.json';save(p.name,proof)
 for p2 in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p2.relative_to(ROOT))]=sha(p2)
 result=dict(id=ID,all_passed=True,frozen=True,actual_complete_cubic_operator_shape=[3,917],all814_811_106_source_bases_and_H174_verified=True,all15_mixed_and153_d3_paths_retained=True,actual_cubic_values_highest_sensitivities_and_tail_OPEN=True,no_new_source_jet_or_M2_M3_evaluation=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20
 save('report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()

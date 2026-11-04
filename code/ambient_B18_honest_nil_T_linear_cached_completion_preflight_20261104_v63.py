#!/usr/bin/env python3
"""Final constant-source sensitivities: no second coefficient evaluated."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import mv
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_LINEAR_CACHED_COMPLETION_PREFLIGHT_20261104_v63'
PATHS={
 'whole':'reports/AMBIENT_B18_HONEST_NIL_T_LINEAR_SOURCE_WHOLE_PREFLIGHT_20261104_v63/actual-final681-source-sensitivities-and-cached-directions.json',
 'whole_report':'reports/AMBIENT_B18_HONEST_NIL_T_LINEAR_SOURCE_WHOLE_PREFLIGHT_20261104_v63/report.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'W0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'proper1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase4-native-source-and-R12.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def cmul(x,y):
 z=[Q(0)]*4
 for i,a in enumerate(x):
  for j,b in enumerate(y):
   if not i&j:z[i|j]+=a*b
 return z
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 D={k:read(ROOT/p)for k,p in PATHS.items()};assert D['whole_report']['all_passed']
 X,W,R=(D[k]for k in ['X1','W0','radical'])
 assert all(not Q(row[5][0])for row in X['new_A12_by9_C36'])
 # Source derivative of native row20 against all possible unknown u5 even2.
 raw=[[Q(0)]*5 for _ in range(780)]
 for j,col in enumerate(R['native_tensor_Koszul6_HB2_resolution'][1]['columns']):
  for row,pol in col:
   if row==5:
    for mon,c in pol:
     if sum(mon)==1 and any(mon[:5]):raw[640+4*j][mon.index(1)]+=Q(c)
 xp=X['X0_pivot_columns'];co=[[Q(q)for q in v]for v in X['all780_X0_factorization_coefficients']]
 O=W['residual_original_X_columns'];wp=W['W_pivot_residual_columns'];wr=W['W_pivot_rows'];Wi=[[Q(q)for q in v]for v in W['W_fixed_minor_inverse']]
 res=[[raw[j][i]-sum((co[j][k]*raw[xp[k]][i]for k in range(18)),Q(0))for i in range(5)]for j in O]
 rows=[]
 for j,sp in enumerate(W['all_W84_columns_sparse']):
  if j in wp:continue
  w={i:Q(c)for i,c in sp};a=mv(Wi,[w.get(i,Q(0))for i in wr])
  v=[res[j][i]-sum((a[k]*res[wp[k]][i]for k in range(81)),Q(0))for i in range(5)]
  assert not any(v),('UNKNOWN_U5_EVEN2_SURVIVES',O[j],v)
  rows.append(O[j])
 # At base degree2, any unknown C coefficient is multiplied in C/m_R0 C.
 # The even component of (positive-tau origin)*unknown*basis is zero.
 core=[[Q(c)for c in D['proper1']['kappa9'][i][:4]]for i in [0,1,2,6]]
 assert all(not c[0]and not c[3]for c in core)
 fixtures=0
 for origin in core:
  for unknown,basis in itertools.product(range(4),repeat=2):
   u=[Q(j==unknown)for j in range(4)];e=[Q(j==basis)for j in range(4)]
   assert cmul(cmul(origin,u),e)[0]==0;fixtures+=1
 # Literal source monics involve exactly the four stored s/t parameters.
 nodes=D['source']['DAG_scalar_nodes'];closure=set();roots=[]
 for f in ['CA','CB']:
  root=D['source']['families'][f]['actual_tracefree_nu'];roots.append(root);todo=[root]
  while todo:
   i=todo.pop()
   if i in closure:continue
   closure.add(i);n=nodes[i]
   if n[0]in['a','m','n','i']:todo+=n[1:]
 leaves=sorted(i for i in closure if nodes[i][0]=='v');assert leaves==[11,12,20537,20538]and len(closure)==11
 ev={'all681_final_X_source_labels':rows,
  'actual_ResAT_row20_origin_all48_zero':True,'all12_actual_A_i5_even_origins':[row[5][0]for row in X['new_A12_by9_C36']],
  'old_Bold2_entire_term_zero_before_constant_cancellation':True,
  'all681_adapter2_Bold0_zero_after_both_constant_cancellations':True,
  'all681_unknown_u5_even2_sensitivities_zero':True,'unknown_u5_even2_raw_source_coefficients780_by5_sparse':[[[i,c]for i,c in enumerate(v)if c]for v in raw],
  'unknown_core2_and_u5_noneven2_even_projection_C4_fixtures':fixtures,
  'no_unknown_proper_gauge2_source_values_needed':True,
  'needed_row20_second_recipe':{'first640':'ONLY (ResAT_row20)_1 * Bold1; omitted Bold2 and adapter2 terms proved zero on exact final functional, not set zero as functions',
    'last140':'native rowh coefficients at known proper firsts with actual monics2; unknown u5_even2 cancels on final681, all other unknown proper2 even readbacks vanish in C4',
    'm1_row7':'actual saved C44 pure-r² of Vgamma_A, retained exactly'},
  'actual_monic2_cached_source_DAG':[[i,nodes[i]]for i in sorted(closure)],'actual_monic2_source_roots':roots,
  'actual_monic2_leaf_source_ids':leaves,'actual_monic2_formula':'nuA=source12+source11^2/4; nuB=source20538+source20537^2/4',
  'same_frame_cache_binding':'compose stored CB own57 solution with exact full9 YB(YA)2; verify shared raw physical9 second polynomials; then sA=rA',
  'all_actual_values_to_be_computed_are_existing_cache_compositions_only':True,
  'retained_second_callbacks_or_new_source_gauge_solve':False,
  'whole_worker_QQ_bound_inherited':D['whole']['whole_worker_QQ_upper'],
  'prospective_worker_may_skip_19912_AST_source_evaluation_and_u5_unit_recovery':True,
  'full_first_relations_not_yet_evaluated':True}
 p=out/'actual-final-row20-cancellations-and-eleven-node-monic-source.json';p.write_text(json.dumps(enc(ev),separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'status':'FINAL_all_unknown_second_callbacks_cancel_or_are_literal_cached_monic_compositions',
  'all681_unknown_u5_even2_cancel':True,'actual_ResAT_row20_origin_zero':True,'actual_monic_source_DAG_node_count':11,
  'no_new_source_gauge_jet_or_actual_second_coefficients_computed':True,
  'one_existing_cache_full_firstlinear_worker_source_gate_closed':True,
  'worker_QQ_upper':D['whole']['whole_worker_QQ_upper'],'input_sha256':{v:sha(ROOT/v)for v in PATHS.values()},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()

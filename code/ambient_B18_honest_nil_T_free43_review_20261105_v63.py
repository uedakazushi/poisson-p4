#!/usr/bin/env python3
"""Independent full661 free43 row from exact source functions, no solve."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
from ambient_B18_cached_source3_B11_review_helpers_20261104_v63 import ordinary3,coefficient
from ambient_B18_Told_Tnil_zero3_source_review_20261105_v63 import add,sc,sp,mv,mp,eye
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_HONEST_NIL_T_QUADRATIC_FREE43_CACHED_20261105_v63';STEM=A+'_REVIEW'
P={
 'primary':f'reports/{A}/report.json','values':f'reports/{A}/actual-one-by661-free43-quadratic-values-and-source-sensitivities.json',
 'carry':'reports/AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63/actual81-pivot-first-carries-and-six-H1-source-adjoints.json',
 'carry_report':'reports/AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63/report.json',
 'operator':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_661_OPERATOR_20261104_v63/actual-complete-six-by661-quadratic-source-operator.json',
 'quad':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json',
 'X':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'W':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'linear':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'copy1':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json',
 'copy2':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase2-existing-cache-B11-copy-values.json',
 'fourth':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
 'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'atom1':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
 'g2':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/actual-gdiff-second-and-complete619-pure-cubic-values.json',
 'g2_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json',
 'zero3_peer':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_SOURCE_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def enc(p):return [[i,str(c)] for i,c in sorted(p.items())]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'INDEPENDENT_FREE43_SOURCE_AND661_CONSTANT_FOLDS','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()}
 assert all(D[k]['all_passed'] for k in ['primary','carry_report','g2_peer','zero3_peer'])
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 assert D['quad']['exact_quotient_duals_rows_sparse'][3]==[[71,'1']]
 assert D['carry']['all6_actual_effective_H1_callbacks_row84_var8_sparse'][3]==[]
 X,W=D['X'],D['W'];xp=X['X0_pivot_columns'];factors=list(map(sp,X['all780_X0_factorization_coefficients']))
 for j in xp:
  column={int(i):sp(p) for i,p in X['all780_source_columns44_sparse_dual9'][j]};assert not column.get(43,{})
 O=W['residual_original_X_columns'];wp=W['W_pivot_residual_columns'];wr=W['W_pivot_rows'];Wi=list(map(sp,W['W_fixed_minor_inverse']));wc=list(map(sp,W['all_W84_columns_sparse']))
 wm=[{j:wc[c].get(r,Q(0)) for j,c in enumerate(wp) if wc[c].get(r,Q(0))} for r in wr]
 assert mp(wm,Wi)==mp(Wi,wm)==eye(81)
 def residual(j):return add({j:Q(1)},{xp[k]:-c for k,c in factors[j].items()})
 remaining=[i for i in range(762) if i not in wp];assert [O[i] for i in remaining]==D['linear']['actual_final_source_labels681']
 final=[]
 for k in remaining:
  coef=mv(Wi,{i:wc[k][r] for i,r in enumerate(wr) if r in wc[k]})
  final.append(add(residual(O[k]),*(sc(residual(O[wp[j]]),-a) for j,a in coef.items())))
 weights=[]
 for column in D['operator']['all661_source_columns_in_original849']:
  weights.append(add(*(sc(final[j],Q(c)) for j,c in column if j<681)))
 assert list(map(enc,weights))==D['values']['all661_constant_original780_source_combinations']
 rs=D['radical']['new_Cscalar_and_R0component_nodes'];pcols=D['radical']['actual_P8_by640_columns'];ids=sorted({col[7] for col in pcols if col[7]})
 assert len(ids)==94 and all(rs[i][0]=='R0_original_honest_C0_scalar' for i in ids)
 roots=[rs[i][1] for i in ids];readbacks,sourceproof=ordinary3(D,roots)
 g2={i:Q(c) for i,c in D['g2']['actual_eight_regular_g_second_values']};rawvalues={};rawsens={};rootrecords=[]
 for i,root in zip(ids,roots):
  known,sens=coefficient(readbacks[root],2);unresolved={}
  for (leaf,tag),c in sens.items():
   assert tuple(tag)==(7,7)
   if leaf in g2:known+=c*g2[leaf]
   else:unresolved[leaf]=unresolved.get(leaf,Q(0))+c
  unresolved={i:c for i,c in unresolved.items() if c};rawvalues[i]=known;rawsens[i]=unresolved
  rootrecords.append([i,root,str(known),enc(unresolved)])
 # Match source values after independent replay, allowing the producer to
 # retain its provenance label for the four literal g-difference roots.
 assert all([a,b,c,d]==[i,root,str(Q(v)),[[j,str(Q(c))] for j,c in sensitivity]] for (a,b,c,d),(i,root,v,sensitivity,label) in zip(rootrecords,D['values']['all94_raw_Pg_roots_and_exact_second_values_or_sensitivities']))
 raw=[(rawvalues[col[7]],rawsens[col[7]]) if col[7] else (Q(0),{}) for col in pcols]+[(Q(0),{})]*140
 actual=[];residual=[]
 for weights0 in weights:
  value=sum(c*raw[j][0] for j,c in weights0.items());sens=add(*(sc(raw[j][1],c) for j,c in weights0.items()))
  actual.append(value);residual.append(sens)
 assert all(not p for p in residual)
 assert list(map(str,actual))==D['values']['all661_actual_quadratic_dual_values']
 nonzero=[i for i,c in enumerate(actual) if c];assert nonzero==D['values']['nonzero_actual_value_indices'] and len(nonzero)==23
 assert actual[570]==Q(356092408505725,80649213849504)
 ep=o/'independent94-source-roots-all661-free43-values.json'
 dump(o,ep.name,{'inputs':inputs,'actual_minimal_target1_free43_from_independent_ZERO3_peer':True,
  'literal_dual3_180':[[71,'1']],'all18_Xp_free43_fullfirsts_zero':True,'actual_H1_carry_full_source_dual_zero':True,
  'all661_original780_source_weights_independently_reconstructed':list(map(enc,weights)),
  'all94_roots_independent_source2_readbacks':rootrecords,'source_equation_and_copy_proof':sourceproof,
  'all661_true_honest_second_sensitivities_zero':True,'all661_actual_free43_quadratic_values':list(map(str,actual)),
  'nonzero_indices':nonzero,'nonzero_count':23,'example570':str(actual[570]),
  'actual_known6_quotient_rank_at_least1':True,'no_exact_Tnil_HF2_total_length_or_nilcube_claim':True,
  'new_source_or_jet_or_unknown_solve':False})
 note=ROOT/'research/notes'/f'{STEM}.md';report={'id':STEM,'all_passed':True,
  'all94_source_roots_and661_actual_free43_values_independently_accepted':True,'all661_honest_second_sensitivities_zero':True,
  'actual_one_row_nonzero_count':23,'actual_HF2_or_nilcube_not_determined':True,
  'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()

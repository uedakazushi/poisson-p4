#!/usr/bin/env python3
"""One actual quadratic row of complete6x661, using immutable root2 caches.
Unknown honest-second leaves stay independent until the full source fold.
"""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import enc,sparse,nqq
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_QUADRATIC_FREE43_CACHED_20261105_v63'
PATHS={
 'carry':'reports/AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63/actual81-pivot-first-carries-and-six-H1-source-adjoints.json',
 'carry_report':'reports/AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63/report.json',
 'operator':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_661_OPERATOR_20261104_v63/actual-complete-six-by661-quadratic-source-operator.json',
 'quad':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'W':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'first':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'cache':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD/phase3-oldB-common-source2-values-and-symbols.json',
 'cache_report':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD/report.json',
 'cache_peer':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD_REVIEW/report.json',
 'g2':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/actual-gdiff-second-and-complete619-pure-cubic-values.json',
 'g2_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json'}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(d,k,v):
 if v:
  d[k]=d.get(k,Q(0))+v
  if not d[k]:del d[k]
def scaledadd(d,v,c):
 for k,a in v.items():put(d,k,c*a)
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 (out/'preflight.json').write_text(json.dumps({'whole_saved_plus_live_QQ_upper':60000,'seconds':60,'address_bytes':2<<30,'output_bytes':32<<20,'scope':'ONLY old-P_g rA² cached roots plus literal four g-differences; opaque honest seconds retain true IDs','new_source_gauge_or_higher_jet':False},indent=2)+'\n')
 D={k:read(ROOT/p)for k,p in PATHS.items()}
 assert all(D[k]['all_passed']for k in ['carry_report','cache_report','cache_peer','g2_peer'])
 (out/'phase.json').write_text(json.dumps({'phase':'all_immutable_inputs_bound_before_source_fold','sha256':{p:sha(ROOT/p)for p in PATHS.values()}},indent=2)+'\n')
 assert D['quad']['exact_quotient_duals_rows_sparse'][3]==[[71,'1']]
 assert D['carry']['all6_actual_effective_H1_callbacks_row84_var8_sparse'][3]==[]
 X,W=D['X1'],D['W'];xp=X['X0_pivot_columns'];co=[[Q(c)for c in row]for row in X['all780_X0_factorization_coefficients']]
 for j in xp:
  assert all(not Q(c)for row,sp in X['all780_source_columns44_sparse_dual9'][j]if row==43 for k,c in sp)
 O=W['residual_original_X_columns'];wp=W['W_pivot_residual_columns'];wr=W['W_pivot_rows'];Wi=[[Q(c)for c in row]for row in W['W_fixed_minor_inverse']]
 final=[i for i in range(762)if i not in wp];assert [O[i]for i in final]==D['first']['actual_final_source_labels681']
 # Full constant source maps, kept factored/sparse rather than780x661 dense.
 first={j:{j:Q(1),**{p:-c for p,c in zip(xp,co[j])if c}}for j in O}
 fweights=[]
 for k in final:
  w={i:Q(c)for i,c in W['all_W84_columns_sparse'][k]};nz=[(i,w[r])for i,r in enumerate(wr)if r in w]
  a=[sum((row[i]*c for i,c in nz),Q(0))for row in Wi]
  v=dict(first[O[k]])
  for p,c in enumerate(a):
   if c:scaledadd(v,first[O[wp[p]]],-c)
  fweights.append(v)
 sourceweights=[]
 for col in D['operator']['all661_source_columns_in_original849']:
  v={}
  for j,c in col:
   if j<681:scaledadd(v,fweights[j],Q(c))
  sourceweights.append(v)
 R=D['radical'];rs=R['new_Cscalar_and_R0component_nodes'];pcols=R['actual_P8_by640_columns'];Cnodes=D['C0']['DAG_scalar_nodes']
 ids=sorted({p[7]for p in pcols if p[7]})
 assert all(rs[i][0]=='R0_original_honest_C0_scalar'for i in ids)
 cached={i:({j:Q(c)for j,c in v},{j:Q(c)for j,c in g})for i,v,g in D['cache']['C0_composite_values_and_opaque_second_gradients']}
 mon20=[tuple(m)for m in D['cache']['base20_monomials_original8']];r2=mon20.index((7,7));g2={i:Q(c)for i,c in D['g2']['actual_eight_regular_g_second_values']}
 vals={};missing=[];rootrecords=[]
 for i in ids:
  root=rs[i][1]
  if root in cached:
   v,g=cached[root];z=v.get(r2,Q(0));sens={}
   for leaf,c in g.items():
    if leaf in g2:z+=c*g2[leaf]
    else:put(sens,leaf,c)
   vals[i]=(z,sens);rootrecords.append([i,root,z,sparse_from_dict(sens),'immutable_root_cache_plus_exact8g2'])
  else:
   # All four missing roots are exactly oldCB_g - newCA_g, no unspecified callback.
   n=Cnodes[root]
   if n[0]=='a' and n[1]in g2 and Cnodes[n[2]][0]=='n' and Cnodes[n[2]][1]in g2:
    z=g2[n[1]]-g2[Cnodes[n[2]][1]];vals[i]=(z,{});rootrecords.append([i,root,z,[],'literal_difference_of_actual_g2'])
   else:missing.append([i,root,n])
 assert not missing,('UNAVAILABLE_SOURCE_ROOT',missing)
 assert len(rootrecords)==94 and sum(r[4]=='literal_difference_of_actual_g2'for r in rootrecords)==4
 raw=[vals[p[7]]if p[7]else(Q(0),{})for p in pcols]+[(Q(0),{})]*140
 known=[];sens=[]
 for source in sourceweights:
  z=Q(0);a={}
  for j,c in source.items():
   v,g=raw[j];z+=c*v;scaledadd(a,g,c)
  known.append(z);sens.append(a)
 allzero=not any(sens)
 data={'actual_complete_dual_index0based':3,'dual_label':'newfree43*rA² (target1 monomial35)',
  'all18_Xp_free43_origins_and_full_firsts_literal_zero':True,'actual_first_pivot_carry_dual_literal_zero':True,
  'all94_raw_Pg_roots_and_exact_second_values_or_sensitivities':rootrecords,
  'source_coverage':{'cached_originalC0_roots':90,'literal_g_difference_roots':4,'new_source_callbacks':0},
  'all661_constant_original780_source_combinations':[sparse_from_dict(v)for v in sourceweights],
  'all661_folded_known_free43_second_r_values':known,'all661_remaining_true_honest_second_sensitivities':[sparse_from_dict(v)for v in sens],
  'all_remaining_honest_second_sensitivities_cancel':allzero,
  'all661_actual_quadratic_dual_values':known if allzero else None,
  'nonzero_actual_value_indices':[i for i,c in enumerate(known)if c]if allzero else None,
  'source_proof':'The exact dual is top free43*r². Xp free43 has order>=2, so all degree1-unit correction terms vanish on a constant-source zero-M1 column. Actual C1 dual is zero, including mixed d3 columns. Thus its quadratic value is precisely the constant-source fold of raw P_g second-r coefficients. All raw regular quotients retain exact leaf IDs; only the eight source-proved g2 values are substituted.',
  'no_other_five_dual_or_actual_HF2_rank_total_length_or_nilcube_claim':True,'new_source_gauge_jet_or_X3_m2_evaluation':False}
 assert nqq(data)<150000
 p=out/'actual-one-by661-free43-quadratic-values-and-source-sensitivities.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;p.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'PROVED_actual_one_by661_row'if allzero else'SOURCE_GATE_remaining_honest_seconds_retained','all94_raw_source_roots_bound':True,'all661_unknown_honest_second_sensitivities_cancel':allzero,'actual_one_by661_nonzero_count':sum(bool(c)for c in known)if allzero else None,'remaining_sensitive_columns':sum(bool(g)for g in sens),'new_source_gauge_jet_or_X3_m2_evaluation':False,'visible_QQ_support':nqq(data),'bytes':p.stat().st_size,'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
def sparse_from_dict(d):return [[i,c]for i,c in sorted(d.items())if c]
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Four missing second-pair gauge-input preflight; no new gauge coefficient."""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_native128_first_20261104_v63 import PATHS as OLD_PATHS
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_GAUGE120_SOURCE_PREFLIGHT_20261105_v63'
PATHS={k:OLD_PATHS[k]for k in ['source','interface','proper','second5','preparation']}
PATHS.update({
 'sensitivity':'reports/AMBIENT_B18_HONEST_NIL_T_HONEST_SECOND_SENSITIVITY_20261105_v63/phase3-full661-honest-second-sensitivities.json',
 'source2':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'adapter2':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'source20':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63/phase1-actual105-source2-and-copy.json',
 'gauge1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase3-gauge30-firsts.json',
 'gauge20':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63/phase2-canonical270-zz-gauge.json',
 'proper20':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63/phase3-actual-proper-G8-base20.json',
 'J30':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase2-fixed30-Jacobian.json',
 'higher661':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_HIGHER_SOURCE_PREFLIGHT_20261104_v63/actual661-first-corrections-m2-replacement-and-native-highest-sensitivities.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json'})
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 sm=D['source20']['base20_monomials_original8'];mon20=[tuple(m)for m in sm];missing=[(0,6),(4,6),(5,6),(6,7)];mon25=mon20+[(6,)]+missing
 assert len(set(mon25))==25
 for mon in mon25:
  assert all(mon[:j]+mon[j+1:]in mon25 for j in range(len(mon)))
 ad=D['sensitivity']['adjoint_deltaX2_6_by36_by44_sparse'];mons=[tuple(m)for m in D['sensitivity']['monomial_order36']]
 need={mons[q]for row in ad for q,v in enumerate(row)if v};assert need-set(mon20)==set(missing)
 assert D['higher661']['all661_full44_by8_X1u_equals_Xp0_a1_verified']
 m0=[[Q(c)for c in row]for row in D['m0']['m0_24_by44']];first=[]
 for col in D['X1']['all780_source_columns44_sparse_dual9']:first.append({i:dict((k,Q(c))for k,c in sp)for i,sp in col})
 contracts=[]
 for row in [15,23]:
  contracts.append([[-sum((m0[row][i]*v.get(8,Q(0))for i,v in first[p].items()),Q(0))for p in D['X1']['X0_pivot_columns']]])
 # Exact fixed source-liveness count, before any finite-algebra evaluation.
 S=D['source'];nodes=S['DAG_scalar_nodes'];bind={i for i,v in D['source20']['all105_actual_source2_values']}
 roots=[i for f in S['families'].values()for pair in f['complete_normalized_physical_w_coefficients'][0]for i in pair]+[f['actual_tracefree_nu']for f in S['families'].values()]
 needed=set();todo=roots[:];bad=[]
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i)
  if i in bind:continue
  n=nodes[i];o=n[0]
  if o in ['a','m']:todo.extend(n[1:3])
  elif o in ['n','i']:todo.append(n[1])
  elif o!='c':bad.append([i,n])
 assert not bad
 uses=collections.Counter()
 for i in needed:
  if i not in bind and nodes[i][0]in ['a','m','n','i']:uses.update(nodes[i][1:])
 uses.update(roots);live=peak=0
 for i in sorted(needed):
  live+=25;peak=max(peak,live)
  if i not in bind and nodes[i][0]in ['a','m','n','i']:
   for j in nodes[i][1:]:
    uses[j]-=1
    if not uses[j]:live-=25
 budget={'physical_source_DAG_dense_peak':peak,'source105_base25':105*25,'source_parameter_pairs_and_nu':(2*150*2+2)*25,'gauge30_base25':30*25,'fixed30_inverse_and_rhs':900+30*4,'proper14_C100_and_kappa':23*100,'finite_GL_and_source_work_arrays':40000,'save_and_metadata_upper':15000}
 assert sum(budget.values())<150000,('WHOLE_BUDGET_GATE',budget)
 data={'base25_original8_monomials':mon25,'all15_actual_raw_X2_adjoint_pair_support':sorted(need),'new_second_pairs_only':missing,'new_gauge_second_scalar_count':120,'new_physical_source_second_count':0,
  'known_first_z6_seed':'all30 even/odd first9 entries in gauge1; embed this first coefficient, not0',
  'known_second_seed':'all existing base20 gauge values, unchanged; project exactly back after extension',
  'new_correction_recipe':'For four new pairs evaluate30 scalar constraints with cached physical source2, then correction=-J30_inverse*RHS. The fixed unit Jacobian is the same canonical auxiliary-etaA gauge, not physical-tau differentiation.',
  'all105_physical_second_source_input':'two frozen own57 degree2 arrays plus exact full9 CB_F9_in_CA_F9 adapter; sA=rA. No implicit Newton source solve.',
  'pure_r_last_m1_contraction_replacements_rows15_23':contracts,
  'pure_r_last_m1_elimination_proof':'For every original Xp, restrict d1*m*Xp=0 to z=0. Each sector r^n is a nonzerodivisor, so (mXp)_last=0. First-r coefficient gives m1_last,r Xp0=-m0_last Xp1_r. Full661 X1u=Xp0*a1 replaces every such contraction by saved18 first callbacks. In particular no third-sector Vgamma_r3 is evaluated.',
  'source_DAG_ancestor_count':len(needed),'no_unbound_source_leaves':True,'aggregate_QQ_upper':sum(budget.values()),'budget':budget,
  'proposed_future_gauge120_values_NOT_computed':True,'remaining_work':'Only after a separate parent GO, build base25 pure-source/canonical-gauge coefficients and lower m1/carry from these actual functions; full661 values still open.'}
 def enc(x):
  if isinstance(x,Q):return str(x)
  if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
  if isinstance(x,(list,tuple)):return [enc(v)for v in x]
  return x
 p=out/'four-pair-base25-and-last-m1-contraction-preflight.json';p.write_text(json.dumps(enc(data),separators=(',',':'))+'\n')
 r={'id':STEM,'all_passed':True,'status':'SOURCE_AND_BUDGET_GATE_ONLY_gauge120_not_run','new_required_second_pairs':missing,'proposed_new_gauge_coefficients':120,'new_physical_source_coefficients':0,'aggregate_QQ_upper':sum(budget.values()),'no_new_source_gauge_jet_value':True,'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':main()

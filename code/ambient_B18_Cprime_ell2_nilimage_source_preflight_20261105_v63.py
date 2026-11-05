#!/usr/bin/env python3
"""Source/budget gate for Cprime/ell² from existing Artin33 physical input.
No canonical-gauge values, nil images, rank or connecting coefficients.
"""
import collections,hashlib,json,time
from pathlib import Path
from fractions import Fraction as Q
from ambient_B18_honest_nil_T_gauge120_source_preflight_20261105_v63 import PATHS as GP
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CPRIME_ELL2_NILIMAGE_SOURCE_PREFLIGHT_20261105_v63'
PATHS={k:GP[k]for k in ['source','interface','proper','second5','preparation','gauge1','J30']}
PATHS.update({
 'C2report':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/report.json',
 'C2peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'source33':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'table33':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/exact-existing-raw33-multiplication-table.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'nil2':'reports/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63/actual-E4-nil-two-generators-and-relations.json',
 'nil2peer':'reports/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63_REVIEW/report.json',
 'Cprime':'reports/AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63/actual-eight-nil-images-and-reduced-Cbar-cokernel.json',
 'rank2':'reports/AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63/report.json',
 'nilfinite':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json',
 'firstpreflight':'reports/AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63/actual-native128-first-input-and-mixed-index-profile.json'})
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def main():
 t=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=True);assert not(out/'report.json').exists();pins={p:sha(ROOT/p)for p in PATHS.values()};(out/'phase.json').write_text(json.dumps({'phase':'source-only pinned input inspection; no value worker','input_sha256':pins},indent=2)+'\n');D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 for k in ['C2report','C2peer','nil2peer','rank2','nilfinite']:assert D[k]['all_passed']
 S,T,A=D['source'],D['source33'],D['table33'];assert T['all105_equations_exactly_zero']and T['finite_actual_function_uniqueness'] and A['m5_zero'];assert D['C2']['actual_C2_dimension']==14
 ids=set(T['joint105_unknown_order']);assert len(ids)==105;nodes=S['DAG_scalar_nodes'];roots=[i for f in S['families'].values()for p in f['complete_normalized_physical_w_coefficients'][0]for i in p]+[f['actual_tracefree_nu']for f in S['families'].values()];seen=set();todo=list(roots);bad=[]
 while todo:
  i=todo.pop()
  if i in seen:continue
  seen.add(i)
  if i in ids:continue
  n=nodes[i];o=n[0]
  if o in ['a','m']:todo.extend(n[1:3])
  elif o in ['n','i']:todo.append(n[1])
  elif o!='c':bad.append([i,n])
 assert not bad
 uses=collections.Counter(roots)
 for i in seen:
  if i not in ids and nodes[i][0]in ['a','m','n','i']:uses.update(nodes[i][1:])
 live=peak=0
 for i in sorted(seen):
  live+=33;peak=max(peak,live)
  if i not in ids and nodes[i][0]in ['a','m','n','i']:
   for j in nodes[i][1:]:
    uses[j]-=1
    if not uses[j]:live-=33
 assert peak==13299
 # Same fifteen group positions, same frame5+fiber10 constraints and same
 # scalar even/odd unit Jacobian; no new flag or frame equivalence premise.
 assert len(D['interface']['actual_group15_positions'])==15
 assert len(D['proper']['CA_fiber10_fixed_coefficient_indices'])==10
 assert len(D['proper']['CA_free_fiber9_coefficient_indices'])==9
 assert Q(D['J30']['det_J30'])==Q(4,9)
 assert len(D['J30']['J30_inverse'])==30
 basis=A['basis33'];assert len(basis)==33
 # The parent edge-functional refinement will require the FORMAL plane
 # quotient z7=0,r²=ell²=0. Check that this quotient factors through33.
 plane=[(0,0),(1,0),(0,1),(1,1)];pidx={m:i for i,m in enumerate(plane)};maps=[]
 for z,r,e in basis:maps.append({pidx[(r,e)]:Q(1)}if z==0 and(r,e)in pidx else{})
 def mul(x,y):
  v={}
  for i,c in x.items():
   for j,d in y.items():
    m=tuple(a+b for a,b in zip(plane[i],plane[j]));k=pidx.get(m)
    if k is not None:v[k]=v.get(k,Q(0))+c*d
  return {k:c for k,c in v.items()if c}
 def project(sp):
  v={}
  for k,c in sp:
   for i,d in maps[k].items():v[i]=v.get(i,Q(0))+Q(c)*d
  return {k:c for k,c in v.items()if c}
 assert all(project(A['table_columns_for_each_left_basis'][i][j])==mul(maps[i],maps[j])for i in range(33)for j in range(33))
 budget={'physical_source_DAG_dense_peak':peak,'existing_source105_base33':3465,'physical150_pairs_nus':19866,'gauge30_base33':990,'fixedJ_inverse_RHS':2880,'actual14_and_kappa9_C132':3036,'finite_GL_source_stream_work':52800,'four_gauge_residual_correction_stages':7920,'nil8_images_and_rank_action_receipts':6000,'metadata_and_product_table':8000};whole=sum(budget.values());assert whole<150000
 data={'actual_base':'R0/(IA*IB,ell²), actual CA-F9 z7,rA,sA=rA+ell/5;33basis is source-bound, not the earlier ell0H20 or jet base25',
  'existing_source105_fixed_in_artin33':True,'new_physical_source_or_jet_values_needed':False,'all_physical150_root_ancestors_bound':True,'source_ancestor_count':len(seen),'dense_source_stream_peak_QQ':peak,
  'canonical_gauge_recipe':'Reuse actual group15 positions; group G=I+even15+odd15*tauA, same physical-source CA frame and fixed10 fiber coefficients. Read the cached gauge30 origin, solve30 equations overArtin33 by four fixed J30 inverse corrections at orders1..4, checking residual30=0 exactly. Do not rerun original105 or invent a gauge seed.',
  'new_derived_gauge_coefficient_upper':990,'gauge_origin_unit_J30':Q(4,9),'gauge_origin_and_existing_ell0_overlap_required':True,
  'trace_algebra':'C132=Artin33[tauA,tauB]/(tauA²-nuA,tauB²-nuB), fixed16 block products via sparse33table; do not materialize a dense132multiplication table',
  'trace_unit_inverse_proof':'nuA,nuB are in the base maximalideal; base m5=0 implies trace maximalideal^11=0, so the exact finite geometric inverse uses powers0..10, only at nonzero constant units',
  'actual_nil_frame':'Read actual proper kappa9 from the SAME canonical group and source14 map. Two generators are d*h and d*c*l, with h=a-d+x-2t,c=(d+x)/2,l=(a+d-x)/3; allfive d*u relations retained.',
  'actual_Cprime2_recipe':'Multiply each actual nilgen by four tracebasis elements. Drop scalar component and project tauA toA6,tauB toB6,tauAB to33 via frozen raw45 common quotient maps; then use actual45→C2 projection. Close all8images under the33 actual coefficient-basis actions (264columns), verify invariance under9actualR0actions. Right-exact quotient is genuine Cprime/ell².',
  'expected_Cprime2_interval_from_rank2_and_tors_mu2':[6,8],
  'if_length6':'Both DVR torsion exponents equal1. Hence Cprime[ell] injects canonically into Cprime/ell; its image is represented by ker(ell:Cprime2→Cprime2)/imell and can be extracted with allR0actions. No division byell inArtin33.',
  'if_length7_or8':'Leave genuine Cprime[ell] and all exponent≥3 possibilities open; no finite-kernel identification or14nilcycle assertion.',
  'formal_connecting_formula_not_evaluated':'For an ell-independent FORMAL lift of each ell0cycle cbar, delta=[(partial_ell V_full)|ell0*cbar] inCprime[ell], not a division inArtin33. If torsionell-ann1 is proved, project this expression to actualCprimebar; target m²=0, so only head/firstcycle coefficients and formalell1+S1 source coefficients matter.',
  'formal_plane4_quotient_maps33':[[[i,c]for i,c in sorted(v.items())]for v in maps],'all1089_plane_quotient_products_verified':True,'plane4_scope':'z7=0,rA²=ell²=0 retains the needed formalell/ell*r coefficient functional; this does not make arbitrary ell division legitimate or close the separate cycle-first sensitivity gate.',
  'head_cycle_first_gate_still_required':'The two free-head cycles have c0=(I-n0m0)e42,e43 and c1=-n0*m1*c0. Keep third-sector last-row regular division sensitivities until their projected cancellation/source binding. Keeper wedge heads use c1=n0*d2_linear*b. No m1 callback is assignedzero.',
  'hard_budget':budget,'whole_saved_live_QQ_upper':whole,'hard_seconds':60,'hard_memory_bytes':2<<30,'hard_total_output_bytes':32<<20,'proposed_new_retained_QQ_upper':18000,'future_derived_coefficient_bit_guard':2048,'no_value_or_gauge_job_run':True}
 p=out/'source33-nilimage-Cprime2-and-formal-edge-gate.json';p.write_text(json.dumps(enc(data),separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_PROFILE_ONLY_no_new_gauge_or_Cprime2_values','all_needed_physical_source_already_exact33':True,'same_canonical30_unit_equations_Jacobian':True,'derived_finite_gauge_value_job_not_run':True,'whole_saved_live_QQ_upper':whole,'actual_Cprime2_length_or_connecting_kernel_decided':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()

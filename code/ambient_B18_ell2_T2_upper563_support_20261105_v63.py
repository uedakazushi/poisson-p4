#!/usr/bin/env python3
"""Coefficient-free original-source ancestor inventory for upper563."""
import argparse,json,resource,signal,time
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER563_SUPPORT_PREFLIGHT_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 source=read('reports/AMBIENT_B18_ELL2_T2_UPPER619_SOURCE_CLEARING_20261105_v63/report.json');assert source['all_passed']
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER619_SOURCE_CLEARING_20261105_v63_REVIEW/report.json')['all_passed']
 prim=read('reports/AMBIENT_B18_ELL2_T2_NATIVE70_VALUES_20261105_v63/actual-native70-lower-prefix-and-witness-if-any.json');c=prim['nonzero_witness'];assert c['case619']==563
 assert read('reports/AMBIENT_B18_ELL2_T2_NATIVE70_WITNESS563_20261105_v63_REVIEW/report.json')['all_passed']
 mini='reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json';piv=field(mini,'original27_X_source_columns')
 mu='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json';sel=field(mu,'actual_selected_Y0_source_columns')
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';ker=field(mug,'all741_zero_origin_source_combinations');stars=[ker[i]['original768'] for i in sel];del ker
 B=[dict(v) for v in field(mug,'old640_source_columns')];ops=dict(field(mug,'old640_typed_closure_nodes'));P=field(H.PATHS['C0'],'honest_D0_relation_columns')
 gp='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json';gamma=field(gp,'Gamma_two');rows=sorted({i for g in gamma for i,x in g})
 def inventory(original):
  typed=set();opaque=set();oldids=[];native=[]
  for j,a in original:
   if j>=640:native.append(j-640);continue
   oldids.append(j)
   for row in rows:
    if row<48:
     k=B[j].get(row)
     if k is None:continue
     n=ops[k]
     if n[0]=='c':continue
     if n[0]=='actual_C0_scalar_function':opaque.add(n[1])
     else:typed.add(k)
    else:opaque.add(P[j][row-48+1])
  return {'original_old_source_ids':oldids,'original_native128_source_ids':native,'typed_roots':sorted(typed),'literal_C0_scalar_roots':sorted(opaque)}
 cubic=inventory(c['original768_constant_combination']);starinv=[inventory(v) for v in stars]
 a1plane=[[slot,x] for slot,x in c['actual_a1_27_times9'] if slot%9 in [7,8]];assert a1plane==[[97,'-75/2'],[98,'-21/5']] and piv[10]==664
 pivotinv=inventory([[664,'1']]);needed={i for d in [cubic,pivotinv]+starinv for i in d['typed_roots']};closure=set()
 def visit(i):
  if i in closure:return
  closure.add(i);n=ops[i]
  kids=[j for j,a,b in n[1]] if n[0]=='sum' else n[3] if n[0]=='product' else []
  for j in kids:visit(j)
 for i in needed:visit(i)
 deps={}
 for i in sorted(closure):
  n=ops[i]
  if n[0]=='actual_C0_scalar_function':v={n[1]}
  elif n[0]=='c':v=set()
  else:v=set().union(*(deps[j] for j in ([j for j,a,b in n[1]] if n[0]=='sum' else n[3])))
  deps[i]=v
 for d in [cubic,pivotinv]+starinv:
  labels=set(d['literal_C0_scalar_roots'])
  for i in d['typed_roots']:labels.update(deps[i])
  d['coefficient_free_scalar_ancestor_labels']=sorted(labels)
 candidate={
 'actual_N_rr':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/phase2-actual-N-B11.json',
 'actual_N_rell':'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json',
 'actual_canonical_kappa33':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/phase3-actual-kappa9-and-normalized-nilframe.json',
 'actual_canonical_gauge33':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/phase2-canonical30-exact-values.json'}
 for p in candidate.values():hs[p]=H.sha(ROOT/p)
 retained=count([B,ops,P,stars,c,gamma,starinv,deps])+6000
 formal_upper=retained+6*len(closure)+12*sum(len(d['coefficient_free_scalar_ancestor_labels']) for d in [cubic,pivotinv]+starinv)+10000
 # Future only: do not construct a918-by99 array or any signed root fold.
 ep=out/'one563-coefficient-free-cubic-second-and-unit-carry-ancestor-inventory.json';H.dump(ep,{'case619':563,'original768_constant_combination':c['original768_constant_combination'],'actual_a1_27_times9':c['actual_a1_27_times9'],'plane_a1_only_pivot664':a1plane,'Gamma_two':gamma,'direct_cubic_Xc_source_inventory':cubic,'individual_pivot_second_inventory':pivotinv,'selected99_original768_constant_combinations':stars,'selected99_second_source_inventories':starinv,'typed_closure_nodes':[[i,ops[i]] for i in sorted(closure)],
 'b1_A2_source_operators_remain_factored':'b1=B99*(I-G3_Jlinear*F25*PiY0)*Ybar1 on full918 slots;A2=P27*(Xc2-Xp1a1-Xstar1b1)',
 'all99_b1_coefficients_and_all_scalar_seconds_thirds_unevaluated':True,'native_cubic_retained_terms':['L0S3_r2ell','L1rS2_rell','L1ellS2_rr','L2rrS1ell','L2rellS1r'],'native_L3_drop_only_by_separate_fullC4_S0_kernel':True,
 'candidate_cache_paths_not_availability_claims':candidate,'prospective_formal_coefficient_ring_basis':['1','r','ell','rr','rell','rr*ell'],'no_physical_or_signed_root_fold_or_actual_upper_value':True,'future_formal_fold_numeric_support_forecast':formal_upper,'future_height_RAW_and_full_source_adjoint_gate_pending':True})
 report={'id':STEM,'all_passed':True,'coefficient_free_source_inventory_only':True,'actual_upper_values_or_source_jets_computed':False,'typed_closure_nodes':len(closure),'direct_cubic_scalar_ancestor_count':len(cubic['coefficient_free_scalar_ancestor_labels']),'selected99_second_scalar_label_union':len(set().union(*(set(v['coefficient_free_scalar_ancestor_labels']) for v in starinv))),'retained_numeric_bound':retained,'prospective_formal_support_forecast':formal_upper,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};assert retained<60000 and formal_upper<150000;H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','typed_closure_nodes','direct_cubic_scalar_ancestor_count','selected99_second_scalar_label_union','retained_numeric_bound','prospective_formal_support_forecast']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

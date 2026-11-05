#!/usr/bin/env python3
"""Source-complete2x619 first-Hilbert closure; no new source coefficients."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=H.ROOT;GATE='AMBIENT_B18_ELL2_T2_FIRST_HILBERT_RANK_INPUT_GATE_20261105_v63';STEM='AMBIENT_B18_ELL2_T2_FIRST_HILBERT_RANK_20261105_v63'
OLD='reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63/actual104-source-rell-values-and-two-by619-old-knownparts.json'
NATIVE='reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json'
def gate(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 o=read(OLD);n=read(NATIVE);old=o['all619_two_old_knownparts_including_frozen_pivot_first_correction'];native=n['all619_two_actual_native_parts'];assert len(old)==len(native)==619 and all(len(v)==2 for v in old+native)
 for v in old+native:H.guard({i:F(c) for i,c in enumerate(v) if F(c)})
 assert o['all_honest_second_rell_sensitivities_exactly_zero'] and n['all619_nativeS0_zero_rechecked'] and n['L1ell_Omega_r_retained'] and n['negative_intrinsic_pivot_correction_already_in_old_knownparts_not_added_again']
 paths=['reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63/report.json','reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/report.json','reports/AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json']
 for p in paths:assert read(p)['all_passed']
 known=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/report.json');assert known['known_G3_top_linear_rank']==25 and known['known_bottom_linear_rank']==23 and known['zero_bottom_kernel_count']==619
 count=2*619*2+2*619+2*2*3+3000
 ep=out/'actual-two-by619-first-Hilbert-complete-source-input.json';H.dump(ep,{'actual_old_two_by619':old,'actual_native_two_by619':native,
  'complete_source_formula':'column=oldknownparts+native(L0S2rell+L1rOmegaEll+L1ellOmegaR); negativeintrinsicpivotcorrectionisalreadyold, appliedonce',
  'source_order':'same619fixedzero-bottom/Y0 original768combinations fromknownfirstprofile',
  'full_first_relation_completeness':'minimal6x1060 hasallconstants0. Degree1relations useonlyconstantcoeffs on1060columns. 418G3 giveknownrank25top,642Xgivebottomrank23; their619kernelprojectsviaverifiedfull918adjointstothese2rows.',
  'exact_HF1_formula':'9*6-(25+23+rank2by619)=6-rank2by619',
  'exact_length_mod_m2_formula':'mu6+HF1=12-rank2by619',
  'actual_T2_total_length_and_K_exponents_not_inferred':True,'native_postpeer_pending_at_gate':True,
  'whole_slot_bound':count,'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'QQ_coefficients':150000,'bits':384},'combined_addition_or_rank_started':False})
 report={'id':GATE,'all_passed':True,'combined_addition_or_rank_started':False,'all_actual_source_values_present':True,'native_postpeer_pending':True,
  'whole_slot_bound':count,'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','whole_slot_bound','seconds']}))
def values(out,gpath,nativepeer):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 gp=ROOT/gpath;g=json.loads(gp.read_text());assert g['all_passed'] and not g['combined_addition_or_rank_started'];ip=gp.parent/'actual-two-by619-first-Hilbert-complete-source-input.json';assert H.sha(ip)==g['evidence_sha256'][ip.name];d=json.loads(ip.read_text())
 hs={gpath:H.sha(gp),str(ip.relative_to(ROOT)):H.sha(ip)}
 if nativepeer:
  p=ROOT/nativepeer;r=json.loads(p.read_text());assert r['all_passed'];hs[nativepeer]=H.sha(p)
 for p,h in g['input_sha256'].items():assert H.sha(ROOT/p)==h;hs[p]=h
 cols=[H.guard({i:F(a)+F(b) for i,(a,b) in enumerate(zip(o,n)) if F(a)+F(b)}) for o,n in zip(d['actual_old_two_by619'],d['actual_native_two_by619'])]
 norms={};reprs={};original=[];selected=[];factors=[];receipts=[]
 for j,v0 in enumerate(cols):
  v=dict(v0);a={};steps=[]
  while v and min(v) in norms:
   p=min(v);c=v[p];v=H.guard(add(v,norms[p],-c));a=H.guard(add(a,reprs[p],c));steps.append([p,str(c)])
  if v:
   p=min(v);c=v[p];i=len(original);original.append(v0);selected.append(j);norms[p]=H.guard(scale(v,1/c));reprs[p]=H.guard(scale(add({i:F(1)},a,-1),1/c));a={i:F(1)};receipts.append({'column':j,'pivot':p,'divide_by':str(c),'subtract':steps})
  assert mv(original,a)==v0;factors.append(a)
 rank=len(original);minor=None
 if rank==1:minor={'rows':[next(iter(original[0]))],'columns':selected,'determinant':str(original[0][next(iter(original[0]))])}
 if rank==2:minor={'rows':[0,1],'columns':selected,'determinant':str(original[0].get(0,F(0))*original[1].get(1,F(0))-original[0].get(1,F(0))*original[1].get(0,F(0)))};assert F(minor['determinant'])
 hf=[6,6-rank];ep=out/'actual-two-by619-rank-and-complete-first-Hilbert-layer.json';H.dump(ep,{'all619_actual_two_class_columns':[srec(v) for v in cols],
  'rank':rank,'fixed_rank_receipts':receipts,'fixed_original_nonzero_minor':minor,'all619_original_basis_factorizations':[srec(v) for v in factors],
  'selected_original_columns':[srec(v) for v in original],'known_full_first_relation_rank':48+rank,
  'actual_HF_T2_through1':hf,'actual_length_T2_mod_m2':sum(hf),
  'full_T2_length_or_H0_length_or_K_exponents_not_claimed':True,'new_physical_implicit_or_gauge_coefficients':False})
 count=4*619+sum(map(len,cols))+sum(map(len,factors))+3000;assert count<g['whole_slot_bound']<150000
 report={'id':STEM,'all_passed':True,'actual_two_class_rank':rank,'actual_HF_T2_through1':hf,'actual_length_T2_mod_m2':sum(hf),
  'T2_total_length_or_H0_K_exponents_computed':False,'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'owned_QQ_slot_bound':count,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','actual_two_class_rank','actual_HF_T2_through1','actual_length_T2_mod_m2','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate');p.add_argument('--native-peer');a=p.parse_args()
 if a.gate:values(a.output.resolve(),a.gate,a.native_peer)
 else:gate(a.output.resolve())

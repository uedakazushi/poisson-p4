#!/usr/bin/env python3
"""Coefficient-free direct-third support on the617 lower-kernel recipes."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_SUPPORT_PREFLIGHT_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 base='reports/AMBIENT_B18_ELL2_T2_LOWER_RANK2_HF_TAIL_SOURCE_20261105_v63/';r=read(base+'report.json');assert hs[base+'report.json']=='b5c0ec061ee48c403354e1e9f5f40b7906fb738ce49b2d18cc92048e1dba8cf4' and r['all_passed']
 recipes=field(base+'complete1012-kernel617-and-filtered-r-tail-source-proof.json','all617_kernel_lower_constant_source_recipes');assert len(recipes)==617 and {v['free_case619'] for v in recipes}==set(range(619))-{563,609}
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER619_SOURCE_CLEARING_20261105_v63_REVIEW/report.json')['all_passed']
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER_FREE3_GAUGE_SOURCE_20261105_v63_REVIEW/report.json')['all_passed']
 original=field('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json','all_zero_bottom_and_Y0_source_combinations');C=[v['original768_combination'] for v in original];del original
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';B=[dict(v) for v in field(mug,'old640_source_columns')];ops=dict(field(mug,'old640_typed_closure_nodes'));P=field(H.PATHS['C0'],'honest_D0_relation_columns');G=field('reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json','Gamma_two');rows={i for g in G for i,x in g}
 cases=[];needed=set();opaque=set();oldids=set();nids=set()
 for rec in recipes:
  factors=[[rec['free_case619'],'1']]+[[j,'-'+x] for j,x in zip([563,609],rec['subtract_actual_basis563609'],strict=True) if F(x)]
  # Preserve exact subtraction as a sign/operator; never concatenate a
  # leading minus into an already-negative rational for later evaluation.
  factors=[{'case619':rec['free_case619'],'scale_operator':['literal','1']}]+[{'case619':j,'scale_operator':['negate_literal',x]} for j,x in zip([563,609],rec['subtract_actual_basis563609'],strict=True) if F(x)]
  old={i for f in factors for i,x in C[f['case619']] if i<640};native={i-640 for f in factors for i,x in C[f['case619']] if i>=640};typed=set();literal=set()
  for j in old:
   for row in rows:
    if row<48:
     k=B[j].get(row)
     if k is None or ops[k][0]=='c':continue
     if ops[k][0]=='actual_C0_scalar_function':literal.add(ops[k][1])
     else:typed.add(k)
    else:literal.add(P[j][row-48+1])
  cases.append({'free_case619':rec['free_case619'],'constant_case619_source_factors':factors,'original_old640_support':sorted(old),'original_native128_support':sorted(native),'typed_direct_cubic_roots':sorted(typed),'literal_C0_direct_cubic_roots':sorted(literal)})
  needed.update(typed);opaque.update(literal);oldids.update(old);nids.update(native)
 closure=set()
 def visit(i):
  if i in closure:return
  closure.add(i);n=ops[i]
  for j in ([j for j,a,b in n[1]] if n[0]=='sum' else n[3] if n[0]=='product' else []):visit(j)
 for i in needed:visit(i)
 deps={}
 for i in sorted(closure):
  n=ops[i]
  if n[0]=='actual_C0_scalar_function':v={n[1]}
  elif n[0]=='c':v=set()
  else:v=set().union(*(deps[j] for j in ([j for j,a,b in n[1]] if n[0]=='sum' else n[3])))
  deps[i]=v
 gradients=0;percase=[]
 for c in cases:
  labels=set(c['literal_C0_direct_cubic_roots'])
  for i in c['typed_direct_cubic_roots']:labels.update(deps[i])
  percase.append([c['free_case619'],len(labels)]);gradients+=2*len(labels)
 inputdata={'all617_source_and_direct_third_support_recipes':cases,'original619_768_combinations':C,'Gamma_two':G,'pruned_typed_source_closure':[[i,ops[i]] for i in sorted(closure)],'additional_literal_C0_third_symbols':sorted(opaque),'no_original_source_scalar_factors_folded':True,'no_typed_AD_or_scalar_values':True,
 'b1_A2_M1_and_native_carries_remain_typed':'UseacceptedexactUPPER619_SOURCE_CLEARINGlinearlyonconstantrecipes. All918b1domain;noassociatedgradedexactness. ComplementaryM1setonlybyacceptedcoherentfree3existenceframe,notcanonicalzerofill.',
 'basis_upper_entries_not_zero':'SubtracttwoactualsourcesFIRST;newuppercoordinatesonkerLOWERareactualsourcevectors,notbasisupperreadbacks.',
 'direct_third_potential_gradient_entries_two_channels':gradients,'per_case_twochannel_scalar_dependency_upper':percase,'future_higher_jet_necessity_not_claimed':True}
 retained=count([C,B,ops,P,cases,deps,G])+7000
 future_input=count(inputdata);future=2*(future_input+1000)+4*sum(map(len,deps.values()))+4*gradients+12000
 ep=out/'factored617-lower-kernel-direct-third-support-and-source-frontier.json';H.dump(ep,inputdata)
 report={'id':STEM,'all_passed':True,'coefficient_free_only':True,'all617_constant_sources_factored_before_AD':True,'typed_nodes':len(closure),'literal_C0_third_symbol_candidates':len(opaque),'direct_old_columns':len(oldids),'direct_native_columns':len(nids),'direct_gradient_support_upper_twochannels':gradients,'retained_numeric_bound':retained,'prospective_naive_formal_ownership_upper':future,'future_formal_height_RAW_and_ownership_gate_pending':True,'no_upper_values_or_new_source_coefficient':True,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};assert retained<150000;H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','typed_nodes','literal_C0_third_symbol_candidates','direct_gradient_support_upper_twochannels','retained_numeric_bound','prospective_naive_formal_ownership_upper']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)

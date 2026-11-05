#!/usr/bin/env python3
"""Coefficient-free eighteen-psi support/ancestor inventory; no source AD."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
from ambient_B18_ell2_raw33_quotient_20261105_v63 import srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_NONACTIVE9_UPPER_SUPPORT_PREFLIGHT_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 src=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE171_SOURCE_PREFLIGHT_20261105_v63/exact171-decomposition-162-zero-9-lower-first-source-input.json')
 assert read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_LOWER_VALUES_20261105_v63_REVIEW/report.json')['all_passed']
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'
 mu='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json'
 combs=field(mug,'all741_zero_origin_source_combinations');sel=field(mu,'actual_selected_Y0_source_columns');C=[H.parse(combs[j]['original768']) for j in sel];del combs
 ops=dict(field(mug,'old640_typed_closure_nodes'));B=[dict(v) for v in field(mug,'old640_source_columns')];P=field(H.PATHS['C0'],'honest_D0_relation_columns')
 prev='reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63/actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json'
 Gamma=list(map(H.parse,field(prev,'Gamma56_two_channels')));lam=list(map(H.parse,field(prev,'lambda26_two_channels')))
 retained=count([src,C,ops,B,P,Gamma,lam])+6000;assert retained<60000
 H.dump(out/'phase.json',{'phase':'support_only_before_formal_eighteen_weights','retained_source_numeric_bound':retained,'no_scalar_source_combinations_or_AD_yet':True})
 rows=[];needed=set();opaque=set();usedj=set()
 for r in src['remaining9_sources']:
  br=H.parse(r['b1_r_selected99']);be=H.parse(r['b1_ell_selected99'])
  for name,w in [('A',{2*j:x for j,x in br.items()}),('B',{**{2*j:x for j,x in be.items()},**{2*j+1:x for j,x in br.items()}})]:
   source=[set(),set()];old=set()
   for k,x in w.items():
    j,ch=divmod(k,2);usedj.add(j);source[ch].update(C[j])
   for ch,ss in enumerate(source):
    for j in ss:
     if j>=640:continue
     for row in Gamma[ch]:
      if row<48:
       root=B[j].get(row)
       if root is None:continue
       n=ops[root]
       if n[0]=='c':continue
       key=(1,n[1]) if n[0]=='actual_C0_scalar_function' else (0,root)
      else:key=(1,P[j][row-48+1])
      old.add(key)
   needed.update(i for kind,i in old if kind==0);opaque.update(i for kind,i in old if kind==1)
   rows.append({'original418_source':r['source_index'],'contraction':name,'formal_weights_on99psi2':srec(w),'original768_support_two_channels':[sorted(v) for v in source],
    'possible_old_typed_or_literal_root_labels':[list(v) for v in sorted(old)],'native128_support_two_channels':[sorted(j-640 for j in v if j>=640) for v in source]})
 assert len(rows)==18
 closure=set()
 def visit(i):
  if i in closure:return
  closure.add(i);n=ops[i]
  if n[0]=='sum':
   for j,a,b in n[1]:visit(j)
  elif n[0]=='product':
   for j in n[3]:visit(j)
  else:assert n[0] in ['c','actual_C0_scalar_function']
 for i in needed:visit(i)
 nodes=[[i,ops[i]] for i in sorted(closure)];atoms=sorted({n[1] for i,n in nodes if n[0]=='actual_C0_scalar_function'})
 # Structural unknown-second closure only: no centers/first values are used,
 # and no cancellation or cache availability is inferred from the labels.
 deps={}
 for i,n in nodes:
  if n[0]=='actual_C0_scalar_function':v={n[1]}
  elif n[0]=='c':v=set()
  else:
   kids=[j for j,a,b in n[1]] if n[0]=='sum' else n[3];assert all(j in deps for j in kids);v=set().union(*(deps[j] for j in kids))
  deps[i]=v
 for r in rows:
  labels=set()
  for kind,i in r['possible_old_typed_or_literal_root_labels']:labels.update(deps[i] if kind==0 else {i})
  r['uncancelled_structural_original_second_labels']=sorted(labels)
 # A prospective constant-root fold cannot have more coefficients than these
 # support sets. M1/pivot/Y0 responses remain factored, not expanded matrices.
 source_entries=sum(sum(len(v) for v in r['original768_support_two_channels']) for r in rows)
 roots=sum(len(r['possible_old_typed_or_literal_root_labels']) for r in rows)
 native_entries=sum(sum(len(v) for v in r['native128_support_two_channels']) for r in rows)
 gradients=sum(len(r['uncancelled_structural_original_second_labels']) for r in rows)
 folding_upper=retained+count(rows)+count(nodes)+2*source_entries+2*roots+4*len(usedj)*56+8000
 ad_upper=count(rows)+count(nodes)+count(deps)+4*len(nodes)+sum(map(len,deps.values()))+2*gradients+4*len(atoms)+8000
 assert retained+count(rows)+count(nodes)+count(deps)<60000
 data={'eighteen_formal_contractions':rows,'selected99_source_indices_used':sorted(usedj),'typed_closure_nodes':nodes,'possible_atom_second_labels':atoms,'additional_literal_second_labels':sorted(opaque),
  'no_coefficient_source_folding_or_AD_or_values_executed':True,'positive_supports_ignore_all_signed_cancellation':True,
  'scope':'Source-only forecast. All honestsecond labels remain independent; no ancestor value/cache/H20 descent claimed. Native L1r/OmegaEll andL1ell/OmegaR must remain. M1andpivot terms remain factored; A0cancellation requires separateplaneW2 theorempeer.',
  'prospective_constant_source_entry_support':source_entries,'prospective_old_root_support':roots,'prospective_native_entry_support':native_entries,'prospective_structural_gradient_entries':gradients,
  'prospective_source_fold_ownership_upper':folding_upper,'prospective_formal_AD_ownership_upper':ad_upper,'source_fold_height_and_RAW_gate_still_required':True}
 ep=out/'eighteen-formal-weights-and-coefficient-free-source-ancestor-inventory.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'source_support_only':True,'formal_scalar_weights':18,'typed_nodes':len(nodes),'typed_atom_labels':len(atoms),'possible_literal_second_labels':len(opaque),'used_selected99':len(usedj),'prospective_source_fold_ownership_upper':folding_upper,'prospective_formal_AD_ownership_upper':ad_upper,'scalar_fold_AD_or_actual_values_executed':False,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','typed_nodes','typed_atom_labels','possible_literal_second_labels','used_selected99','prospective_source_fold_ownership_upper','prospective_formal_AD_ownership_upper']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

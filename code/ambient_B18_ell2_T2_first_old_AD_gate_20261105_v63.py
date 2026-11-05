#!/usr/bin/env python3
"""Exact sparse-support forecast for the formal old rell second AD."""
import argparse,json,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_FIRST_OLD_AD_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 f=read('reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63/actual-two-rell-source-folds-and-retained-native-first-operators.json')
 a=read('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json');origins={i:H.parse(v).get(0,F(0)) for i,v in a['all79_atom_values']}
 values={};deps={}
 for i,n in f['needed_typed_closure_nodes']:
  if n[0]=='c':v=F(n[1],n[2]);s=set()
  elif n[0]=='actual_C0_scalar_function':v=origins[n[1]];s={n[1]}
  elif n[0]=='sum':
   v=sum((F(a,b)*values[j] for j,a,b in n[1]),F(0));s=set().union(*(deps[j] for j,a,b in n[1] if a))
  else:
   assert n[0]=='product';js=n[3];v=F(n[1],n[2]);s=set()
   for j in js:v*=values[j]
   for k,j in enumerate(js):
    c=F(n[1],n[2])
    for l,u in enumerate(js):
     if l!=k:c*=values[u]
    if c:s.update(deps[j])
  values[i]=v;deps[i]=s
 rowbound=0;P=set();active=0
 for row in f['all619_source_folds']:
  for func in row['two_old_typedB_or_originalP_rell_root_functionals']:
   symbols=set()
   for (kind,i),c in func:
    if kind:symbols.add(i);P.add(i)
    else:symbols.update(deps[i])
   rowbound+=len(symbols);active+=bool(func)
 slots={'frozen_source_fold_records':f['folded_coefficient_count'],
  'typed_constant_r_ell_known_rell_slots':4*len(values),'exact_structural_gradient_upper':sum(map(len,deps.values())),
  'all1248_source_second_sensitivity_upper':rowbound,'cached_atomfirst_and_original84closure':79*10+84*4,
  'coherent_P_source_root_catalog_and_metadata_reserve':12000}
 total=sum(slots.values());assert total<150000
 ep=out/'actual-old-rell-AD-sparse-dependencies-and-budget.json'
 H.dump(ep,{'typed_node_constants':[[i,str(values[i])] for i in sorted(values)],
  'typed_source_second_dependency_labels':[[i,sorted(deps[i])] for i in sorted(deps)],
  'original_P_or_direct_C0_second_labels':sorted(P),'active_old_functionals':active,
  'next_stage_unknowns':'all originalC0 honest second_rell symbols retained; coherent negativealias/copy ancestor transport requiredbeforeanyclaimvalues',
  'next_AD_formula':'product_rell=a0*b_rell+b0*a_rell+a_r*b_ell+a_ell*b_r; ordinarycoefficient nofactor2',
  'new_function_values_or_source_jets':False,'AD_started':False,'whole_sparse_AD_slot_bound':total,'slot_breakdown':slots,
  'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'QQ_coefficients':150000,'bits':384}})
 report={'id':STEM,'all_passed':True,'whole_sparse_AD_slot_bound':total,'active_old_functionals':active,
  'typed_nodes':len(values),'structural_gradient_upper':sum(map(len,deps.values())),
  'all_source_sensitivity_upper':rowbound,'P_direct_label_count':len(P),'AD_started':False,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','whole_sparse_AD_slot_bound','structural_gradient_upper','all_source_sensitivity_upper','P_direct_label_count','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

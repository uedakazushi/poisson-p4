#!/usr/bin/env python3
"""Formal rell second AD; honest second coefficients stay literal symbols."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_FIRST_OLD_AD_20261105_v63'
def main(out,gate):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'authorizedformaloldAD_no_secondvalues'})
 hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 gr=read(gate);assert gr['all_passed'] and not gr['AD_started'];ip=str(Path(gate).parent/'actual-old-rell-AD-sparse-dependencies-and-budget.json');budget=read(ip);assert hs[ip]==gr['evidence_sha256'][Path(ip).name]
 source=read('reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63/actual-two-rell-source-folds-and-retained-native-first-operators.json')
 first=read('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json');old={i:H.parse(p) for i,p in first['all79_atom_values']}
 d=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json');memo={b['original_node']:H.parse(b['origin_and_ds']) for b in d['all_exact_cached_leaf_bindings_origin_and_ds']}
 for i,n in d['old_exact84_original_arithmetic_nodes']:
  if n[0]=='c':v={0:F(n[1],n[2])} if n[1] else {}
  elif n[0]=='a':v=add(memo[n[1]],memo[n[2]])
  elif n[0]=='n':v=scale(memo[n[1]],-1)
  else:v=H.pairmul(memo[n[1]],memo[n[2]])
  memo[i]=H.guard(v)
 atoms={i:[old[i].get(0,F(0)),old[i].get(8,F(0)),memo[i].get(1,F(0))/5,F(0)] for i in source['needed_original_oldB_atom_labels']}
 vals={};grads={};orig=dict((i,F(c)) for i,c in budget['typed_node_constants']);deps=dict(budget['typed_source_second_dependency_labels'])
 def plus(a,b):return [x+y for x,y in zip(a,b)]
 def times(a,b):return [a[0]*b[0],a[0]*b[1]+b[0]*a[1],a[0]*b[2]+b[0]*a[2],a[0]*b[3]+b[0]*a[3]+a[1]*b[2]+a[2]*b[1]]
 for i,n in source['needed_typed_closure_nodes']:
  if n[0]=='c':v=[F(n[1],n[2]),F(0),F(0),F(0)];g={}
  elif n[0]=='actual_C0_scalar_function':v=atoms[n[1]];g={n[1]:F(1)}
  elif n[0]=='sum':
   v=[F(0)]*4;g={}
   for j,a,b in n[1]:c=F(a,b);v=plus(v,[c*x for x in vals[j]]);g=add(g,grads[j],c)
  else:
   assert n[0]=='product';v=[F(n[1],n[2]),F(0),F(0),F(0)];g={}
   for j in n[3]:g=add(scale(g,vals[j][0]),grads[j],v[0]);v=times(v,vals[j])
  H.guard({k:c for k,c in enumerate(v) if c});H.guard(g);assert v[0]==orig[i] and set(g)<=set(deps[i]);vals[i]=v;grads[i]=g
 results=[];support=0
 for row in source['all619_source_folds']:
  parts=[]
  for channel,func in enumerate(row['two_old_typedB_or_originalP_rell_root_functionals']):
   known=F(row['two_known_oldpivot_and_L0S1_corrections'][channel]);unknown={}
   for (kind,i),coefficient in func:
    c=F(coefficient)
    if kind:unknown=add(unknown,{i:c})
    else:known+=c*vals[i][3];unknown=add(unknown,grads[i],c)
   H.guard(unknown);H.guard({0:known} if known else {});support+=len(unknown)+bool(known)
   parts.append({'known_lower_first_response_only':str(known),'literal_original_C0_second_rell_sensitivity':srec(unknown)})
  results.append({'index':row['index'],'two_formal_old_functionals':parts})
 total=support+sum(map(len,grads.values()))+4*len(vals)+len(memo)*2+source['folded_coefficient_count']+10000;assert total<=budget['whole_sparse_AD_slot_bound']<150000
 ep=out/'actual-two-by619-old-formal-rell-second-sensitivities.json';H.dump(ep,{'all619_results':results,
  'typed_node_known4_coefficients':[[i,list(map(str,v))] for i,v in sorted(vals.items())],
  'typed_node_original_second_gradients':[[i,srec(g)] for i,g in sorted(grads.items())],
  'all_unavailable_honest_second_coefficients_retained':True,
  'native_terms_unchanged_external_source_fold':True,'no_actual_second_values_or_HF_claim':True,
  'unknown_second_label_list':sorted({i for r in results for p in r['two_formal_old_functionals'] for i,c in p['literal_original_C0_second_rell_sensitivity']}),
  'owned_QQ_slot_bound':total})
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'typed_nodes':len(vals),'formal_sensitivity_nonzeros':sum(len(p['literal_original_C0_second_rell_sensitivity']) for r in results for p in r['two_formal_old_functionals']),
  'unknown_second_label_count':len({i for r in results for p in r['two_formal_old_functionals'] for i,c in p['literal_original_C0_second_rell_sensitivity']}),
  'actual_second_values_or_HF_computed':False,'new_physical_implicit_gauge_coefficients':False,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'owned_QQ_slot_bound':total,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','formal_sensitivity_nonzeros','unknown_second_label_count','owned_QQ_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);a=p.parse_args();main(a.output.resolve(),a.gate)

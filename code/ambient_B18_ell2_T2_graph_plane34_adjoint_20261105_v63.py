#!/usr/bin/env python3
"""Known12-output plane adjoints, retaining all unknown F2 nuisance."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_GRAPH_PLANE34_ADJOINT_20261105_v63'
def solve(cols,targets):
 norm={};reprs={};receipts=[]
 for j,orig in enumerate(cols):
  v=dict(orig);a={j:F(1)};steps=[]
  while v and min(v) in norm:
   i=min(v);c=v[i];v=H.guard(add(v,norm[i],-c));a=H.guard(add(a,reprs[i],-c));steps.append([i,str(c)])
  if v:
   i=min(v);c=v[i];norm[i]=H.guard(scale(v,1/c));reprs[i]=H.guard(scale(a,1/c));assert mv(cols,reprs[i])==norm[i];receipts.append({'source_column':j,'pivot':i,'divide_by':str(c),'subtract':steps})
 results=[]
 for k,target in enumerate(targets):
  v=dict(target);rho={}
  while v and min(v) in norm:
   i=min(v);c=v[i];v=H.guard(add(v,norm[i],-c));rho=H.guard(add(rho,reprs[i],c))
  if not v:assert mv(cols,rho)==target;results.append({'invariant':k,'representation_exists':True,'rho':srec(rho),'all104_coefficients_F1_and_F2_exact':True})
  else:
   i=min(v);dual={i:F(1)}
   for q in sorted(norm,reverse=True):
    c=-sum((a*dual.get(j,F(0)) for j,a in norm[q].items() if j!=q),F(0))
    if c:dual[q]=c
   assert all(sum((a*dual.get(j,F(0)) for j,a in col.items()),F(0))==0 for col in cols);value=sum((a*dual.get(j,F(0)) for j,a in target.items()),F(0));assert value
   results.append({'invariant':k,'representation_exists':False,'obstruction104':srec(H.guard(dual)),'nonzero_value':str(value)})
 return results,receipts,norm,reprs
def main(out,gate,sourcepeer):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 p=ROOT/gate;g=json.loads(p.read_text());assert g['all_passed'] and not g['twelve_adjoint_solve_started'];ip=p.parent/'actual-plane-G1-degree23-and-twelve-adjoint-source-input.json';assert H.sha(ip)==g['evidence_sha256'][ip.name];d=json.loads(ip.read_text());pp=ROOT/sourcepeer;peer=json.loads(pp.read_text());assert peer['all_passed']
 cols=list(map(H.parse,d['twelve_source_adjoint_columns104']));targets=list(map(H.parse,d['two_invariant_functionals_on104_with_F2_zero']));first,_,_,_=solve(cols[:6],targets)
 if all(r['representation_exists'] for r in first):results,receipts,norm,reprs=solve(cols[:6],targets);used=6
 else:results,receipts,norm,reprs=solve(cols,targets);used=12
 for r in results:
  if r['representation_exists']:assert mv(cols,H.parse(r['rho']))==targets[r['invariant']]
 ep=out/'actual-two-plane-source-adjoint-certificates.json';H.dump(ep,{'degree3_only_results':first,'complete_results':results,'source_columns_used':used,
  'fixed_receipts':receipts,'normalized_basis':[[i,srec(v)] for i,v in sorted(norm.items())],
  'normalized_basis_source_representations':[[i,srec(v)] for i,v in sorted(reprs.items())],
  'all52_unknown_F2_nuisance_coordinates_retained':True,'actual_invariant_Phi3_Phi4_M1_M2_values_not_computed':True})
 count=sum(map(len,cols+targets))+sum(map(len,norm.values()))+sum(map(len,reprs.values()))+6000;assert count<d['whole_twelve_adjoint_slot_bound']<20000
 report={'id':STEM,'all_passed':True,'degree3_only_representation_exists':[r['representation_exists'] for r in first],
  'degree34_representation_exists':[r['representation_exists'] for r in results],'source_columns_used':used,
  'actual_Phi3_Phi4_M1_M2_or_invariant_values_computed':False,'input_sha256':{gate:H.sha(p),str(ip.relative_to(ROOT)):H.sha(ip),sourcepeer:H.sha(pp)},
  'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:H.sha(ep)},
  'owned_QQ_slot_bound':count,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','degree3_only_representation_exists','degree34_representation_exists','source_columns_used','owned_QQ_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);p.add_argument('--source-peer',required=True);a=p.parse_args();main(a.output.resolve(),a.gate,a.source_peer)

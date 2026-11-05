#!/usr/bin/env python3
"""Known mixed-source adjoints; no actual new physical coefficients."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_graph_plane34_adjoint_20261105_v63 import solve
from ambient_B18_ell2_raw33_quotient_20261105_v63 import mv,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_GRAPH_MIXED_ADJOINT_20261105_v63'
def main(out,gate,sourcepeer):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'authorized_known186adjoint_no_sourcevalues'})
 p=ROOT/gate;g=json.loads(p.read_text());assert g['all_passed'] and not g['adjoint_solve_started'];ip=p.parent/'actual-full-G1-and186-mixed-adjoint-source-input.json';assert H.sha(ip)==g['evidence_sha256'][ip.name];d=json.loads(ip.read_text());pp=ROOT/sourcepeer;peer=json.loads(pp.read_text());assert peer['all_passed']
 cols=list(map(H.parse,d['actual186_source_adjoint_columns286']));targets=list(map(H.parse,d['two_invariant_functionals286']))
 results,receipts,norm,reprs=solve(cols[:180],targets);first=results
 if all(r['representation_exists'] for r in first):used=180
 else:
  del receipts,norm,reprs
  results,receipts,norm,reprs=solve(cols,targets);used=186
 for r in first+results:
  r.pop('all104_coefficients_F1_and_F2_exact',None)
  if r['representation_exists']:r['all286_F1_F2_coefficients_exact']=True
  elif 'obstruction104' in r:r['obstruction286']=r.pop('obstruction104')
 for r in results:
  if r['representation_exists']:assert mv(cols,H.parse(r['rho']))==targets[r['invariant']]
 ep=out/'actual-two-full286-mixed-source-adjoint-certificates.json';H.dump(ep,{'degree23_only_results':first,'complete_results':results,'source_columns_used':used,
  'output186_order':d['output186_order'],'fixed_receipts':receipts,
  'normalized_basis':[[i,srec(v)] for i,v in sorted(norm.items())],
  'normalized_basis_source_representations':[[i,srec(v)] for i,v in sorted(reprs.items())],
  'all234_F1_and52_F2_nuisance_coordinates_retained':True,'actual_new_Phi3_Phi4_invariant_M1_M2_values_not_computed':True})
 count=sum(map(len,cols+targets))+sum(map(len,norm.values()))+sum(map(len,reprs.values()))+sum(len(r['subtract']) for r in receipts)+10000;assert count<d['whole_known_array_slot_bound']<=150000
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'degree23_representation_exists':[r['representation_exists'] for r in first],
  'degree234_representation_exists':[r['representation_exists'] for r in results],'source_columns_used':used,
  'actual_new_source_invariant_M1_M2_values_computed':False,
  'input_sha256':{gate:H.sha(p),str(ip.relative_to(ROOT)):H.sha(ip),sourcepeer:H.sha(pp)},
  'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'owned_QQ_slot_bound':count,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','degree23_representation_exists','degree234_representation_exists','source_columns_used','owned_QQ_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);p.add_argument('--source-peer',required=True);a=p.parse_args();main(a.output.resolve(),a.gate,a.source_peer)

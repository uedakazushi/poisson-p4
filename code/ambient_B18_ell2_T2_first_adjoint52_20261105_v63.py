#!/usr/bin/env python3
"""Point-free source adjoints for the two remaining T2 linear classes.

Only known G2-linear arrays; no actual X2 values or source coefficients.
"""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=H.ROOT;GATE='AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_INPUT_GATE_20261105_v63';STEM='AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63'
def parse(v):return H.parse(v)
def preflight(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 d=read('reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json')
 r=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json');assert r['rank']==25
 g=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json')
 G2=list(map(parse,g['full_G2_linear234_by102_columns']));PiY=list(map(parse,d['PiY3_from102_columns']));Pi2=list(map(parse,r['Pi_remaining_top_linear_columns']))
 psi=[{} for _ in range(2)]
 for i in range(102):
  for v in range(9):
   for t,a in PiY[i].items():
    for k,b in Pi2[t*9+v].items():psi[k]=add(psi[k],{i*9+v:a*b})
 pairs=[(7,7),(7,8)];cols=[{} for _ in range(52)]
 for j,col in enumerate(G2):
  for k,c in col.items():
   row,v=divmod(k,9)
   for w in range(9):
    p=tuple(sorted((v,w)))
    if p in pairs:cols[2*row+pairs.index(p)]=add(cols[2*row+pairs.index(p)],{j*9+w:c})
 count=sum(map(len,cols+psi))+sum(map(len,PiY+Pi2+G2));bound=count+3*52*52+2*918+6000;assert bound<20000
 evidence={'source_domain_order':'Y1[source102,variable9], flat9*source+variable; ordinarycoefficients, ell²=0',
  'target_pairs_only':[['rA','rA'],['rA','ell']],'lambda52_order':'row26major,pair(rr,rell)minor',
  'actual_two_PiY_Pi2_composed_functionals918':[srec(v) for v in psi],
  'actual52_adjoint_columns_in918_domain':[srec(v) for v in cols],
  'all_full_G2linear_source_columns':[srec(v) for v in G2],
  'source_theorem':'Forzero-bottom642kernel, Z=Xgamma-Xp*a1 startsdegree2. ActualY0=0. If lambda*G2linear=psi onALL918 domaincoords, then projectedY1=lambda*M0*(Xgamma2-Xp1*a1); associatedgradedexactnessisnotassumed.',
  'whole_slot_bound':bound,'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'coefficients':150000,'bits':384},
  'rank_or_adjoint_solve_started':False,'physical_X2_M1_or_gauge_queries':False}
 ep=out/'actual52-source-adjoint-input-and-two-full918-functionals.json';H.dump(ep,evidence)
 report={'id':GATE,'all_passed':True,'whole_slot_bound':bound,'rank_or_adjoint_solve_started':False,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','whole_slot_bound','seconds']}))
def solve(out,gate):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'authorized52knownarrayadjoint_no_X2'})
 p=ROOT/gate;g=json.loads(p.read_text());assert g['all_passed'] and not g['rank_or_adjoint_solve_started'];ip=p.parent/'actual52-source-adjoint-input-and-two-full918-functionals.json';assert H.sha(ip)==g['evidence_sha256'][ip.name];d=json.loads(ip.read_text())
 cols=list(map(parse,d['actual52_adjoint_columns_in918_domain']));psi=list(map(parse,d['actual_two_PiY_Pi2_composed_functionals918']));norm={};reprs={};receipts=[]
 for j,col in enumerate(cols):
  v=dict(col);a={j:F(1)};steps=[]
  while v and min(v) in norm:
   k=min(v);c=v[k];v=H.guard(add(v,norm[k],-c));a=H.guard(add(a,reprs[k],-c));steps.append([k,str(c)])
  if v:
   k=min(v);c=v[k];norm[k]=H.guard(scale(v,1/c));reprs[k]=H.guard(scale(a,1/c));assert mv(cols,reprs[k])==norm[k]
   receipts.append({'source_column':j,'pivot':k,'divide_by':str(c),'subtract':steps})
 results=[]
 for channel,target in enumerate(psi):
  v=dict(target);a={}
  while v and min(v) in norm:
   k=min(v);c=v[k];v=H.guard(add(v,norm[k],-c));a=H.guard(add(a,reprs[k],c))
  if not v:
   assert mv(cols,a)==target;results.append({'channel':channel,'exact_rr_rell_adjoint_exists':True,'lambda52':srec(a),'all918_source_coefficients_verified':True})
  else:
   q=min(v);dual={q:F(1)}
   for k in sorted(norm,reverse=True):
    c=-sum((a*dual.get(i,F(0)) for i,a in norm[k].items() if i!=k),F(0))
    if c:dual[k]=c
   assert all(sum((c*dual.get(i,F(0)) for i,c in col.items()),F(0))==0 for col in cols)
   value=sum((c*dual.get(i,F(0)) for i,c in target.items()),F(0));assert value
   results.append({'channel':channel,'exact_rr_rell_adjoint_exists':False,'restricted_support_nonexistence_dual918':srec(H.guard(dual)),
    'dual_value_on_requested_psi':str(value),'all52_dual_annihilations_verified':True,'scope':'onlyrr/rellsupportrejected; no X2-only/globalnegativeclaim'})
 count=sum(map(len,cols+psi))+sum(map(len,norm.values()))+sum(map(len,reprs.values()))+sum(len(r['subtract']) for r in receipts)+sum(len(r.get('lambda52',r.get('restricted_support_nonexistence_dual918',[]))) for r in results)+5000;assert count<20000
 ep=out/'actual-two-full918-source-adjoint-certificates-or-obstructions.json';H.dump(ep,{'results':results,'fixed_known_adjoint_rank':len(norm),
  'fixed_rank_receipts':receipts,'normalized_basis':[[i,srec(v)] for i,v in sorted(norm.items())],
  'normalized_basis52_source_representations':[[i,srec(v)] for i,v in sorted(reprs.items())],
  'actual_X2_or_M1_values_not_computed':True,'owned_QQ_slot_bound':count})
 report={'id':STEM,'all_passed':True,'both_rr_rell_adjoints_exist':all(r['exact_rr_rell_adjoint_exists'] for r in results),
  'per_channel_exists':[r['exact_rr_rell_adjoint_exists'] for r in results],
  'actual_X2_M1_physical_or_gauge_queries':False,'input_sha256':{gate:H.sha(p),str(ip.relative_to(ROOT)):H.sha(ip)},
  'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'owned_QQ_slot_bound':count,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','both_rr_rell_adjoints_exist','per_channel_exists','owned_QQ_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-gate');a=p.parse_args()
 if a.source_gate:solve(a.output.resolve(),a.source_gate)
 else:preflight(a.output.resolve())

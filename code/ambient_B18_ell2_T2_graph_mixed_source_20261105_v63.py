#!/usr/bin/env python3
"""Full known G1 and mixed low-order adjoint source gate only."""
import argparse,itertools,json,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,srec,mv
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_GRAPH_MIXED_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 graph=read('reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json')
 plane=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_PLANE34_INPUT_GATE_20261105_v63/actual-plane-G1-degree23-and-twelve-adjoint-source-input.json')
 t=read(H.PATHS['transport']);pref=read(H.PATHS['active']);res=read(H.PATHS['resolution']);physical=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')
 A=[list(map(H.parse,v)) for v in physical['C2_coordinate_actions9_sparse_columns']];lifts=[]
 for col in res['physical_lift_m0_columns']:
  p={}
  for i,pol in col:
   for (r,e),c in pol:
    m=[0]*9;m[7]=r;m[8]=e
    if F(c):p[i,tuple(m)]=F(c)
  lifts.append(p)
 full=[]
 for rawj in t['active_F1_free17']:
  block,j=divmod(rawj,14);var=8 if block==0 else pref['active_variables'][block-1];p={}
  for (target,m),c in lifts[j].items():
   n=list(m);n[var]+=1
   if n[8]<2:p[target,tuple(n)]=c
  for i,c in A[var][j].items():p=add(p,lifts[i],-c)
  full.append(p)
 for z in pref['zero_actions']:
  for target in range(3):m=[0]*9;m[z]=1;full.append({(target,tuple(m)):F(1)})
 assert len(full)==26
 G1lin=list(map(H.parse,graph['actual_G1linear3_by26_columns']))
 for j,p in enumerate(full):
  lin={target*9+m.index(1):c for (target,m),c in p.items() if sum(m)==1};assert lin==G1lin[j]
  q={(target,(m[7],m[8])):c for (target,m),c in p.items() if not any(m[:7])};expected={(i,tuple(m)):F(c) for i,m,c in plane['exact_plane_G1_polynomial_columns'][j]};assert q==expected
  assert all(sum(m)<=3 and sum(m[:7])<=1 for target,m in p)
 # All26 polynomial relations annihilate the accepted actualC2 module.
 eta=res['physical_target3_basis_positions']
 for p in full:
  q={}
  for (target,m),c in p.items():
   v={eta[target]:c}
   for k,n in enumerate(m):
    for _ in range(n):v=mv(A[k],v)
   q=add(q,v)
  assert not q
 def mon(word):return tuple(word.count(i) for i in range(9))
 mons2=[mon(w) for w in itertools.combinations_with_replacement(range(9),2) if w!=(8,8)];mons3=[mon((z,7,7)) for z in range(7)]+[mon((z,7,8)) for z in range(7)]+[mon((7,7,7)),mon((7,7,8))];mons4=[mon((7,7,7,7)),mon((7,7,7,8))]
 outputs=[(degree,target,m) for degree,mons in [(2,mons2),(3,mons3),(4,mons4)] for target in range(3) for m in mons];oi={key:i for i,key in enumerate(outputs)};assert len(outputs)==186
 columns=[{} for _ in outputs]
 for source,p in enumerate(full):
  for (target,m),c in p.items():
   for var in range(9):
    n=list(m);n[var]+=1;key=(sum(n),target,tuple(n))
    if key in oi:columns[oi[key]]=add(columns[oi[key]],{source*9+var:c})
   for pair,word in enumerate([(7,7),(7,8)]):
    n=list(m)
    for v in word:n[v]+=1
    key=(sum(n),target,tuple(n))
    if key in oi:columns[oi[key]]=add(columns[oi[key]],{234+2*source+pair:c})
 # Any omitted F2 mixed-z coefficient produces>=twoz atselecteddegree3,
 # andany omittedF3 is multipliedbyG1lin(z), absentfromplane4.
 psi=list(map(H.parse,graph['two_invariant_functionals_on_all234_F1_coefficients']))
 slots={'sparse_G1polynomial':sum(map(len,full)),'sparse186_adjoint_columns':sum(map(len,columns)),
  'max_normalized_rank_basis':186*286,'max_source_representations':186*186,'max_receipt_subtractions':186*185//2,
  'two_targets_and_outputs':2*286+2*186,'constant_action_inputs_and_metadata':12000}
 total=sum(slots.values());assert total<150000
 ep=out/'actual-full-G1-and186-mixed-adjoint-source-input.json';H.dump(ep,{'exact_full_G1_polynomial_columns':[[[i,list(m),str(c)] for (i,m),c in sorted(p.items())] for p in full],
  'all26_actual_C2_module_relation_identities_verified':True,'all26_linear_and_plane_overlaps_exact':True,
  'domain286_order':'F1 source26*all9vars(234),thenF2 source26*(rr,rell)(52)',
  'output186_order':[[d,t,list(m)] for d,t,m in outputs],
  'actual186_source_adjoint_columns286':[srec(v) for v in columns],
  'two_invariant_functionals286':[srec(v) for v in psi],
  'full_nuisance_proof':'Phi2=G1linF1. SelectedPhi3=G1_2F1+G1linF2plane; mixedF2haszsoG1lin(z)gives>=2z, outsidezi*r²/zi*r*ell. PlanePhi4=G1_3F1+G1_2F2plane; G1linF3haszandvanishes. All234F1and52F2nuisancesretained.',
  'source_callback_inventory':'Phi2all44cachedB56. SelectedPhi3 mixed14 ordinarycoeffs requireexactfunctionalcacheoverlapwithselected135/source33beforevalues; pureplane3cachedgenuinecommon105. Phi4r4/r3ellformal; noordinaryzero-fillfromtruncation.',
  'no_adjoint_solve_or_new_source_values':True,'whole_known_array_slot_bound':total,'slot_breakdown':slots,
  'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'QQ_coefficients':150000,'bits':384}})
 report={'id':STEM,'all_passed':True,'source_domain':286,'source_output':186,'adjoint_solve_started':False,'new_source_M1_M2_values':False,
  'whole_known_array_slot_bound':total,'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','source_domain','source_output','whole_known_array_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

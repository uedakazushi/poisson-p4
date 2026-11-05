#!/usr/bin/env python3
"""Known-ring input for final combined rho descent; no Phi values/dot."""
import argparse,itertools,json,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import mv,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 raw=read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json');read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json')
 adj=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_MIXED_ADJOINT_20261105_v63/actual-two-full286-mixed-source-adjoint-certificates.json');read('reports/AMBIENT_B18_ELL2_T2_GRAPH_MIXED_ADJOINT_20261105_v63_REVIEW/report.json')
 assert all(r['representation_exists'] for r in adj['complete_results'])
 A=[list(map(H.parse,c)) for c in raw['coordinate_actions9_sparse_columns']];basis=raw['actual33_basis'];mons=[];NF=[]
 for degree in range(5):
  for word in itertools.combinations_with_replacement(range(9),degree):
   if word.count(8)>1:continue
   m=tuple(word.count(k) for k in range(9));v={0:F(1)}
   for k in word:v=mv(A[k],v)
   H.guard(v);mons.append(m);NF.append(v)
 assert len(mons)==660;mi={m:i for i,m in enumerate(mons)};bm=[]
 for i,(s,r,e) in enumerate(basis):
  m=[0]*9;m[7]=r;m[8]=e
  if s:m[s-1]=1
  bm.append(mi[tuple(m)]);assert NF[bm[-1]]=={i:F(1)}
 weights=[];bars=[]
 for result in adj['complete_results']:
  components=[{} for _ in range(3)]
  for k,c in H.parse(result['rho']).items():
   degree,target,m=adj['output186_order'][k];assert sum(m)==degree;components[target][mi[tuple(m)]]=c
  weights.append(components);bars.append([{i:w.get(j,F(0)) for i,j in enumerate(bm) if w.get(j,F(0))} for w in components])
 # Independently verify the stronger scalar plane9 quotient, without reading
 # a coefficient of any implicit source function.
 pmon=[(r,e) for d in range(5) for e in range(min(1,d)+1) for r in [d-e]];pi={m:i for i,m in enumerate(pmon)};images=[]
 for s,r,e in basis:images.append({pi[r,e]:F(1)} if not s and (r,e) in pi else {})
 def pmul(a,b):
  v={}
  for i,c in a.items():
   for j,d in b.items():
    m=tuple(x+y for x,y in zip(pmon[i],pmon[j]))
    if m in pi:v[pi[m]]=v.get(pi[m],F(0))+c*d
  return {i:c for i,c in v.items() if c}
 plane_ok=True
 for k,action in enumerate(A):
  x={pi[1,0]:F(1)} if k==7 else ({pi[0,1]:F(1)} if k==8 else {})
  for j,v in enumerate(action):plane_ok &= mv(images,v)==pmul(x,images[j])
 slot={'all660_NF_dense_upper':660*33,'all6_functional_identities':2*3*660,
  'raw9_actions_dense_source_upper':9*33*33,'candidate_bars_and_weights':2*3*33+100,
  'source_metadata_and_plane9proof':10000};total=sum(slot.values());assert total<60000
 ep=out/'actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json';H.dump(ep,{'monomials660_degree_le4':[list(m) for m in mons],
  'all660_actual_raw33_normal_forms':[srec(v) for v in NF],'actual33_basis_monomial_indices':bm,
  'two_three_component_original_Taylor_functionals':[[srec(v) for v in w] for w in weights],
  'candidate_two_three_component_raw33_covectors':[[srec(v) for v in w] for w in bars],
  'raw33_to_plane9_ring_quotient_verified':bool(plane_ok),'plane9_monomials':pmon,'raw33_plane9_map_columns':[srec(v) for v in images],
  'source_sufficiency':'Ifbar_rho*NF(mon)=rho(mon) onALL660ordinarymonomials forALL6componentfunctionals, thencombinedinvariantfunctionaldescendsactualraw33. Higherterms>=5 killedbyprovenm5raw33=0; no sourcecoefficient4 guessed or individuallycomputed.',
  'next_values_only_if_factor_pass':'UseactualprimitiveV*Wv inraw33 fromfrozenphysical12(delta,bA,bB),monicsnu; nofreshsourceeval. A23graphcallbacksremainseparate.',
  'factor_test_or_invariant_values_not_started':True,'whole_known_factor_slot_bound':total,'slot_breakdown':slot,
  'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'QQ_coefficients':150000,'bits':384}})
 report={'id':STEM,'all_passed':True,'raw33_to_plane9_quotient_verified':bool(plane_ok),'factor_test_started':False,'actual_invariant_values_computed':False,
  'whole_known_factor_slot_bound':total,'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','raw33_to_plane9_quotient_verified','whole_known_factor_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())

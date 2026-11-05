#!/usr/bin/env python3
"""Independent known-ring NF/factor-input gate; no factor verdict or Phi dot."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,itertools,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','gate':f'reports/{A}/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json',
'raw':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
'raw_peer':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json',
'adjoint':'reports/AMBIENT_B18_ELL2_T2_GRAPH_MIXED_ADJOINT_20261105_v63/actual-two-full286-mixed-source-adjoint-certificates.json',
'adjoint_peer':'reports/AMBIENT_B18_ELL2_T2_GRAPH_MIXED_ADJOINT_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert D['primary']['all_passed'] and D['raw_peer']['all_passed'] and D['adjoint_peer']['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 R=D['raw'];E=D['gate'];actions=[list(map(sp,p)) for p in R['coordinate_actions9_sparse_columns']];mons=[tuple(m) for m in E['monomials660_degree_le4']];expected=[tuple(w.count(i) for i in range(9)) for d in range(5) for w in itertools.combinations_with_replacement(range(9),d) if w.count(8)<2];assert mons==expected and len(mons)==660
 NF=[]
 for m in mons:
  p={0:Q(1)}
  for k in reversed(range(9)):
   for _ in range(m[k]):p=mv(actions[k],p)
  NF.append(p)
 assert list(map(enc,NF))==E['all660_actual_raw33_normal_forms'];index={m:i for i,m in enumerate(mons)};basisindices=[]
 for j,(s,r,e) in enumerate(R['actual33_basis']):
  m=[0]*9;m[7]=r;m[8]=e
  if s:m[s-1]=1
  k=index[tuple(m)];basisindices.append(k);assert NF[k]=={j:Q(1)}
 assert basisindices==E['actual33_basis_monomial_indices']
 weights=[];bars=[]
 for result in D['adjoint']['complete_results']:
  w=[{} for _ in range(3)]
  for k,c in sp(result['rho']).items():
   degree,target,m=D['adjoint']['output186_order'][k];assert sum(m)==degree;w[target][index[tuple(m)]]=c
  weights.append(w);bars.append([{j:p.get(k,Q(0)) for j,k in enumerate(basisindices) if p.get(k,Q(0))} for p in w])
 assert [[enc(p) for p in v] for v in weights]==E['two_three_component_original_Taylor_functionals'] and [[enc(p) for p in v] for v in bars]==E['candidate_two_three_component_raw33_covectors']
 pmon=[(d-e,e) for d in range(5) for e in range(min(1,d)+1)];pi={m:i for i,m in enumerate(pmon)};images=[{pi[r,e]:Q(1)} if not s and (r,e) in pi else {} for s,r,e in R['actual33_basis']]
 def pmul(a,b):
  p={}
  for i,c in a.items():
   for j,d in b.items():
    x=tuple(v+w for v,w in zip(pmon[i],pmon[j]))
    if x in pi:p=add(p,{pi[x]:c*d})
  return p
 assert [list(m) for m in pmon]==E['plane9_monomials'] and list(map(enc,images))==E['raw33_plane9_map_columns']
 for k,M in enumerate(actions):
  v={pi[1,0]:Q(1)} if k==7 else {pi[0,1]:Q(1)} if k==8 else {}
  for j,p in enumerate(M):assert mv(images,p)==pmul(v,images[j])
 for i,(s,r,e) in enumerate(R['actual33_basis']):
  for j in range(33):
   p={j:Q(1)}
   for k in ([s-1] if s else [])+[7]*r+[8]*e:p=mv(actions[k],p)
   assert mv(images,p)==pmul(images[i],images[j])
 assert E['raw33_to_plane9_ring_quotient_verified'] and sum(E['slot_breakdown'].values())==45839==E['whole_known_factor_slot_bound']<150000 and E['factor_test_or_invariant_values_not_started']
 ep=o/'independent660-NF-and-raw33-plane9-factor-input.json';dump(o,ep.name,{'inputs':inputs,'all660_raw33_NF_and_all33_literal_basis_units_verified':True,'all6_candidate_component_covectors_from_exact_mixed_rho_support_verified':True,'all297_coordinate_intertwiners_and1089_plane9_ring_products_verified':True,'whole_known_factor_slot_bound':45839,
'source_factor_sufficiency':'BothoriginalTaylorfunctionalsvanishm5. Actualraw33alsokillsm5. Checkingbar_rho*NF(mon)=rho(mon) forall660ordinarymonomials with ell²0 thereforeprovesfactorizationonALLformalfunctions, componentwise. Candidatebar_rhoisuniquebecause33literalbasismonomialsareindependent.',
'plane9_source_proof':'Settingz7=0 sendsIAintoorder3andIBintoorder2, soIAIBintom5. Thus the actualraw33 quotientcanretainbothplanequarticcoefficients; thisisnowverifiedbyitsentireringtable.',
'no_factor_verdict_Phi_dot_M1_or_invariant_values_computed':True,'no_new_source_implicit_or_gauge_coefficients':True})
 report={'id':STEM,'all_passed':True,'known660_factor_input_and_actual_plane9_quotient_accepted':True,'actual_factor_verdict_or_invariant_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()

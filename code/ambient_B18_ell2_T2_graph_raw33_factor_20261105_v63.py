#!/usr/bin/env python3
"""Universal raw33 factor test for final two rho; no implicit values."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import mv,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63'
def main(out,gate):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 p=ROOT/gate;g=json.loads(p.read_text());assert g['all_passed'] and not g['factor_test_started'];ip=p.parent/'actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json';assert H.sha(ip)==g['evidence_sha256'][ip.name];d=json.loads(ip.read_text());NF=list(map(H.parse,d['all660_actual_raw33_normal_forms']));weights=[list(map(H.parse,v)) for v in d['two_three_component_original_Taylor_functionals']];bars=[list(map(H.parse,v)) for v in d['candidate_two_three_component_raw33_covectors']]
 failures=[];checked=0
 for k in range(2):
  for target in range(3):
   for i,nf in enumerate(NF):
    value=sum((c*nf.get(j,F(0)) for j,c in bars[k][target].items()),F(0));requested=weights[k][target].get(i,F(0));checked+=1
    if value!=requested:failures.append({'invariant':k,'physical_target':target,'monomial':d['monomials660_degree_le4'][i],'bar_on_NF':str(value),'requested':str(requested)})
 assert checked==3960
 rawpath='reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json';raw=json.loads((ROOT/rawpath).read_text());A=[list(map(H.parse,v)) for v in raw['coordinate_actions9_sparse_columns']];basis=raw['actual33_basis'];images=list(map(H.parse,d['raw33_plane9_map_columns']));pm=[tuple(v) for v in d['plane9_monomials']];pi={m:i for i,m in enumerate(pm)}
 def pmul(a,b):
  v={}
  for i,c in a.items():
   for j,e in b.items():
    m=tuple(x+y for x,y in zip(pm[i],pm[j]))
    if m in pi:v[pi[m]]=v.get(pi[m],F(0))+c*e
  return {i:c for i,c in v.items() if c}
 for i,(s,r,e) in enumerate(basis):
  for j in range(33):
   v={j:F(1)}
   for k in ([s-1] if s else [])+[7]*r+[8]*e:v=mv(A[k],v)
   assert mv(images,v)==pmul(images[i],images[j])
 included=set(tuple(v) for v in d['monomials660_degree_le4']);omitted=[]
 for degree in range(5):
  for word in itertools.combinations_with_replacement(range(9),degree):
   m=tuple(word.count(k) for k in range(9))
   if m in included:continue
   assert m[8]>=2;v={0:F(1)}
   for k in word:v=mv(A[k],v)
   assert not v;omitted.append(list(m))
 assert len(omitted)==55 and len(included)==660
 ep=out/'actual3960-rho-factor-identities-or-obstructions.json';H.dump(ep,{'combined_two_rho_factor_raw33':not failures,
  'all3960_monomial_component_checks':checked,'failures':failures,
  'actual_two_three_component_bar_rho33':d['candidate_two_three_component_raw33_covectors'],
  'all1089_raw33_to_plane9_products_verified':True,
  'all715_degree_le4_monomials_accounted':True,'included660_have_ell_degree_le1':True,
  'omitted55_ell_degree_ge2_raw33_zero_and_rho_identically_zero':omitted,
  'formal_allorder_factor_proof':'rho vanishesdegree>=5 andell²; actualraw33m5=0,ell²=0. Equalityonall660basisTaylorcoefficients plusomitted55 impliesbar_rho(Phi modraw33)=rho(Phi) forallformalPhi3vectors.',
  'actual_invariant_values_not_computed':True,'A23_graph_callbacks_not_assumed_zero':True})
 count=sum(map(len,NF))+sum(len(v) for row in bars+weights for v in row)+len(failures)*4+10000;assert count<d['whole_known_factor_slot_bound']<60000
 report={'id':STEM,'all_passed':True,'combined_two_rho_factor_raw33':not failures,'factor_failure_count':len(failures),
  'all3960_identities_checked':True,'all715_monomials_accounted':True,'actual_invariant_values_computed':False,
  'input_sha256':{gate:H.sha(p),str(ip.relative_to(ROOT)):H.sha(ip),rawpath:H.sha(ROOT/rawpath)},
  'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'owned_QQ_slot_bound':count,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','combined_two_rho_factor_raw33','factor_failure_count','owned_QQ_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);a=p.parse_args();main(a.output.resolve(),a.gate)

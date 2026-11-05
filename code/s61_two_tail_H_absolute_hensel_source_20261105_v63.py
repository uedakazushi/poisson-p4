#!/usr/bin/env python3
"""Source input/content proof and finite absolute-Hensel design; NO factor/root/series values."""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(kind,out):
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p,h=None):
  v=sha(p)
  if h:assert v==h
  pins[str(p.relative_to(ROOT))]=v;return p
 peer=ROOT/'reports/S61_TWO_TAIL_H_CHANGED15_ROOT_NEWTON_REVIEW_20261105_v63'
 r=json.loads(pin(peer/'report.json','2dcca597d6c2b49d0f6a51a6ab12a0647ac732366dbefff273365ce3d74353fa').read_text())
 assert r['all_passed'] and r['zero_combined_inconclusive'] and r['all_four_Newton_stages_and_CayleyHamilton_checked']
 for s,h in r['input_sha256'].items():pin(ROOT/s,h)
 for s,h in r['evidence_sha256'].items():pin(peer/s,h)
 pin(ROOT/'research/notes/S61_TWO_TAIL_H_CHANGED15_ROOT_NEWTON_REVIEW_20261105_v63.md')
 evidence={}
 if kind=='source':
  p=pin(ROOT/'reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63/factor-sparse.txt','241488fafa9b45e5353f0b707f75902bc5bc32e85ee093b9f2e0a5f07c7d571c')
  lines=p.read_text().splitlines();assert int(lines[0])==94;H=[tuple(map(int,s.split())) for s in lines[1:]];F={};groups={}
  for a,e,d,v in H:F[a,e]=F.get((a,e),0)+v;groups.setdefault((a,e),[]).append([d,str(v)])
  F={k:v for k,v in F.items() if v};assert max(a for a,e in F)==6 and max(e for a,e in F)==9
  leading={e:v for (a,e),v in F.items() if a==6};assert leading=={0:966796875,1:1160156250}
  P=[F.get((0,e),0) for e in range(10)]
  # Formal integer pseudo-division by6e+5, no new parameter sample.
  rem=P[:];quo=[0]*9
  for k in range(9,0,-1):
   t=rem[k];rem=[6*x for x in rem];quo=[6*x for x in quo];quo[k-1]+=t
   rem[k]-=6*t;rem[k-1]-=5*t;assert rem[k]==0
  R=rem[0];assert all(x==0 for x in rem[1:]) and R and R%32003
  # Literal full coefficient identity 6^9 P=(6e+5)quo+R.
  back=[0]*10
  for e,v in enumerate(quo):back[e]+=5*v;back[e+1]+=6*v
  back[0]+=R;assert back==[6**9*v for v in P]
  content=dict(F_d1_input_only=True,canonical_integer_coefficients=[[[a,e],str(v)] for (a,e),v in sorted(F.items())],
   original_H_term_groups=[dict(a=a,e=e,literal_d_terms=v) for (a,e),v in sorted(groups.items())],
   alpha=17578125,L=[55,66],leading_a6=leading,source_a_degree=6,source_e_degree=9,
   a0=P,pseudodivisor=[5,6],scale=6**9,quotient=quo,remainder=R,remainder_mod32003=R%32003,
   formal_content_identity_verified=True,primitive_in_a_over_Fp_e_and_closure=True,
   no_factor_root_Rabin_Hensel_or_recombination_values_executed=True)
  (out/'fixed-d1-input-and-formal-content.json').write_text(json.dumps(content,indent=2)+'\n');evidence['fixed-d1-input-and-formal-content.json']=sha(out/'fixed-d1-input-and-formal-content.json')
  rep=dict(absolute_Fp2_design=True,one_fixed_d1_and_p32003=True,precision_e=10,degree_bounds_a6_e9=True,
   six_root_seed_verified_before_lifts=True,unordered_proper_subsets=31,LC_allocations=2,full_products=62,
   no_factor_API_completeness_assumption=True,base_product_source_upper=1300000,proposed_hard_products=2000000,
   proposed_hard_owned=24000,proposed_RAW=4194304,proposed_field_temp_bits=30,
   native_worker_not_prepared_or_run=True,fresh_code_and_custody_required=True)
  (out/'finite-six-root-proof-and-budget.json').write_text(json.dumps(rep,indent=2)+'\n');evidence['finite-six-root-proof-and-budget.json']=sha(out/'finite-six-root-proof-and-budget.json')
  pin(ROOT/'research/notes/S61_TWO_TAIL_H_ABSOLUTE_HENSEL_SOURCE_20261105_v63.md')
 pin(Path(__file__));rep=dict(id=out.name,all_passed=True,frozen=True,source_or_publication_only=True,kind=kind,
  changed15_independent_zero_peer_accepted=True,no_factor_or_Hensel_values_executed=True,no_geometric_irreducibility_claim_yet=True,
  input_sha256=pins,evidence_sha256=evidence)
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(dict(SHA256=sha(p),kind=kind,content_remainder=R if kind=='source' else None),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--kind',choices=['peer','source'],required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.kind,a.output.resolve())

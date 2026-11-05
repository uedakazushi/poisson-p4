#!/usr/bin/env python3
"""Publish completed negative fixed Rabin result without any arithmetic replay."""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p,h=None):
  s=sha(p)
  if h:assert s==h
  pins[str(p.relative_to(ROOT))]=s;return p
 base=ROOT/'reports/S61_TWO_TAIL_H_FIXED_RABIN_ENTRY_20261105_v63'
 p=pin(base/'report.json','8628312b3b266d5dd12107db3f996492eb98f0607b781411c21225bc838fee58');r=json.loads(p.read_text())
 assert r['computations_completed'] and not r['all_passed'] and not r['irreducible_over_QQ_a_e_d_verified']
 cert=json.loads(pin(base/'fixed6-Rabin-certificate.json',r['evidence_sha256']['fixed6-Rabin-certificate.json']).read_text())
 assert cert['state']['stage']=='completed' and cert['x_p6_equals_x'] and cert['fixed_specialization_degree']==6
 assert len(cert['Rabin_gcds']['3'])==5 and len(cert['Rabin_gcds']['2'])==7
 for s,h in r['input_sha256'].items():pin(ROOT/s,h)
 pin(ROOT/'reports/S61_TWO_TAIL_H_FIXED_RABIN_CONVERSION_SCOPE_20261105_v63/report.json','e88bb64ab174109c60cddf461dd47291aa720442e310343823c752050d860f4a')
 pin(ROOT/'research/notes/S61_TWO_TAIL_H_FIXED_RABIN_RESULT_20261105_v63.md');pin(Path(__file__))
 result=dict(id=out.name,all_passed=True,frozen=True,publication_only=True,
  completed_negative_fixed_candidate_accepted=True,guard_failure=False,
  no_QQ_reducibility_or_irreducibility_claim=True,no_new_prime_specialization_products_GCD_or_factors=True,
  no_Fitting_globalrank_or_exhaustion_claim=True,input_sha256=pins)
 p=out/'report.json';p.write_text(json.dumps(result,indent=2)+'\n');print(sha(p))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())

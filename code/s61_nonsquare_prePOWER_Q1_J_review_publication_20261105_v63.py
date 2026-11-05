#!/usr/bin/env python3
"""Metadata-only publication of independent Q1/J ROOT d-Horner review."""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_Q1_J_ROOT_DHORNER_REVIEW_PUBLICATION_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st):return json.loads(pin(ROOT/'reports'/st/'report.json').read_text())
 peer=load('S61_NONSQUARE_PREPOWER_Q1_J_ROOT_DHORNER_REVIEW_20261105_v63');assert peer['all_passed']
 primary=load('S61_NONSQUARE_PREPOWER_Q1_J_HASSE_ENTRY_20261105_v63');assert primary['all_passed']and primary['result']['J_valuation']==0
 pin(ROOT/'research/notes/S61_NONSQUARE_PREPOWER_Q1_J_HASSE_RESULT_20261105_v63.md')
 pin(ROOT/'research/notes/S61_NONSQUARE_PREPOWER_Q1_J_ROOT_DHORNER_REVIEW_20261105_v63.md');pin(Path(__file__))
 b=dict(exact_ordinary_Jvaluation0=True,all606_coefficients_numericD5_18_symbolica18=True,
  different_arithmetic='original4811/descendingdHorner; noproducerimports/binomial/rootpowers',
  no_Jquotient_otherfactor_GCD_badlocus_or_complement_claim=True)
 (out/'publication-binding.json').write_text(json.dumps(b,indent=2)+'\n')
 r=dict(id=ID,all_passed=True,frozen=True,publication_only=True,input_sha256=pins,
  evidence_sha256={'publication-binding.json':sha(out/'publication-binding.json')})
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

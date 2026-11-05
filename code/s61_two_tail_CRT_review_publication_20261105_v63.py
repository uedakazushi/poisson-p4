#!/usr/bin/env python3
"""Metadata-only publication of CRT source audit and exact-factor primary; no arithmetic rerun."""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(kind,out):
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st):return json.loads(pin(ROOT/'reports'/st/'report.json').read_text())
 if kind=='source':
  r=load('S61_TWO_TAIL_CRT_LIFT_ROOT_SOURCE_REVIEW_20261105_v63');assert r['all_passed']and r['pre_arithmetic_review_only']
  pin(ROOT/'research/notes/S61_TWO_TAIL_CRT_LIFT_ROOT_SOURCE_REVIEW_20261105_v63.md')
  b=dict(source_code_only=True,actual_value_verification_not_claimed_by_sourceaudit=True,history_ATTEMPT1_2_preserved=True)
 else:
  r=load('S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63');assert r['all_passed']and r['result']['common_char0_factor_proved']and r['result']['verified_originals']==[1,1]
  for n,h in r['evidence_sha256'].items():p=pin(ROOT/'reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63'/n);assert sha(p)==h
  pin(ROOT/'research/notes/S61_TWO_TAIL_CRT_LIFT_RESULT_20261105_v63.md')
  b=dict(primary_exact_both_integer_identities=True,H94_A01795_A2724=True,
   independent_identity_peer_NOT_claimed_here=True,no_maximality_irreducibility_Fitting_or_fullrank_claim=True)
 pin(Path(__file__));(out/'publication-binding.json').write_text(json.dumps(b,indent=2)+'\n')
 report=dict(id=out.name,all_passed=True,frozen=True,publication_only=True,kind=kind,input_sha256=pins,
  evidence_sha256={'publication-binding.json':sha(out/'publication-binding.json')})
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--kind',choices=['source','primary'],required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.kind,a.output.resolve())

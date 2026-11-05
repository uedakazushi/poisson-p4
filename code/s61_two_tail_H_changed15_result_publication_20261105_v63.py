#!/usr/bin/env python3
"""Publish completed-zero changed15 and accepted source audit; no value replay."""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p,h=None):
  v=sha(p)
  if h:assert v==h
  pins[str(p.relative_to(ROOT))]=v;return p
 base=ROOT/'reports/S61_TWO_TAIL_H_CHANGED15_AXIS_ENTRY_20261105_v63'
 r=json.loads(pin(base/'report.json','2b93dce6522f3f08300270b278bbb860f0cb29e7cbfcc31c9a777ed4bed275ee').read_text())
 assert r['all_passed'] and r['computations_completed'] and r['zero_combined_inconclusive']
 assert r['detPair_remainder']==r['combined_remainder']==[0]*4
 for s,h in r['input_sha256'].items():pin(ROOT/s,h)
 for s,h in r['evidence_sha256'].items():pin(base/s,h)
 p=pin(ROOT/'reports/S61_TWO_TAIL_H_CHANGED15_ROOT_SOURCE_REVIEW_20261105_v63/report.json','95d94bfb53418cc6094e8d2bdacdb704496cacf627d98f7389830007306cae74');assert json.loads(p.read_text())['all_passed']
 pin(ROOT/'research/notes/S61_TWO_TAIL_H_CHANGED15_AXIS_RESULT_20261105_v63.md');pin(Path(__file__))
 rep=dict(id=out.name,all_passed=True,frozen=True,publication_only=True,completed_zero_changed15_is_inconclusive=True,
  no_char0_divisibility_or_actual_bad_Fitting_rank_claim=True,no_other_minor_axis_prime_or_retry=True,
  no_new_polynomial_values_or_arithmetic_replay=True,input_sha256=pins)
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n');print(sha(p))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())

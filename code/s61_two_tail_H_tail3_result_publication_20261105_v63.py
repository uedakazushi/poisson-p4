#!/usr/bin/env python3
"""Pin completed inconclusive j3 and zero-kappa state diagnostic; no polynomial values replay."""
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
 base=ROOT/'reports/S61_TWO_TAIL_H_TAIL3_AXIS_ENTRY_20261105_v63'
 r=json.loads(pin(base/'report.json','4a35df1237b235b7acf62bc9eeb7595aeb79de171ecf157bc11763792cddc23b').read_text())
 assert r['all_passed'] and r['computations_completed'] and r['zero_combined_remainder_inconclusive']
 s=json.loads(pin(base/'complete-fixed-source-ring-certificate.json',r['evidence_sha256']['complete-fixed-source-ring-certificate.json']).read_text())
 assert s['T3_remainder']==s['combined_remainder']==[0]*4 and s['nodes'][151:155]==[[0]*4 for _ in range(4)]
 for p,h in r['input_sha256'].items():pin(ROOT/p,h)
 root=ROOT/'reports/S61_TWO_TAIL_H_TAIL3_ROOT_SOURCE_REVIEW_20261105_v63/report.json'
 pin(root,'99f11d4863e349fc677ac4c037f0001344414b7eeb338603d0db31b7e92efc84')
 for p in (ROOT/'reports/S61_TWO_TAIL_H_TAIL3_AXIS_CUSTODY_20261105_v63').iterdir():
  if p.is_file():pin(p)
 pin(ROOT/'research/notes/S61_TWO_TAIL_H_TAIL3_AXIS_RESULT_20261105_v63.md');pin(Path(__file__))
 rep=dict(id=out.name,all_passed=True,frozen=True,publication_only=True,completed_inconclusive_j3=True,
  saved_primary_four_kappa_remainders_zero=True,no_claim_H_divides_T3_or_kappa_over_char0=True,
  no_other_tail_or_axis_retry=True,next_alternative_minor_SOURCE_only=True,no_bad_Fitting_fullrank_or_exhaustion_claim=True,input_sha256=pins)
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n');print(sha(p))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())

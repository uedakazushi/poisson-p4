#!/usr/bin/env python3
"""Publish fixed D12/Haxis source audit or positive primary; no value replay."""
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
 if kind=='source':
  p=pin(ROOT/'reports/S61_TWO_TAIL_H_D12_AXIS_ROOT_SOURCE_REVIEW_20261105_v63/report.json','c57b649c05838c3ae6df6a61f7f8884ae11aa69d06e38a11414f971ce78db16d')
  r=json.loads(p.read_text());assert r['all_passed'] and r['full_frozen_code_source_custody_and_normalization_read']
  for s,h in r['input_sha256'].items():pin(ROOT/s,h)
  pin(ROOT/'research/notes/S61_TWO_TAIL_H_D12_AXIS_ROOT_SOURCE_REVIEW_20261105_v63.md')
 else:
  base=ROOT/'reports/S61_TWO_TAIL_H_D12_AXIS_QUOTIENT_ENTRY_20261105_v63'
  r=json.loads(pin(base/'report.json','23646500355809a05d19af1774fa33ce75a8ec56f0e31151f52fa56ed76063c8').read_text())
  assert r['all_passed'] and r['computations_completed'] and r['H_nondivides_D12_PRIMARY']
  for s,h in r['input_sha256'].items():pin(ROOT/s,h)
  for s,h in r['evidence_sha256'].items():pin(base/s,h)
  pin(ROOT/'research/notes/S61_TWO_TAIL_H_D12_AXIS_QUOTIENT_RESULT_20261105_v63.md')
 pin(Path(__file__));rep=dict(id=out.name,all_passed=True,frozen=True,publication_only=True,kind=kind,
  no_value_verification_rerun=True,independent_value_peer_NOT_claimed=True,no_actual_bad_Fitting_or_rank_claim=True,input_sha256=pins)
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n');print(sha(p))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--kind',choices=['source','primary'],required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();main(a.kind,a.output.resolve())

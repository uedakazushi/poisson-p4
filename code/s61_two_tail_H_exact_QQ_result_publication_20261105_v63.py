#!/usr/bin/env python3
"""Publish directQQ primary/opaque scope and frozen source audit; no factor arithmetic replay."""
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
 b=ROOT/'reports/S61_TWO_TAIL_H_EXACT_QQ_ENTRY_20261105_v63'
 r=json.loads(pin(b/'report.json','835f6a6404958bf947efa5e1809f49bbe95872ce2033fd89bdab0d574c640a5b').read_text())
 assert r['all_passed'] and r['exact_QQ_factorization_PRIMARY'] and r['primitive_positive_factor_count']==2
 for s,h in r['input_sha256'].items():pin(ROOT/s,h)
 for s,h in r['evidence_sha256'].items():pin(b/s,h)
 pin(ROOT/'reports/S61_TWO_TAIL_H_EXACT_QQ_ROOT_CODE_REVIEW_20261105_v63/report.json','b9be184d3a7bfab6490d1bc1f996794ec51c1dbe0601e319fe7bd86297fc51cd')
 pin(ROOT/'reports/S61_TWO_TAIL_H_EXACT_QQ_PREPARATION_HISTORY_20261105_v63/report.json','781df7e29f8b97806b8b7bc1ab38f0b7425a3887237afb71e1df880c79a4fd2c')
 pin(ROOT/'research/notes/S61_TWO_TAIL_H_EXACT_QQ_RESULT_20261105_v63.md');pin(Path(__file__))
 rep=dict(id=out.name,all_passed=True,frozen=True,publication_only=True,PRIMARY_full_integer_factor_product_only=True,
  independent_backproduct_and_base_irreducibility_peer_NOT_claimed=True,native_internal_growth_remains_opaque=True,
  no_factor_replay_or_Fitting_fullrank_exhaustion_claim=True,input_sha256=pins)
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n');print(sha(p))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())

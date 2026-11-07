#!/usr/bin/env python3
"""Publish completed fixed-image factor result; no root/factor arithmetic replay."""
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
 base=ROOT/'reports/S61_TWO_TAIL_H_ABSOLUTE_SIX_ROOT_ENTRY_20261105_v63'
 r=json.loads(pin(base/'report.json','3028597f392c2c16d3a1aa51065b59cb5702401df3ed2f93f045742dc2ad4fa6').read_text())
 assert r['all_passed'] and r['computations_completed'] and r['fixed_image_factor_found_only'] and r['no_char0_reducibility_claim']
 for s,h in r['input_sha256'].items():pin(ROOT/s,h)
 for s,h in r['evidence_sha256'].items():pin(base/s,h)
 s=json.loads((base/'complete-six-root-certificate.json').read_text());m=s['matched_polynomial_factor_pair']
 assert m['mask']==3 and m['LC_allocation']==0 and len(s['failures'])==2
 assert all(v[1]==0 for k in ['G','Q'] for v in m[k])
 pin(ROOT/'reports/S61_TWO_TAIL_H_ABSOLUTE_SIX_ROOT_ROOT_CODE_REVIEW_20261105_v63/report.json','c68e2870124f99016edcda9d7a4b21b2e2b8946095a53a63480d9a929add117e')
 pin(ROOT/'research/notes/S61_TWO_TAIL_H_ABSOLUTE_SIX_ROOT_RESULT_20261105_v63.md');pin(Path(__file__))
 rep=dict(id=out.name,all_passed=True,frozen=True,publication_only=True,completed_fixed_modular_factor_pair_only=True,
  no_char0_reducibility_or_geometric_irreducibility_claim=True,no_factor_root_or_series_replay=True,
  no_Fitting_fullrank_or_exhaustion_claim=True,input_sha256=pins)
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n');print(sha(p))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())

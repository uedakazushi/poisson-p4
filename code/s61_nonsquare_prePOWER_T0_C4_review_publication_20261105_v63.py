#!/usr/bin/env python3
"""Metadata-only publication of ROOT fixed C4 long-division review; no polynomial arithmetic."""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_T0_C4_ROOT_LONGDIVISION_REVIEW_PUBLICATION_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st):return json.loads(pin(ROOT/'reports'/st/'report.json').read_text())
 peer=load('S61_NONSQUARE_PREPOWER_T0_C4_ROOT_LONGDIVISION_REVIEW_20261105_v63')
 assert peer['all_passed']and peer['independent_actual_quotient_reconstruction']
 primary=load('S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ENTRY_20261105_v63');assert primary['all_passed']and primary['result']['quotient_terms']==4500
 pin(ROOT/'research/notes/S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_RESULT_20261105_v63.md')
 pin(ROOT/'research/notes/S61_NONSQUARE_PREPOWER_T0_C4_ROOT_LONGDIVISION_REVIEW_20261105_v63.md');pin(Path(__file__))
 b=dict(ROOT_peer_sha256=sha(ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C4_ROOT_LONGDIVISION_REVIEW_20261105_v63/report.json'),
  exact_identity='T0=a4*C4*J0/D197; ALL4500 coefficients and full remainder0',
  independence='original7208 QQ literal; fixed25-term C4 leadinga /1004; nofinite-series/Rpowers/10020 engine',
  no_otherfactor_tail_GCD_badlocus_or_global_completion=True)
 (out/'publication-binding.json').write_text(json.dumps(b,indent=2)+'\n')
 r=dict(id=ID,all_passed=True,frozen=True,publication_only=True,input_sha256=pins,
  evidence_sha256={'publication-binding.json':sha(out/'publication-binding.json')})
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

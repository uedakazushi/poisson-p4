#!/usr/bin/env python3
"""Metadata-only Japanese source/value ROOT review publications; no arithmetic rerun."""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(kind,out):
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st):return json.loads(pin(ROOT/'reports'/st/'report.json').read_text())
 stem='S61_NONSQUARE_PREPOWER_Q1_B4_ROOT_'+('SOURCE_REVIEW'if kind=='source'else'MULTINOMIAL_REVIEW')+'_20261105_v63'
 peer=load(stem);assert peer['all_passed'];pin(ROOT/f'research/notes/{stem}.md')
 if kind=='source':
  prepared=load('S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_WORKER_CUSTODY_20261105_v63');assert prepared['all_passed']
  binding=dict(source_code_whole_audit_published=True,no_actual_value_verification_claim_by_this_source_audit=True,failed_comment_guard_preserved=True)
 else:
  actual=load('S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_ENTRY_20261105_v63');assert actual['all_passed']and actual['result']['quotient_terms']==2317
  assert peer['all_original4811_coefficients_full_multiplyback']and peer['independent_closed_multinomial126']
  pin(ROOT/'research/notes/S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_RESULT_20261105_v63.md')
  binding=dict(all4811_original_coefficients_andD197=True,quotient_terms2317=2317,
   independent_method='closed126composition multinomial B4 and171458 ordinaryconvolution products; nolex/nativefactor engine',no_GCD_newfactor_or_complement_claim=True)
 pin(Path(__file__));(out/'publication-binding.json').write_text(json.dumps(binding,indent=2)+'\n')
 r=dict(id=out.name,all_passed=True,frozen=True,publication_only=True,kind=kind,input_sha256=pins,
  evidence_sha256={'publication-binding.json':sha(out/'publication-binding.json')})
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--kind',choices=['source','value'],required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.kind,a.output.resolve())

#!/usr/bin/env python3
"""Publication custody only; no old note/report or mathematical value is changed."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_Q1_F_ROOT_HORNER_REVIEW_PUBLICATION_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 peer=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_F_ROOT_HORNER_REVIEW_20261105_v63/report.json'
 primary=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_F_HASSE_ENTRY_20261105_v63/report.json'
 oldnote=ROOT/'research/notes/S61_NONSQUARE_PREPOWER_Q1_F_HASSE_ENTRY_20261105_v63.md'
 note=ROOT/'research/notes/S61_NONSQUARE_PREPOWER_Q1_F_ROOT_HORNER_REVIEW_20261105_v63.md'
 p=json.loads(peer.read_text());r=json.loads(primary.read_text())
 assert p['all_passed']and p['independent_actual_value_reconstruction']and r['all_passed']and r['result']['F_valuation']==0
 assert sha(peer)=='3d4b5fb36b922aac2714a3fb25327b3c1f884e75e64bd3eba885d409b49bc6b5'
 rel=str(oldnote.relative_to(ROOT))
 if rel in p.get('input_sha256',{}):assert p['input_sha256'][rel]==sha(oldnote)
 paths=[peer,primary,oldnote,note,Path(__file__)]
 d=dict(independent_v_F_QC_zero_ACCEPTED=True,all665_and_common_denominator=True,
  original_primary_note_report_UNCHANGED=True,new_separate_Japanese_REVIEW_note=True,
  no_new_coefficient_or_F_division=True,no_shared_ledger_edit=True,
  source_scope='QQ[a,e,d] fixed actual Q_C; no fullbadlocus/global57/239 claim')
 (out/'publication-binding.json').write_text(json.dumps(d,indent=2)+'\n')
 report=dict(id=ID,all_passed=True,frozen=True,publication_only=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
  evidence_sha256={'publication-binding.json':sha(out/'publication-binding.json')})
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(sha(out/'report.json'))
if __name__=='__main__':main()

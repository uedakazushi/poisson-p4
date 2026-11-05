#!/usr/bin/env python3
"""Immutable T0/C independent-review publication only, no mathematical rerun."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_T0_C_ROOT_JET_HORNER_REVIEW_PUBLICATION_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 peer=ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C_ROOT_JET_HORNER_REVIEW_20261105_v63/report.json'
 primary=ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C_HASSE_ENTRY_20261105_v63/report.json'
 oldnote=ROOT/'research/notes/S61_NONSQUARE_PREPOWER_T0_C_HASSE_ENTRY_20261105_v63.md'
 note=ROOT/'research/notes/S61_NONSQUARE_PREPOWER_T0_C_ROOT_JET_HORNER_REVIEW_20261105_v63.md'
 p=json.loads(peer.read_text());r=json.loads(primary.read_text());assert p['all_passed']and r['all_passed']and r['result']['C_valuation']==4
 assert sha(peer)=='60be3df6763cbfb294727209de1d2fd0bfd672b21488bd0d6753a47a83f2140f'
 rel=str(oldnote.relative_to(ROOT))
 if rel in p.get('input_sha256',{}):assert p['input_sha256'][rel]==sha(oldnote)
 paths=[peer,primary,oldnote,note,Path(__file__)]
 d=dict(independent_vC_T0_exact4_ACCEPTED=True,allfourzero_layers_and646_firstlayer_and_denominator=True,
  separate_Japanese_REVIEW_note=True,old_primary_note_report_UNCHANGED=True,no_new_arithmetic_or_ledger_edit=True,
  no_badlocus_otherfactor_or_quotient_claim=True)
 (out/'publication-binding.json').write_text(json.dumps(d,indent=2)+'\n')
 r=dict(id=ID,all_passed=True,frozen=True,publication_only=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
  evidence_sha256={'publication-binding.json':sha(out/'publication-binding.json')})
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(sha(out/'report.json'))
if __name__=='__main__':main()

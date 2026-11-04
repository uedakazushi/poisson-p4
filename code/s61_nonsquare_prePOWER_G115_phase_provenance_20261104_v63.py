#!/usr/bin/env python3
"""Preserve both phase bytes after a metadata-only final-phase overwrite."""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_G115_PHASE_PROVENANCE_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 out.mkdir(parents=True,exist_ok=False);p=ROOT/'reports/S61_NONSQUARE_PREPOWER_G115_T19_LITERAL_PAIR_PREFLIGHT_20261104_v63';r=json.loads((p/'report.json').read_text());old=dict(phase='source-only same G115 target19; no scalar coefficient work');dump(out/'phase-before-final.json',old);assert sha(out/'phase-before-final.json')==r['evidence_sha256']['phase.json']
 (out/'phase-final-unchanged.json').write_bytes((p/'phase.json').read_bytes());assert json.loads((p/'phase.json').read_text())==dict(phase='FINAL source-only; scalar job NOT launched')
 for name,h in r['evidence_sha256'].items():
  if name!='phase.json':assert sha(p/name)==h
 for rel,h in r['input_sha256'].items():assert sha(ROOT/rel)==h
 note=ROOT/f'research/notes/{ID}.md';report=dict(id=ID,all_passed=True,frozen=True,metadata_only=True,original_preflight_report_sha256=sha(p/'report.json'),original_report_phase_digest_matches_preserved_before_final=True,current_phase_sha256=sha(p/'phase.json'),
  issue='Source preflight driver wrote final phase after computing report.evidence_sha256; only phase.json was replaced. No math, operator, bounds, coefficients, source file or original report changed.',
  frozen_primary_files_untouched=True,all_other_primary_evidence_and_input_hashes_verified=True,preserved_files_sha256={x.name:sha(x)for x in out.iterdir()},code_sha256=sha(Path(__file__)),note_sha256=sha(note))
 dump(out/'report.json',report);print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

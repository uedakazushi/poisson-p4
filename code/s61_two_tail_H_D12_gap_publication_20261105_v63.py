#!/usr/bin/env python3
"""Publish independent Newton acceptance and unconditional gap source, metadata only."""
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
 base=ROOT/'reports/S61_TWO_TAIL_H_D12_ROOT_NEWTON_REVIEW_20261105_v63'
 r=json.loads(pin(base/'report.json','041d9234ccf704286eacd0da6e0f6161d8b268b27dfaa678e2080657f5f32193').read_text())
 assert r['all_passed'] and r['Cayley_Hamilton_matrix_zero'] and r['H_nondivides_D12_in_BASE6_independently_proved']
 for s,h in r['input_sha256'].items():pin(ROOT/s,h)
 for s,h in r['evidence_sha256'].items():pin(base/s,h)
 pin(ROOT/'reports/S61_TWO_TAIL_H_ROOT_SQUAREFREE_20261105_v63/report.json','49990187b3b0b8edaf8939e4d13ac297f4f2ecba7f8b9f5326639089fa043ce2')
 pin(ROOT/'reports/S61_TWO_TAIL_COMMON_FACTOR_SOURCE_20261105_v63/report.json','9a49a8f61d41faf0a2d2d68bfce6817d8c5d35b6d24aa570c33d54d703d844d9')
 for n in ['S61_TWO_TAIL_H_D12_ROOT_NEWTON_REVIEW_20261105_v63.md','S61_TWO_TAIL_H_D12_NONEMPTY_GAP_20261105_v63.md']:pin(ROOT/'research/notes'/n)
 pin(Path(__file__))
 rep=dict(id=out.name,all_passed=True,frozen=True,publication_and_source_only=True,
  accepted_independent_H_nondivides_D12_in_complete_six_factor_base=True,
  nonempty_dimension_two_parameter_open_in_H_intersect_FULLU_D12=True,
  SOME_component_only_not_all_H_components=True,selected_two_tail_sufficient_criterion_gap_only=True,
  no_actual_bad_Fitting_or_rank_drop_or_exhaustion_claim=True,third_tail_values_NOT_executed=True,
  no_arithmetic_rerun=True,input_sha256=pins)
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n');print(sha(p))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())

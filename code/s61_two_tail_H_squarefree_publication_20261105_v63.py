#!/usr/bin/env python3
"""Publish accepted ROOT squarefree proof and codimension-one source consequence only."""
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
 base=ROOT/'reports/S61_TWO_TAIL_H_ROOT_SQUAREFREE_20261105_v63'
 r=json.loads(pin(base/'report.json','49990187b3b0b8edaf8939e4d13ac297f4f2ecba7f8b9f5326639089fa043ce2').read_text())
 assert r['all_passed'] and r['gcd_is_one'] and r['exact_H_squarefree_over_QQ_a_e_d'] and not r['irreducibility_proved']
 for s,h in r['input_sha256'].items():pin(ROOT/s,h)
 for s,h in r['evidence_sha256'].items():pin(base/s,h)
 pin(ROOT/'reports/S61_TWO_TAIL_COMMON_FACTOR_SOURCE_20261105_v63/report.json','9a49a8f61d41faf0a2d2d68bfce6817d8c5d35b6d24aa570c33d54d703d844d9')
 pin(ROOT/'research/notes/S61_TWO_TAIL_H_ROOT_SQUAREFREE_20261105_v63.md');pin(Path(__file__))
 rep=dict(id=out.name,all_passed=True,frozen=True,publication_and_source_only=True,
  accepted_independent_H_squarefree=True,selected_two_tail_codimension_one_part_is_reduced_H=True,
  cofactor_ideal_if_proper_height_at_least_two=True,nonempty_H_divisor_on_six_factor_base_only=True,
  D12_additional_open_and_actual_bad_Fitting_UNRESOLVED=True,no_arithmetic_verification_rerun=True,input_sha256=pins)
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n');print(sha(p))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())

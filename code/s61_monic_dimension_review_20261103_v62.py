#!/usr/bin/env python3
"""Independent exact-input and scope audit of the monic dimension detector."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 directory=ROOT/'reports/S61_AFFINE_MONIC_DIMENSION_CRITERION_20261103_v62'
 path=directory/'report.json';d=json.loads(path.read_text());fresh=json.loads(args.replay.read_text())
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),path]
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 for file,h in d['evidence_sha256'].items():p=directory/file;assert sha(p)==h,file;paths.append(p)
 previous=ROOT/'reports/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62'
 original=(previous/'wronski-exact.sing').read_text().split('ideal I=',1)[1].split(';\n',1)[0].split(',')[:134]
 info=json.loads((directory/'monic-model.json').read_text());order=info['coefficient_coordinate_order']
 assert len(original)==134 and sorted(order)==list(range(45))and order[0]==32
 substitution='cc^45-'+ '-'.join(f'cc^{j}*x{i}'for j,i in enumerate(order)if j)
 assert info['substitution']==substitution
 m2=(directory/'monic-six-exact.m2').read_text();sing=(directory/'monic-six-exact.sing').read_text()
 assert m2.split('I=ideal(',1)[1].split(');\n',1)[0].split(',')==original
 assert sing.split('ideal I=',1)[1].split(';\n',1)[0].split(',')==original
 assert f'x32={substitution};'in m2 and f'poly x32={substitution};'in sing
 with gzip.open(previous/'restricted-original134.json.gz','rt')as stream:rows=json.load(stream)
 counts=Counter()
 for row in rows:
  for mon,c in row:counts.update(set(mon))
 assert min(range(45),key=lambda i:(counts[i],i))==32 and counts[32]==607
 log=(directory/'compatibility-probe.log').read_text();assert 'MONIC_DETECTOR_FINITE_CURVE_AND_NILPOTENT_TRACES_PASS'in log and'error:'not in log
 assert not d['actual_full_unit_result_obtained']and not d['full_affine_dimension_zero_claimed']
 paths += [ROOT/'research/notes/S61_AFFINE_MONIC_DIMENSION_CRITERION_20261103_v62_REVIEW.md',
  ROOT/'reports/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62_REVIEW.json']
 result=dict(id='S61_AFFINE_MONIC_DIMENSION_CRITERION_20261103_v62_REVIEW',all_passed=True,
  fresh_QQ_function_field_M2_toy_unit_nonunit_traces_match=True,all134_original_equations_preserved=True,
  exact_linear_q32_elimination_verified=True,coordinate_permutation_uses_all45_q=True,
  monic_polynomial_unit_iff_dimension_at_most0_lemma_accepted=True,
  proof='Over a geometric integral component, unit gives f*g=h(c). The constant-root prime factors c-a of h distribute between f and g without requiring B to be UFD. Monicity forces every coefficient q_i constant, so no positive-dimensional component survives. Conversely the monic multiplication norm on any finite algebra is nonzero.',
  fixed_pair_finiteness_or_numeric_specialization_assumed=False,full_model_UNIT_obtained=False,
  finite_length_or_five_point_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS exact six-parameter monic unit criterion; no full solver result')

if __name__=='__main__':main()

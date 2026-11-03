#!/usr/bin/env python3
"""Bounded independent QQ span and 13-plane checks; never runs cone GB."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path

import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--prepare-replay',type=Path,required=True)
 parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args()
 directory=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63'
 paths=[Path(__file__),directory/'extraction.json',directory/'coordinate-plane.json',directory/'quadrics.json.gz']
 extraction,plane=[json.loads(p.read_text())for p in paths[1:3]]
 for data in [extraction,plane]:
  for name,expected in data['input_sha256'].items():
   p=ROOT/name;assert sha(p)==expected,name;paths.append(p)
 assert sha(directory/'quadrics.json.gz')==sha(args.prepare_replay/'quadrics.json.gz')
 data=json.loads(gzip.decompress((directory/'quadrics.json.gz').read_bytes()))
 mons=list(itertools.combinations_with_replacement(range(35),2));index={m:i for i,m in enumerate(mons)}
 M=s.zeros(350,len(mons));active=set(plane['active_variable_indices'])
 assert len(active)==13 and [data['variables'][i]for i in plane['active_variable_indices']]==plane['active_variables']
 for row,terms in enumerate(data['all350_quadratic_rows']):
  for i,j,coefficient in terms:
   c=s.Rational(coefficient);assert c and (i,j)in index
   M[row,index[i,j]]+=c
   assert not(i in active and j in active),(row,i,j)
   assert int(c.q)%32003
 assert DomainMatrix.from_Matrix(M).rank()==57
 selected=data['selected_original_rows'];columns=data['pivot_monomial_indices']
 assert len(selected)==len(columns)==57
 minor=M[selected,columns]
 # Independent modular elimination of a selected QQ minor proves p-saturation.
 p=32003;A=[[int(v.p)*pow(int(v.q),-1,p)%p for v in row]for row in minor.tolist()];det=1
 for j in range(57):
  pivot=next(i for i in range(j,57)if A[i][j])
  if pivot!=j:A[j],A[pivot]=A[pivot],A[j];det=-det%p
  value=A[j][j];det=det*value%p;inverse=pow(value,-1,p)
  for i in range(j+1,57):
   factor=A[i][j]*inverse%p
   for k in range(j,57):A[i][k]=(A[i][k]-factor*A[j][k])%p
 assert det
 paths += [directory/'README.md',ROOT/'research/notes/Q_ORDINARY_QUADRATIC_CONE_20261103_v63_REVIEW.md',
           ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63_REVIEW.json']
 report=dict(id='Q_ORDINARY_QUADRATIC_CONE_20261103_V63_REVIEW',all_passed=True,
   fresh_prepare_only_all350_physical_twojet_replay_matches=True,
   independent_QQ_quadratic_row_span_rank=57,nonzero_rows=sum(bool(row)for row in data['all350_quadratic_rows']),
   selected57_by57_minor_mod32003=det,selected_rows_generate_p_integral_QQ_span=True,
   all350_quadrics_vanish_identically_on_listed_QQ_coordinate_plane=True,
   quadratic_obstruction_cone_dimension_lower_bound=13,
   full_actual_model_13_dimensional_branch_claimed=False,
   full_initial_ideal_equals_quadratic_ideal_claimed=False,
   maximum_coordinate_plane_enumeration_repeated=False,cone_GB_repeated=False,
   input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
 print('PASS QQ rank57, p-unit minor, exact13-plane; no higher-branch claim')


if __name__=='__main__':main()

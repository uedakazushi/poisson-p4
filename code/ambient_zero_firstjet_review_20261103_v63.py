#!/usr/bin/env python3
"""Exact physical tensors, recognition matrices, and zero-fiber scope audit."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--replay',type=Path,required=True)
 parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args()
 path=ROOT/'reports/AMBIENT_ZERO_FIRSTJET_FIBERS_20261103_v63/report.json'
 d=json.loads(path.read_text());fresh=json.loads(args.replay.read_text())
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),path]
 for name,expected in d['input_sha256'].items():
  p=ROOT/name;assert sha(p)==expected,name;paths.append(p)
 x,y,z,w=s.symbols('x y z w');Y=[x,y,z,w];pairs=list(itertools.combinations(range(4),2))
 mons=[s.prod(Y[i]for i in ids)for ids in itertools.combinations_with_replacement(range(4),2)]
 J=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
 r=s.Symbol('r');k=s.symbols('k0:9');local={str(v):v for v in [*Y,r,*k]}
 for name,body in d['families'].items():
  entries=[s.sympify(v,locals=local)for v in body['raw_coefficients']];Q=s.zeros(4)
  for a,(i,j)in enumerate(pairs):
   Q[i,j]=sum(entries[10*a+b]*mon for b,mon in enumerate(mons));Q[j,i]=-Q[i,j]
  P=J+Q
  for i,j,h in itertools.combinations(range(4),3):
   actual=sum(P[i,l]*s.diff(P[j,h],Y[l])+P[j,l]*s.diff(P[h,i],Y[l])+P[h,l]*s.diff(P[i,j],Y[l])for l in range(4))
   assert s.expand(actual)==0,(name,i,j,h)
  div=s.Matrix([s.expand(sum(s.diff(Q[i,j],Y[i])for i in range(4)))for j in range(4)])
  assert div==s.Matrix([s.sympify(v,locals=local)for v in body['divergence']])
  assert (div+J*s.Matrix([s.diff(Q[0,1]+Q[2,3],v)for v in Y])).applyfunc(s.expand)==s.zeros(4,1)
  rec=d['full_fixed_line_recognition'][name];C=s.Matrix(rec['full_matrix']);M=s.Matrix(rec['family_embedding'])
  n=9 if name=='R0'else 10
  assert C*M==s.zeros(C.rows,n)and C.rank()==60-n and M.rank()==n
  minor=body['differential_witness'];assert s.Matrix(minor['matrix']).det()==s.Rational(minor['determinant'])!=0
  assert minor['rank']==n+3
 for name,body in d['legacy_crosswalk'].items():assert s.Rational(body['linear_part_determinant'])!=0
 assert [3-s.Matrix(body['matrix']).rank()for body in d['factor_fiber_examples'].values()]==[0,1,3]
 diagnostic=d['rankdrop_diagnostic'];assert s.Matrix(diagnostic['K']).rank()==2
 assert diagnostic['contraction_groebner_basis']==['a']and diagnostic['off_chart_kernel_flag_excluded']
 paths += [ROOT/'research/notes/AMBIENT_ZERO_FIRSTJET_FIBERS_20261103_v63.md',
   ROOT/'research/notes/AMBIENT_ZERO_FIRSTJET_FIBERS_20261103_v63_REVIEW.md',
   ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63_REVIEW.json']
 result=dict(id='AMBIENT_ZERO_FIRSTJET_FIBERS_20261103_V63_REVIEW',all_passed=True,
   fresh_full_QQ_producer_replay_matches=True,independent_physical_J_plus_Q_Jacobi_all_parameters=True,
   full60_coefficient_recognition_matrices_and_embeddings_checked=True,
   exact_source_dimensions=[13,13,12],legacy_M0_M1_R0_identification_accepted=True,
   proper_image_exhausts_underlying_reduced_zero_fibers_accepted=True,
   scheme_source_flag_fibers_point_P1_P3_classification_accepted=True,
   reduced_M0_intersection_M1_equals_R0_accepted=True,
   explicit_irrelevant_saturation_elimination_prime_presentation_accepted=True,
   target_scheme_zero_fiber_reducedness_claimed=False,whole_S61_zero_fiber_exhausted=False,
   input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
 print('PASS whole-family physical tensors and reduced zero-fiber/prime-image scope')


if __name__=='__main__':main()

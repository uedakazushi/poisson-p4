#!/usr/bin/env python3
"""Independent finite-linear certificates for the closed CA/CB incidence."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import sympy as s

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2))
MONS=list(itertools.combinations_with_replacement(range(5),2))
BASIS=list(itertools.product(PAIRS,MONS))
INDEX={v:i for i,v in enumerate(BASIS)}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def physical_representation(g):
 out=s.zeros(150)
 for col,((i,j),mon) in enumerate(BASIS):
  def add(k,l,m,c):
   if k==l:return
   if k>l:k,l,c=l,k,-c
   out[INDEX[(k,l),tuple(sorted(m))],col]+=c
  for k in range(5):
   add(k,j,mon,g[k,i]);add(i,k,mon,g[k,j])
  for removed in range(2):
   derivative=mon[removed]
   for k in range(5):add(i,j,(mon[1-removed],k),-g[derivative,k])
 return out


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--replay',type=Path,required=True)
 parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args()
 source=ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json'
 d=json.loads(source.read_text());fresh=json.loads(args.replay.read_text())
 assert {k:v for k,v in d.items() if k!='seconds'}=={k:v for k,v in fresh.items() if k!='seconds'}
 paths=[Path(__file__),source,
        ROOT/'research/notes/AMBIENT_FULL_FIRSTJET_FIBERS_20261103_v63.md',
        ROOT/'research/notes/AMBIENT_FULL_FIRSTJET_FIBERS_20261103_v63_REVIEW.md',
        ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json',
        ROOT/'reports/KNOWN51_FIRSTJET_IMAGE_DIMENSIONS_20261103_v62_REVIEW.json']
 for name,expected in d['input_sha256'].items():
  p=ROOT/name;assert sha(p)==expected,name;paths.append(p)
 parabolic=[s.Matrix(g) for g in d['parabolic_lie_matrices']]
 assert s.Matrix.hstack(*(g.reshape(16,1) for g in parabolic)).rank()==7
 J=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
 for g in parabolic:
  assert g*J+J*g.T==s.zeros(4)
  assert all(g[i,1]==0 for i in [0,2,3])
 for name,family in d['families'].items():
  M=s.Matrix(family['raw_direction_matrix']);assert M.rank()==16
  for g,action in zip(parabolic,family['parabolic_parameter_actions']):
   full=s.zeros(5);full[1:5,1:5]=g
   assert physical_representation(full)*M==M*s.Matrix(action),(name,'physical Lie action')
  for stem,dimension in [('beta',6),('fiber',10)]:
   matrix=s.Matrix(family['beta_embedding'] if stem=='beta' else family['q_direction_matrix'])
   inverse=s.Matrix(family[stem+'_left_inverse'])
   constraints=s.Matrix(family[stem+'_closed_constraints'])
   assert inverse*matrix==s.eye(dimension)
   assert constraints*matrix==s.zeros(constraints.rows,dimension)
   assert constraints.rank()==matrix.rows-dimension
 A=d['families']['CA'];B=d['families']['CB']
 MA,MB=s.Matrix(A['q_direction_matrix']),s.Matrix(B['q_direction_matrix'])
 assert A['beta_embedding']==B['beta_embedding'] and MA[:,1:]==MB[:,1:]
 assert MA.row_join(MB).rank()==11
 solution=s.Matrix([s.sympify(v) for v in d['intersection_parameter_difference']])
 target=s.Matrix([s.sympify(v) for v in B['q_shift']])-s.Matrix([s.sympify(v) for v in A['q_shift']])
 assert (MA.row_join(-MB[:,:1])*solution-target).applyfunc(s.expand)==s.zeros(45,1)
 assert solution[0]==0 and s.expand(solution[10]-s.Rational(25,4)*s.Symbol('t')**2)==0
 chart=s.Matrix([[s.sympify(v) for v in row]for row in d['line_chart_section']])
 assert chart.T*J*chart==J and s.expand(chart.det())==1
 witness=s.Matrix(d['rank3_kernel_witness']);assert witness.rank()==3
 assert witness.nullspace()==[s.Matrix([0,1,0,0])]
 report=dict(id='AMBIENT_FIBER_BUNDLE_20261103_V63_REVIEW',all_passed=True,
             fresh_full_normalizer_replay_matches=True,
             independent_raw150_physical_Lie_action_all7_by_all16_directions=True,
             all_closed_embedding_left_inverse_and_constraint_ranks_checked=True,
             exact_affine9_intersection_and_common_beta6_checked=True,
             explicit_symplectic_line_chart_section_checked=True,
             proper_closed_incidence_and_rank3_kernel_inverse_accepted=True,
             whole_reduced_CA_CB_rank3_fibers_affine10_accepted=True,
             scheme_intersection_of_reduced_components_affine9_accepted=True,
             proper_smooth_birational_incidence_resolutions_accepted=True,
             unknown_full_S61_fibers_or_rankdrop_fibers_exhausted=False,
             input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
 print('PASS independent physical parabolic actions and closed fiber/intersection certificates')


if __name__=='__main__':main()

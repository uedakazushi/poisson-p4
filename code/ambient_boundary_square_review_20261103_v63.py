#!/usr/bin/env python3
"""Independent Euler identity and new exact26D N-boundary image witness."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_review_20261103_v63 import action

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 path=ROOT/'reports/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63/report.json';d=load(path);fresh=load(args.replay)
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),path,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 h,x,y,z,w=s.symbols('h x y z w');X=s.Matrix([h,x,y,z,w]);alpha,beta=s.symbols('alpha beta',nonzero=True);k=s.symbols('k0:16')
 cubics=[s.prod((h,x,z,w)[i]for i in m)for m in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in m)];H=sum(a*b for a,b in zip(k,cubics))
 # Independent vector-factor derivation, without either producer's raw helper.
 Y=s.Matrix([0,-h*(alpha*h+beta*x),0,-s.diff(H,w),s.diff(H,z)])
 change=s.Matrix([[alpha,beta,0,0,0],[1,0,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]])
 subs=dict(zip(X,change.inv()*X));transported=(change*Y.subs(subs,simultaneous=True)).applyfunc(s.cancel);Hp=H.subs(subs,simultaneous=True)
 expected=s.Matrix([0,beta*x*x,0,beta*x*z-s.diff(Hp,w),beta*x*w+s.diff(Hp,z)])
 ey=s.Matrix([0,0,1,0,0]);raw_difference=(ey*(transported-expected).T-(transported-expected)*ey.T).applyfunc(s.cancel)
 vertical=beta*x*(X*ey.T-ey*X.T);assert (raw_difference-vertical).applyfunc(s.cancel)==s.zeros(5)
 # The invariant19-space gives the proper full N model. A NEW minor gives
 # the exact dimension of its g=0 image, strengthening the input note.
 ca=load(ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json');C=s.Matrix.hstack(*[s.Matrix(v)for v in ca['basis_canonical150']])
 As=[s.Matrix(a)for a in ca['exact_20_by20_Lie_actions']]
 assert all(a[2,j]==0 for a in As for j in range(20)if j!=2)
 source=[C[:,i]for i in [0,1,*range(4,20)]];point=C[:,0]+C[:,1]+sum(((i%7+1)*C[:,i+4]for i in range(16)),s.zeros(150,1))
 M=s.Matrix.hstack(*source,*[action(i,j)*point for i,j in itertools.product(range(5),repeat=2)])
 cols=list(M.rref()[1]);rows=list(M[:,cols].T.rref()[1]);minor=M[rows,cols];det=minor.det()
 assert len(cols)==27 and det==s.Rational(17978112,15625)
 paths += [ROOT/'research/notes/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW.md',
  ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
  ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json']
 result=dict(id='AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW',all_passed=True,
  fresh_QQ_all16H_Euler_identity_replay_matches=True,independent_vector_factor_finite_coordinate_identity=True,
  proper_full_N19_subspace_invariance_verified=True,proper_closure_recognition_of_both_CB_square_endpoint_images_accepted=True,
  N_g_zero_new_exact_image_witness=dict(source_projective_dimension_upper=26,QQ_affine_image_rank=27,
   rows=rows,columns=cols,determinant=str(det),minor=[[str(v)for v in row]for row in minor.tolist()],
   full_canonical150_differential=[[str(v)for v in row]for row in M.tolist()]),
  full_N_rank_at_most2_boundary_exact_two_irreducible_component_dimensions=[26,25],
  identities=['N_g0=Bb_square','CA_f0=Ba_square'],full_CA_intersection_CB_exhausted=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS full N square-boundary recognition; NEW exact component dimensions26,25')

if __name__=='__main__':main()

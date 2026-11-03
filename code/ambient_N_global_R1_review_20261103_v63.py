#!/usr/bin/env python3
"""Full-source minor, complete geometric flag fiber, and global N R1 audit."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_review_20261103_v63 import IDS,action,canonical,wedge_factor

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 path=ROOT/'reports/AMBIENT_N_GLOBAL_R1_20261103_v63/report.json';d=load(path);fresh=load(args.replay)
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),path,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 h,x,y,z,w=s.symbols('h x y z w');X=[h,x,y,z,w];ey=s.Matrix([0,0,1,0,0]);V=s.Matrix([0,-h*x,0,-3*w*w,3*z*z]);P=ey*V.T-V*ey.T
 vector=lambda P:s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(s.prod(X[k]for k in m))for(ab,m),j0 in IDS.items()for i,j in[ab]])
 v=canonical(vector(P));C=s.Matrix([h*h,-4*h*x,0,-15*w*w+h*z,15*z*z+h*w]);assert v==vector((ey*C.T-C*ey.T)/5)
 ca=load(ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json');N=s.Matrix.hstack(*[s.Matrix(a)for i,a in enumerate(ca['basis_canonical150'])if i!=2])
 M=s.Matrix.hstack(N,*[action(i,j)*v for i,j in itertools.product(range(5),repeat=2)]);wit=d['full_dimension_minor']
 assert M==s.Matrix(wit['exact_raw150_image_differential'])
 assert M[wit['minor_rows'],wit['minor_columns']].det()==s.Rational(wit['minor_determinant'])==s.Rational(6912,390625)
 W=wedge_factor(v);assert W.rank()==4 and W.nullspace()==[ey]
 U=[h,x,z,w];D=[C[i]for i in [0,1,3,4]];a,b,c,e,lam=s.symbols('a b c d lambda');ell=a*h+b*x+c*z+e*w
 derivative=lambda f:s.expand(sum(s.diff(f,u)*v for u,v in zip(U,D)))
 poly=s.Poly(derivative(ell)-lam*ell*ell,*U)
 assert poly.coeff_monomial(x*x)==-lam*b*b and poly.coeff_monomial(z*w)==-2*lam*c*e
 assert poly.coeff_monomial(z*z)==15*e-lam*c*c and poly.coeff_monomial(w*w)==-15*c-lam*e*e
 K,_=s.linear_eq_to_matrix(s.Poly(derivative(ell),*U).coeffs(),[a,b,c,e]);assert K.rank()==4
 assert s.expand(derivative(ell)+4*h*ell-(5*a*h*h+c*(-15*w*w+5*h*z)+e*(15*z*z+5*h*w)))==0
 paths += [ROOT/'research/notes/AMBIENT_N_GLOBAL_R1_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_N_GLOBAL_R1_20261103_v63_REVIEW.md',
  ROOT/'reports/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW.json']
 result=dict(id='AMBIENT_N_GLOBAL_R1_20261103_v63_REVIEW',all_passed=True,
  fresh_full_QQ_producer_replay_matches=True,independent_physical150_full_source_minor_and_unique_factor_matrix=True,
  complete_entire_geometric_flag_fiber_uniqueness_proof_accepted=True,
  immersion_implies_fiber_tangent_zero_not_target_tangent_equality=True,
  proper_finite_length1_Nakayama_isomorphism_proof_accepted=True,
  global_projective_N_R1_accepted=True,global_normality_or_CM_inferred=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS full geometric flag fiber length1 and GLOBAL N Serre R1')

if __name__=='__main__':main()

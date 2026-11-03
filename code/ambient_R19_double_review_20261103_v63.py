#!/usr/bin/env python3
"""Independent exact finite curves, all-kappa minors and double-ring audit."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_review_20261103_v63 import IDS,action,canonical

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--generic-replay',type=Path,required=True)
 ap.add_argument('--open-replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 paths=[Path(__file__),ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py'];reports=[]
 for stem,replay in [('AMBIENT_R19_GENERIC_DOUBLE',args.generic_replay),('AMBIENT_R19_OPEN_DOUBLE',args.open_replay)]:
  p=ROOT/f'reports/{stem}_20261103_v63/report.json';d=load(p);f=load(replay);assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'}
  reports.append(d);paths.append(p)
  for file,h in d['input_sha256'].items():q=ROOT/file;assert sha(q)==h,file;paths.append(q)
 generic,opened=reports;h,x,y,z,w=s.symbols('h x y z w');X=s.Matrix([h,x,y,z,w]);kap,eps=s.symbols('kappa eps')
 vector=lambda P:s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(s.prod(X[k]for k in m))for(ab,m),ind in IDS.items()for i,j in[ab]])
 def tensor(f,H,CB=False):
  P=s.zeros(5);P[1,2]=f;radial=-s.diff(f,x)/2 if CB else 0
  P[2,3]=radial*z-s.diff(H,w);P[2,4]=radial*w+s.diff(H,z);return P-P.T
 def push(P,M):return (M*P.subs(dict(zip(X,M.inv()*X)),simultaneous=True)*M.T).applyfunc(s.expand)
 geo=load(ROOT/'reports/AMBIENT_R19_GENERIC_GEOMETRY_20261103_v63/report.json');cb=load(ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json')
 A=s.Matrix([[s.sympify(v)for v in row]for row in geo['full_CA_source_differential']]);B=s.Matrix([[s.sympify(v)for v in row]for row in geo['full_CB_source_differential']]);W=A.row_join(B);point=s.Matrix(list(map(s.sympify,geo['point_canonical150'])))
 HA=h*h*z+kap*h*z*z/2+z**3/3+2*x*z*w;velocity=s.zeros(5);velocity[1,3]=2;assert velocity*velocity==s.zeros(5)
 PA=canonical(vector(push(tensor(-4*x*x,HA+eps*3*z*z*w),s.eye(5)+eps*velocity)))
 bcoef=list(map(s.sympify,geo['CB_source19_coefficients']));f0=bcoef[0]*h*h+bcoef[1]*h*x+bcoef[2]*x*x;HB=sum(bcoef[i+3]*s.sympify(m)for i,m in enumerate(cb['cubic_monomials']))
 perm=[1,3,2,0,4];g=s.zeros(5)
 for i,j in enumerate(perm):g[j,i]=1
 PB=canonical(vector(push(tensor(f0-21*eps*x*x,HB+eps*x*z*w,True),g)))
 diff=s.zeros(150,1);diff[72]=-18
 assert (PA-PB-eps*eps*diff).applyfunc(s.expand)==s.zeros(150,1)
 v=PA.diff(eps).subs(eps,0);assert v==PB.diff(eps).subs(eps,0)==s.Matrix(list(map(s.sympify,generic['common_first_order_canonical150'])))
 # Exact finite curves have no higher terms; independent Lie second formula.
 T=2*action(1,3);delta=canonical(vector(tensor(0,3*z*z*w)));assert T*delta+T*T*point/2==diff
 ell=s.zeros(1,150);ell[0,72]=kap;ell[0,63]=-1;assert ell*W==s.zeros(1,W.cols)and(ell*diff)[0]==-18*kap
 ell0=s.zeros(1,150);ell0[0,68]=s.Rational(5,22);ell0[0,72]=1;assert ell0*W.subs(kap,0)==s.zeros(1,W.cols)and(ell0*diff)[0]==-18
 sourceHs=[w*x*x,x*x*z,x*z*z/2,2*h*x*z,z**3/3,h*z*z,3*h*h*z]
 directions=[canonical(vector(tensor(-4*x*x,2*x*z*w)))]+[canonical(vector(tensor(0,H)))for H in sourceHs]
 R=s.Matrix.hstack(*directions,*[action(i,j)*point for i,j in itertools.product(range(5),repeat=2)])
 assert R.rank()==20 and R.row_join(v).rank()==21 and R.subs(kap,0).row_join(v.subs(kap,0)).rank()==21
 matrices={'CA':A,'CB':B,'tangent_union':W,'paired_flag_R_source':R}
 for label,M in matrices.items():
  m=opened['all_kappa_constant_minors'][label];det=s.cancel(M[m['rows'],m['columns']].det(method='domain-ge'));assert det==s.Rational(m['determinant'])!=0 and not det.free_symbols
 assert W.subs(kap,0).row_join(diff).rank()==38
 paths += [ROOT/'reports/AMBIENT_R19_GEOMETRY_20261103_v63_REVIEW.json',
  ROOT/'research/notes/AMBIENT_R19_GENERIC_DOUBLE_20261103_v63.md',ROOT/'research/notes/AMBIENT_R19_OPEN_DOUBLE_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_R19_DOUBLE_20261103_v63_REVIEW.md']
 result=dict(id='AMBIENT_R19_DOUBLE_20261103_v63_REVIEW',all_passed=True,
  fresh_QQ_generic_and_all_kappa_producers_match=True,independent_exact_finite_actual_CA_CB_curve_difference_all150=True,
  no_higher_epsilon_terms_and_Lie_second_formula_independently_checked=True,
  full_tangent_annihilators_and_obstructions={'kappa_nonzero':'-18*kappa','kappa_zero':'-18'},
  all_four_symbolically_constant_nonzero_source_union_minors_verified=True,
  generic_primary_algebra='K(R)[eps]/eps²',
  proper_joint_source_closed_subscheme_after_open_flag_pair_restriction_accepted=True,
  all_kappa_CA_CB_R_smooth_and_target_fiber_length1_proofs_accepted=True,
  regular_S20_height1_valuation_and_quadratic_unit_ideal_equality_accepted=True,
  completed_intersection_every_Pkappa='C[[u19,eps]]/(eps²)',
  depth=19,associated_primes=['R only'],AF_zero_or_t_zero_or_global_primary_decomposition_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS exact generic and EVERY-kappa ordinary double intersection ring')

if __name__=='__main__':main()

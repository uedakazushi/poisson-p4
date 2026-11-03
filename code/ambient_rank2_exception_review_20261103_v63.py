#!/usr/bin/env python3
"""Independent direct CB membership and intrinsic global-N exclusion witness."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_review_20261103_v63 import IDS,canonical,wedge_factor

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 path=ROOT/'reports/AMBIENT_RANK2_EXCEPTION_NOT_N_20261103_v63/report.json';d=load(path);fresh=load(args.replay)
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),path,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 h,x,y,z,w=s.symbols('h x y z w');X=[h,x,y,z,w];U=[h,x,z,w]
 V=s.Matrix([0,4*x*x,0,-2*x*z,3*h*h+z*z+2*x*w]);ey=s.Matrix([0,0,1,0,0]);P=ey*V.T-V*ey.T
 vector=lambda P:s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(s.prod(X[k]for k in m))for(ab,m),index in IDS.items()for i,j in[ab]])
 raw=vector(P);cv=canonical(raw)
 C=s.Matrix([-s.Rational(8,5)*h*x,s.Rational(12,5)*x*x,-s.Rational(18,5)*x*z,3*h*h+z*z+s.Rational(2,5)*x*w])
 field=s.Matrix([C[0],C[1],0,C[2],C[3]]);Q=ey*field.T-field*ey.T;assert vector(Q)==cv
 # Literal direct membership in the TWO actual canonical family images.
 ca=load(ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json');caC=s.Matrix.hstack(*[s.Matrix(v)for v in ca['basis_canonical150']])
 ca_rows=list(caC.T.rref()[1]);ca_coeff=caC[ca_rows,:].inv()*cv[ca_rows,0];assert caC*ca_coeff==cv and ca_coeff[3]==0
 cb=load(ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json');cbC=s.Matrix.hstack(*[s.Matrix(v)for v in cb['enlarged_canonical_basis']])
 endpoint=cbC[:,[*range(3),*range(6,22)]];coeff=s.Matrix([0,6,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0]);before=endpoint*coeff
 perm=[1,3,2,0,4];after=s.zeros(150,1)
 for(ab,m),index in IDS.items():
  ij=[perm[i]for i in ab];after[IDS[(tuple(sorted(ij)),tuple(sorted(perm[i]for i in m)))]]+=(1 if ij[0]<ij[1]else-1)*before[index]
 assert after==cv
 W=wedge_factor(cv);assert W.rank()==4 and W.nullspace()==[ey]
 ell5=s.Matrix(s.symbols('ell0:5'));contract=Q*ell5
 mono=[s.prod(X[i]for i in m)for m in itertools.combinations_with_replacement(range(5),2)]
 equations=[s.Poly(v,*X).coeff_monomial(m)for v in contract for m in mono];K,_=s.linear_eq_to_matrix(equations,list(ell5));assert K.rank()==5
 a,b,c,e,lam,alpha=s.symbols('a b c d lambda alpha');ell=a*h+b*x+c*z+e*w
 derivative=lambda f:s.expand(sum(s.diff(f,u)*v for u,v in zip(U,C)))
 first=s.Poly(derivative(ell)-lam*ell**2,*U)
 assert first.coeff_monomial(w*w)==-lam*e*e and first.coeff_monomial(h*h)==3*e-lam*a*a and first.coeff_monomial(z*z)==e-lam*c*c
 matrix,_=s.linear_eq_to_matrix(s.Poly(derivative(ell),*U).coeffs(),[a,b,c,e]);assert matrix.rank()==4
 second=derivative(ell)+s.Rational(48,5)*x*ell+alpha*x*x
 S,_=s.linear_eq_to_matrix(s.Poly(second,*U).coeffs(),[a,b,c,e,alpha]);assert S.nullspace()==[s.Matrix([0,-s.Rational(1,12),0,0,1])]
 paths += [ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
  ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
  ROOT/'reports/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW.json',
  ROOT/'research/notes/AMBIENT_RANK2_EXCEPTION_NOT_N_20261103_v63_REVIEW.md']
 result=dict(id='AMBIENT_RANK2_EXCEPTION_NOT_N_20261103_v63_REVIEW',all_passed=True,
  fresh_QQ_intrinsic_witness_producer_replay_matches=True,raw_physical_vector_field=list(map(str,V)),
  exact_CA_canonical_basis_coefficients=list(map(str,ca_coeff)),
  independent_CB_b_zero_membership=dict(permutation_old_to_new=perm,canonical_basis_coefficients=list(map(str,coeff)),all150_entries_match=True),
  unique_intrinsic_constant_factor_rank=4,intrinsic_factor_kernel=[0,0,1,0,0],constant_covector_contraction_rank=5,
  necessary_global_N_g0_two_covector_conditions_impossible=True,witness_outside_global_N=True,
  actual_global_CA_intersection_CB_point_outside_N_certified=True,
  extra_family_dimension_or_global_exhaustion_from_this_packet_alone_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS actual rank2 CA/CB common point outside global N, direct150 certificate')

if __name__=='__main__':main()

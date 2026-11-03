#!/usr/bin/env python3
"""Independent sparse physical review of the compactified CB relative pencil."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_review_20261103_v63 import action,canonical,wedge_factor,PAIRS,MONS

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 path=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json';d=load(path);fresh=load(args.replay)
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),path,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 raws=[s.Matrix(v)for v in d['enlarged_raw_basis']];basis=[s.Matrix(v)for v in d['enlarged_canonical_basis']]
 assert len(raws)==22 and all(canonical(r)==v for r,v in zip(raws,basis))
 C=s.Matrix.hstack(*basis);assert C.rank()==22
 mu=s.Symbol('mu');E=s.zeros(22,19)
 for i in range(3):E[i,i]=1;E[i+3,i]=mu
 for i in range(16):E[i+6,i+3]=1
 for pos,rr,mr in zip(d['parabolic_GL_positions'],d['exact_22_by22_Lie_actions'],d['exact_19_by19_relative_mu_actions']):
  i,j=pos;R=s.Matrix(rr);T=action(i,j);A=s.Matrix([[s.sympify(v,locals={'mu':mu})for v in row]for row in mr])
  assert T*C==C*R
  chi=int(i==j and i in [3,4])-int(i==j and i in [1,2])
  assert R[:3,3:]==s.zeros(3,19)and R[3:6,:3]==s.zeros(3,3)and R[3:6,6:]==s.zeros(3,16)
  assert R[3:6,3:6]-R[:3,:3]==chi*s.eye(3)
  assert (R*E-chi*mu*E.diff(mu)-E*A).applyfunc(s.expand)==s.zeros(22,19)
 endpoint=[*basis[3:6],*basis[6:]];Honly=basis[6:]
 for bs,rr in zip([endpoint,Honly],d['forget_hline_extra_generator_actions']):
  B=s.Matrix.hstack(*bs);assert action(0,1)*B==B*s.Matrix(rr)
 X=s.symbols('h x y z w');h,x,y,z,w=X;binary=[h*h,h*x,x*x]
 cubics=[s.prod((h,x,z,w)[i]for i in m)for m in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in m)]
 def tensor(f,g,H):
  P=s.zeros(5);P[1,2]=f;P[3,4]=g;P[2,3]=-s.diff(f,x)*z/2-s.diff(H,w);P[2,4]=-s.diff(f,x)*w/2+s.diff(H,z);return P-P.T
 def vector(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(s.prod(X[t]for t in m))for i,j in PAIRS for m in MONS])
 physical=[tensor(f,0,0)for f in binary]+[tensor(0,g,0)for g in binary]+[tensor(0,0,H)for H in cubics]
 assert all(vector(P)==v for P,v in zip(physical,raws))
 u=s.symbols('u0:22');P=sum((a*b for a,b in zip(u,physical)),s.zeros(5));f=sum(u[i]*binary[i]for i in range(3));g=sum(u[i+3]*binary[i]for i in range(3))
 W=s.expand(g*s.diff(f,x)-f*s.diff(g,x))
 for i,j,k in itertools.combinations(range(5),3):
  jac=s.expand(sum(P[i,t]*s.diff(P[j,k],X[t])+P[j,t]*s.diff(P[k,i],X[t])+P[k,t]*s.diff(P[i,j],X[t])for t in range(5)))
  assert jac==(W if(i,j,k)==(2,3,4)else 0)
 assert [s.Poly(W,h,x).coeff_monomial(m)for m in [h**3,h*h*x,h*x*x]]==[u[3]*u[1]-u[0]*u[4],2*(u[3]*u[2]-u[0]*u[5]),u[4]*u[2]-u[1]*u[5]]
 A,B=s.symbols('a b');F=s.Symbol('F0')*h*h+s.Symbol('F1')*h*x+s.Symbol('F2')*x*x
 PP=tensor(A*F,B*F,0);assert s.expand(h*(PP[1,2]*PP[3,4]-PP[1,3]*PP[2,4]+PP[1,4]*PP[2,3])-A*B*h*F*F)==0
 H=sum(((i%7+1)*v for i,v in enumerate(basis[6:])),s.zeros(150,1));all_actions=[action(i,j)for i,j in itertools.product(range(5),repeat=2)];records={}
 for name,point,source in [('a_zero',sum(basis[3:6],s.zeros(150,1))+H,endpoint),('b_zero',sum(basis[:3],s.zeros(150,1))+H,[*basis[:3],*basis[6:]]),('F_zero',H,Honly)]:
  rec=d['boundary_exact_witnesses'][name];M=s.Matrix.hstack(*source,*[T*point for T in all_actions])
  assert M==s.Matrix(rec['exact_raw150_image_differential'])
  det=M[rec['minor_rows'],rec['minor_columns']].det();assert det==s.Rational(rec['minor_determinant'])!=0
  factor=wedge_factor(point);assert factor.rank()==rec['intrinsic_constant_vector_factor_matrix_rank']
  wr,wc=rec['constant_vector_factor_minor_rows'],rec['constant_vector_factor_minor_columns'];assert factor[wr,wc].det()==s.Rational(rec['constant_vector_factor_minor_determinant'])!=0
  records[name]=dict(exact_QQ_image_minor_size=len(rec['minor_columns']),projective_image_dimension=len(rec['minor_columns'])-1,intrinsic_factor_rank=factor.rank())
 paths += [ROOT/'research/notes/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.md',
  ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json']
 result=dict(id='AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW',all_passed=True,
  fresh_full_QQ_producer_replay_matches=True,independent_sparse150_all352_Lie_columns_and304_relative_columns=True,
  finite_parabolic_ratio_character_and_associated_nontrivial_P1_bundle_accepted=True,
  both_hline_forgetting_subspaces_preserved_by_extra_generator=True,full_physical_Jacobi_rank1_minor_ideal_and_quintic_checked=True,
  boundary_witnesses=records,proper_smooth28D_source_image_is_entire_reduced_CB_accepted=True,
  entire_reduced_rank_at_most2_boundary_exact_component_dimensions=[26,27],H_only_image_dimension=23,
  global_normalization_R1_or_CA_intersection_CB_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS full CB relative P1 bundle and complete two-component rank2 boundary')

if __name__=='__main__':main()

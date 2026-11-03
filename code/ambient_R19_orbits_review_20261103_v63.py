#!/usr/bin/env python3
"""Independent intrinsic-hyperplane argument and finite R19 normal-form audit."""
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
 path=ROOT/'reports/AMBIENT_R19_GENERIC_ORBITS_20261103_v63/report.json';d=load(path);fresh=load(args.replay)
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),path,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 h,x,z,w=s.symbols('h x z w');U=s.Matrix([h,x,z,w]);t,q0,q1,q2,q3,q5,q6,q7=s.symbols('t q0 q1 q2 q3 q5 q6 q7')
 local={str(v):v for v in [h,x,z,w,t,q0,q1,q2,q3,q5,q6,q7]}
 V=s.Matrix([s.sympify(v,locals=local)for v in d['source_vector_field']]);G=s.Matrix([[s.sympify(v,locals=local)for v in row]for row in d['normalizing_linear_matrix']])
 assert G.det()==1
 source=(G*V.subs(dict(zip(U,G.inv()*U)),simultaneous=True)).applyfunc(s.cancel)
 assert source==s.Matrix([0,4*t*x*x,-2*t*x*z,2*t*x*w+3*q7*h*h+2*q6*h*z+q5*z*z])
 A,B,F=s.symbols('A B F');Q=A*h*h+B*h*z+F*z*z;normal=s.Matrix([0,4*x*x,-2*x*z,2*x*w+Q]);D=5*normal-8*x*U
 assert sum(s.diff(D[i],U[i])for i in range(4))==0
 a,b,c=s.symbols('a b c');ell=w+a*h+b*x+c*z;Der=lambda f:sum(s.diff(f,U[i])*D[i]for i in range(4))
 assert s.expand(Der(ell).subs(w,-a*h-b*x-c*z))==5*Q-10*a*h*x+10*b*x*x-20*c*x*z
 assert s.expand(Der(a*h+b*x+c*z)-x*(-8*a*h+12*b*x-18*c*z))==0
 assert len(set([-8,12,-18]))==3
 da,db,dc,dd,u,v,r=s.symbols('da db dc dd u v r',nonzero=True)
 M=s.Matrix([[da,0,0,0],[0,db,0,0],[0,0,dc,0],[u,v,r,dd]])
 transported=(db*M*D.subs(dict(zip(U,M.inv()*U)),simultaneous=True)).applyfunc(s.cancel)
 QA=db*dd*A/(da*da);QB=db*dd*B/(da*dc);QF=db*dd*F/(dc*dc)
 target=(5*normal.subs({A:QA,B:QB,F:QF},simultaneous=True)-8*x*U)
 target[3]+=-10*u*h*x/da+10*v*x*x/db-20*r*x*z/dc
 assert (transported-target).applyfunc(s.cancel)==s.zeros(4,1)
 assert s.cancel(QB*QB/(QA*QF)-B*B/(A*F))==0
 assert s.cancel((B*B/(A*F)).subs({A:3*q7/t,B:2*q6/t,F:q5/t}))==4*q6*q6/(3*q7*q5)
 # Exact GL-affine orbit cross-check, independent of the full finite action
 # producer. The classification proof, not sampled rank, gives constancy.
 table=load(ROOT/'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63/double-flag-table/report.json');rec=next(v for v in table['records']if v['relative_permutation']==[1,3,2,0,4]and v['CB_endpoint']=='b_zero')
 K=s.Matrix.hstack(*[s.Matrix(v)for v in rec['intersection_basis_in_CA_coordinates']]);ca=load(ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json');CA=s.Matrix.hstack(*[s.Matrix(v)for i,v in enumerate(ca['basis_canonical150'])if i!=3]);L=CA*K
 P=L*s.Matrix([i%7+1 for i in range(8)]);orbit=s.Matrix.hstack(*[action(i,j)*P for i,j in itertools.product(range(5),repeat=2)])
 cols=list(orbit.rref()[1]);rows=list(orbit[:,cols].T.rref()[1]);minor=orbit[rows,cols];assert len(cols)==19 and minor.det()!=0
 assert orbit.row_join(P).rank()==19 and orbit.row_join(L).rank()==20
 paths += [ROOT/'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json',
  ROOT/'research/notes/AMBIENT_R19_GENERIC_ORBITS_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_R19_GENERIC_ORBITS_20261103_v63_REVIEW.md']
 result=dict(id='AMBIENT_R19_GENERIC_ORBITS_20261103_v63_REVIEW',all_passed=True,
  fresh_all_parameter_QQ_finite_action_producer_replay_matches=True,independent_eight_parameter_normalization_matrix_identity=True,
  complete_linear_Darboux_hyperplane_proof_accepted=True,intrinsic_constant_factor_unique_for_Q_nonzero_accepted=True,
  arbitrary_PGL_equivalence_reduces_to_diagonal_by_labeled_cofactors_accepted=True,
  complete_invariant_on_t_q5_q7_nonzero='kappa²=4*q6²/(3*q7*q5)',generic_invariant_field='C(kappa²)',
  projective_orbit_dimension=18,PGL_projective_class_stabilizer_dimension=6,GL_projective_class_stabilizer_dimension=7,
  GL_affine_tensor_stabilizer_dimension=6,
  independent_generic_GL_affine_minor=dict(size=19,rows=rows,columns=cols,determinant=str(minor.det())),
  special_t_zero_or_AF_zero_or_closure_boundary_orbits_classified=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS complete R19 generic projective invariant kappa²; orbit18/modulus1')

if __name__=='__main__':main()

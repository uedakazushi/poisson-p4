#!/usr/bin/env python3
"""Independent dual-number powers, full150 differentials and singularity audit."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from ambient_CA_global_flag_review_20261103_v63 import IDS,action

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def sparse(d):
 M=s.zeros(*d['shape'])
 for i,j,v in d['nonzero']:M[i,j]=s.sympify(v)
 return M
def reduce(M):return M.applyfunc(s.cancel)
def wedge(A,B):
 return s.Matrix([A[i,p]*B[j,q]-B[i,p]*A[j,q]+(A[i,q]*B[j,p]-B[i,q]*A[j,p]if p!=q else 0)for((i,j),(p,q))in IDS])

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_T01_SOURCE_DIFFERENTIAL_20261103_v63';public=folder/'report.json';d=load(public);f=load(args.replay/'report.json')
 assert{k:v for k,v in d.items()if k not in('seconds','matrices_sha256')}=={k:v for k,v in f.items()if k not in('seconds','matrices_sha256')}
 packed=folder/'matrices.json.gz';assert sha(packed)==d['matrices_sha256']
 raw=gzip.decompress(packed.read_bytes());assert raw==gzip.decompress((args.replay/'matrices.json.gz').read_bytes());payload=json.loads(raw)
 paths=[Path(__file__),public,packed,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 u=s.Symbol('s');K=s.QQ.frac_field(u);D=lambda M:DomainMatrix.from_Matrix(M).convert_to(K).to_dense()
 N=s.zeros(5)
 for i,v in enumerate([1,u,(u+3)/2,-9/(u+3)]):N[i+1,i]=v
 powers=[s.eye(5)]
 for i in range(1,8):powers.append(reduce(powers[-1]*N))
 a=(4*u+6)/(3*u*(u+3));lam=-2*a/(9*u);mu=-4/(9*u*u*(u+3))
 rawtensor=reduce(lam*wedge(powers[1],powers[4])+mu*wedge(powers[2],powers[3]));q=s.zeros(150,1)
 for key,v in{((1,4),(0,0)):a,((2,3),(0,0)):-s.Rational(2,9),((2,4),(0,1)):s.Rational(4,3),((3,4),(0,2)):s.Rational(2,3),((3,4),(1,1)):1}.items():q[IDS[key]]=v
 assert reduce(rawtensor-q)==s.zeros(150,1)
 columns=[]
 for i,j in itertools.product(range(5),repeat=2):
  E=s.zeros(5);E[i,j]=1;value=s.eye(5);tangent=s.zeros(5);derivatives=[]
  # Multiplication in dual numbers (N+eps E)^k, eps²=0.  This is
  # independent of the producer's closed Frechet summation formula.
  for k in range(1,5):
   tangent=reduce(tangent*N+value*E);value=reduce(value*N)
   derivatives.append(tangent-s.trace(tangent)*s.eye(5)/5)
  columns.append(reduce(lam*(wedge(derivatives[0],powers[4])+wedge(powers[1],derivatives[3]))+mu*(wedge(derivatives[1],powers[3])+wedge(powers[2],derivatives[2]))))
 powerpairs=list(itertools.combinations(range(1,5),2));columns.extend(reduce(wedge(powers[i],powers[j]))for i,j in powerpairs)
 J=s.Matrix.hstack(*columns);assert J==sparse(payload['jacobian'])
 sigma=-3*(u+3)/(2*u+3);other=reduce(J.subs(u,sigma,simultaneous=True));assert other==sparse(payload['other_jacobian'])
 assert s.cancel(a.subs(u,sigma)-a)==0
 O=s.Matrix.hstack(*[action(i,j)*q for i,j in itertools.product(range(5),repeat=2)]);assert D(O)==D(sparse(payload['strict_GL5_action']))
 def comm(M):
  cols=[]
  for i,j in itertools.product(range(5),repeat=2):
   E=s.zeros(5);E[i,j]=1;v=reduce(E*M-M*E);cols.append(s.Matrix(list(v)+[0]*6))
  return s.Matrix.hstack(*cols)
 C=comm(N);assert D(J)*D(C)==D(O)and D(other)*D(comm(reduce(N.subs(u,sigma,simultaneous=True))))==D(O)
 def rank(M):return len(D(M).rref()[1])
 assert rank(J)==24 and rank(other)==24 and rank(O)==17
 joint=J.row_join(other);assert rank(joint)==30
 for label,M in [('source_jacobian',J),('strict_GL5_action',O),('conjugate_source_joint',joint)]:
  rec=d[label];det=D(M[rec['pivot_rows'],rec['pivot_columns']]).det();assert s.cancel(K.to_sympy(det)-s.sympify(rec['nonzero_minor']))==0
 stab=sparse(payload['strict_stabilizer_basis']);assert stab.shape==(25,8)and rank(stab)==8 and D(O)*D(stab)==D(s.zeros(150,8))
 stabsource=C*stab;assert D(stabsource)==D(sparse(payload['stabilizer_conjugation_kernel']))and rank(stabsource)==3
 coeff=J[:,25:];rows=list(D(coeff.T).rref()[1]);assert len(rows)==6
 inverse=coeff[rows,:].inv();poly=[]
 for r in range(5):
  direction=s.Matrix(list(powers[r]));rhs=J[:,:25]*direction;dc=reduce(-inverse*rhs[rows,:]);vector=direction.col_join(dc)
  assert D(J)*D(vector)==D(s.zeros(150,1));poly.append(vector)
 poly=s.Matrix.hstack(*poly);assert D(poly)==D(sparse(payload['polynomial_generator_kernel']))and rank(poly)==5
 spanning=poly.row_join(stabsource);assert D(spanning)==D(sparse(payload['source_fiber_spanning_directions']))and rank(spanning)==7
 center=s.Matrix.hstack(*[s.Matrix(list(powers[r]))for r in range(5)]);assert rank(center)==5 and D(O)*D(center)==D(s.zeros(150,5))
 G=s.zeros(5);G[4,1]=-9*u/2;assert reduce(G*N-N*G-powers[4])==s.zeros(5)
 qdot=s.zeros(150,1);qdot[IDS[((1,4),(0,0))]]=1;common=O.row_join(qdot)
 assert D(common)==D(sparse(payload['common_GL_orbit_and_Qa_curve']))and rank(common)==18
 curve=s.Matrix([s.diff(v,u)for v in list(N)]+[s.diff(lam,u)if p==(1,4)else s.diff(mu,u)if p==(2,3)else 0 for p in powerpairs])
 assert D(curve)==D(sparse(payload['source_curve_s_derivative']))and D(J)*D(curve)==D(qdot*s.diff(a,u))
 assert rank(J.row_join(common))==rank(other.row_join(common))==24
 paths += [ROOT/'reports/Q_T01_MEMBERSHIP_20261103_v63_REVIEW.json',ROOT/'research/notes/Q_T01_SOURCE_DIFFERENTIAL_20261103_v63.md',ROOT/'research/notes/Q_T01_SOURCE_DIFFERENTIAL_20261103_v63_REVIEW.md']
 report=dict(id='Q_T01_SOURCE_DIFFERENTIAL_20261103_v63_REVIEW',all_passed=True,
  fresh_report_and_DECODED_matrices_match_compressed_mtime_allowed=True,
  independent_dual_number_matrix_power_all150_by31_differential=True,
  independent_all25_actual_GL_tensor_actions_and_both_equivariance_identities=True,
  exact_QQ_s_source_ranks=[24,24],joint_rank=30,common_rank=18,
  all_saved_generic_rank_minors_independently_checked=True,
  independently_solved_all5_polynomial_generator_compensations=True,
  actual_source_fiber_regular7_by_integrable_group_and_generator_families_accepted=True,
  actual_reduced_T01_tangent_at_every_finite_Qa_affine_at_least30_projective_at_least29=True,
  entire_finite_Qa_T01_singular_by_closed_tangent_rank_locus=True,
  target_smoothness_branch_count_normalization_or_special_source_fiber_regular_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
 print('PASS dual-number full150 derivative,24+24/joint30,actual sourcefiber7 and ALL-Qa T01 singularity')

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Independent complete66 double-flag linear spaces and global intersection audit."""
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
def transport(v,perm):
 out=s.zeros(150,v.cols)
 for(ab,m),i in IDS.items():
  ij=[perm[j]for j in ab];j=IDS[(tuple(sorted(ij)),tuple(sorted(perm[k]for k in m)))];sign=1 if ij[0]<ij[1]else-1
  out[j,:]=sign*v[i,:]
 return out
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 directory=ROOT/'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63'
 path=directory/'report.json';d=load(path);fresh=load(args.replay/'report.json')
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 tp=directory/'double-flag-table/report.json';table=load(tp);ft=load(args.replay/'double-flag-table/report.json')
 assert {k:v for k,v in table.items()if k!='seconds'}=={k:v for k,v in ft.items()if k!='seconds'}
 paths=[Path(__file__),path,tp,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py']
 for dd in [d,table]:
  for file,h in dd['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 # Independent double-coset partition of all120 permutations. The parabolic
 # is upper block triangular in the ordered blocks(y),(z,w),(x),(h).
 perms=set(itertools.permutations(range(5)));swap=[0,1,2,4,3];cosets=[]
 while perms:
  p=min(perms);orbit={p,tuple(swap[v]for v in p),tuple(p[swap[i]]for i in range(5)),tuple(swap[p[swap[i]]]for i in range(5))}
  assert orbit<=perms;perms-=orbit;cosets.append((p,orbit))
 assert len(cosets)==33 and sum(len(v)for p,v in cosets)==120
 assert {tuple(r['relative_permutation'])for r in table['records']}=={p for p,v in cosets}
 assert len({(tuple(r['relative_permutation']),r['CB_endpoint'])for r in table['records']})==66
 block=[3,2,0,1,1];positions={(i,j)for i in range(5)for j in range(5)if block[i]<=block[j]};assert len(positions)==16
 ca=load(ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json');CA=s.Matrix.hstack(*[s.Matrix(v)for v in ca['basis_canonical150']]);A=CA[:,[0,1,2,*range(4,20)]]
 cb=load(ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json');CB=s.Matrix.hstack(*[s.Matrix(v)for v in cb['enlarged_canonical_basis']])
 endpoints={'a_zero':CB[:,[*range(3,6),*range(6,22)]],'b_zero':CB[:,[*range(3),*range(6,22)]]};assert A.rank()==19
 same=0;exceptional=[];records=[];bigK=None;bigL=None
 for rec in table['records']:
  perm=rec['relative_permutation'];B=transport(endpoints[rec['CB_endpoint']],perm)
  k=rec['linear_intersection_dimension'];KA=s.Matrix.hstack(*[s.Matrix(v)for v in rec['intersection_basis_in_CA_coordinates']])if k else s.zeros(19,0)
  assert KA.rank()==k and A.row_join(-B).rank()==38-k
  L=A*KA;rows=list(B.T.rref()[1]);KB=B[rows,:].inv()*L[rows,:];assert B*KB==L
  orbit=25-len(positions&{(perm[i],perm[j])for i,j in positions});assert orbit==rec['full_flag_pair_orbit_dimension']
  isn=all(KA[2,j]==0 for j in range(k));assert isn==rec['entire_linear_family_lies_in_N_with_the_CA_flag']
  same+=isn
  if not isn:
   exceptional.append((perm,rec['CB_endpoint'],k))
   if perm==[3,1,2,0,4]:assert k==7 and KB[:2,:]==s.zeros(2,7)and any(KB[2,j]for j in range(7))
   elif perm==[1,3,2,0,4]:
    assert rec['CB_endpoint']=='b_zero'and k==8 and orbit==12;bigK,bigL=KA,L
   else:raise AssertionError(perm)
  records.append(dict(permutation=perm,endpoint=rec['CB_endpoint'],exact_QQ_intersection_dimension=k,exact_QQ_ambient_difference_rank=38-k,paired_flag_orbit_dimension=orbit,N_with_same_flag=isn))
 assert same==64 and len(exceptional)==2
 point=bigL*s.Matrix([i%7+1 for i in range(8)]);D=s.Matrix.hstack(bigL,*[action(i,j)*point for i,j in itertools.product(range(5),repeat=2)])
 minor=D[d['R_dimension_minor_rows'],d['R_dimension_minor_columns']];assert minor.det()==s.Rational(d['R_dimension_minor_determinant'])==3623878656000
 assert D.rank()==20
 h,x,y,z,w=s.symbols('h x y z w');X=[h,x,y,z,w];q=s.Matrix(s.symbols('q0:8'));params=bigK*q
 f=s.expand(params[0]*h*h+params[1]*h*x+params[2]*x*x)
 H=s.expand(sum(params[3+i]*s.sympify(mon)for i,mon in enumerate(ca['cubic_monomials'])))
 assert f==s.sympify(d['R_raw_f'])and H==s.sympify(d['R_raw_H'])
 witness=load(ROOT/'reports/AMBIENT_RANK2_EXCEPTION_NOT_N_20261103_v63_REVIEW.json')
 V=s.Matrix([s.sympify(v)for v in witness['raw_physical_vector_field']]);ey=s.Matrix([0,0,1,0,0]);P=ey*V.T-V*ey.T
 v=s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(s.prod(X[k]for k in m))for(ab,m),ind in IDS.items()for i,j in[ab]])
 assert canonical(v)==bigL*s.Matrix(d['R_witness_parameters'])
 paths += [ROOT/'reports/AMBIENT_RANK2_EXCEPTION_NOT_N_20261103_v63_REVIEW.json',
  ROOT/'research/notes/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_GLOBAL_INTERSECTION_GEOMETRY_REVIEW_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.md']
 result=dict(id='AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW',all_passed=True,
  fresh_all_QQ66_table_and_completion_replay_matches=True,independent_complete120_permutation_partition_double_cosets=33,
  all66_full_linear_intersection_bases_checked_over_QQ=True,all66_exact_QQ_kernel_dimensions_and_pair_orbit_dimensions_checked=True,
  records=records,same_CA_flag_N_cases=64,one_remaining_entire_CB_square_case=True,
  remaining_R_source_upper=19,independent_physical150_image_differential_rank=20,dimension_minor_determinant=str(minor.det()),
  separate_not_N_witness_actual_membership_in_R_verified=True,
  proper_covers_and_double_flag_exhaustion_geometry_accepted=True,
  entire_projective_reduced_intersection_exact_components=['N27','rank2_R19'],
  new_full_Poisson_component_or_global_scheme_primary_decomposition_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS independent complete66 QQ flag spaces; FULL global redintersection=N27 union R19')

if __name__=='__main__':main()

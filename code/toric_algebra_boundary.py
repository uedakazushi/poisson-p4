#!/usr/bin/env python3
"""Exact representatives for all twenty commutative unital C-algebras of rank5.

The completeness of the algebra list is the cited Poonen classification.
Poisson identities follow universally from commuting regular multiplication.
This is a catalog of the algebra-model boundary, not all of the T01 boundary.
"""
import argparse,hashlib,json
from pathlib import Path
from itertools import combinations_with_replacement,product
import sympy as s
from s61_model import S61,X,coefficients,quotient_coefficients,vertical_rref
from weight1_structure import J,COORDS,coefficient_map
ROOT=Path(__file__).resolve().parents[1]
x,y,z,w=s.symbols('x y z w')

def local_catalog():
 return [
  ('A1',1,[x],[x],[1]),
  ('A2',2,[x],[x*x],[1,x]),
  ('A3_curv',3,[x],[x**3],[1,x,x*x]),
  ('A3_sq',3,[x,y],[x*x,x*y,y*y],[1,x,y]),
  ('A4_curv',4,[x],[x**4],[1,x,x*x,x**3]),
  ('A4_21_r1',4,[x,y],[x*x,x*y,y**3],[1,x,y,y*y]),
  ('A4_21_r2',4,[x,y],[x*x,y*y],[1,x,y,x*y]),
  ('A4_sq',4,[x,y,z],[a*b for a,b in combinations_with_replacement([x,y,z],2)],[1,x,y,z]),
  ('A5_curv',5,[x],[x**5],[1,x,x*x,x**3,x**4]),
  ('A5_211_s2',5,[x,y],[x*x,x*y,y**4],[1,x,y,y*y,y**3]),
  ('A5_211_s1',5,[x,y],[x*x+y**3,x*y],[1,x,y,y*y,y**3]),
  ('A5_22_a',5,[x,y],[x*y,x**3,y**3],[1,x,y,x*x,y*y]),
  ('A5_22_b',5,[x,y],[x*x,x*y*y,y**3],[1,x,y,x*y,y*y]),
  ('A5_31_r1',5,[x,y,z],[x*x,y*y,x*y,x*z,y*z,z**3],[1,x,y,z,z*z]),
  ('A5_31_r2',5,[x,y,z],[x*x,y*y,z*z,x*y,x*z],[1,x,y,z,y*z]),
  ('A5_31_r3',5,[x,y,z],[x*y,x*z,y*z,x*x+y*y,x*x+z*z],[1,x,y,z,x*x]),
  ('A5_sq',5,[x,y,z,w],[a*b for a,b in combinations_with_replacement([x,y,z,w],2)],[1,x,y,z,w])]

def quotient_matrices(item):
 name,n,vs,ideal,basis=item
 gb=s.groebner(ideal,*vs,order='grevlex')
 assert gb.is_zero_dimensional
 pol=[s.Poly(b,*vs) for b in basis];mons=[p.monoms()[0] for p in pol]
 leading=[p.LM(order=gb.order).exponents for p in gb.polys]
 standards=[]
 for powers in product(range(n+1),repeat=len(vs)):
  if not any(all(a>=b for a,b in zip(powers,lead)) for lead in leading): standards.append(powers)
 assert len(standards)==n
 mons=sorted(standards,key=lambda powers:(sum(powers),powers))
 basis=[s.prod(a**b for a,b in zip(vs,powers)) for powers in mons]
 matrices=[]
 for a in basis:
  cols=[]
  for b in basis:
   rem=gb.reduce(s.sympify(a*b))[1];poly=s.Poly(rem,*vs)
   v=s.Matrix([poly.coeff_monomial(mon) for mon in mons])
   assert s.expand(rem-sum((v[i]*basis[i] for i in range(n)),s.S.Zero))==0
   cols.append(v)
  matrices.append(s.Matrix.hstack(*cols))
 assert matrices[0]==s.eye(n)
 assert all(a*b==b*a for a in matrices for b in matrices)
 return matrices

def partitions(items,total,minimum=0):
 if total==0: yield []
 for i in range(minimum,len(items)):
  if items[i][1]<=total:
   for rest in partitions(items,total-items[i][1],i):yield [i]+rest

def derivation_dimension(matrices,unit):
 n=len(matrices);v=s.symbols('d0:'+str(n*n));D=s.Matrix(n,n,v)
 equations=list(D*unit)
 for i in range(n):
  for j in range(i,n):
   equations.extend(D*matrices[i][:,j]-sum((D[k,i]*matrices[k][:,j]+D[k,j]*matrices[k][:,i] for k in range(n)),s.zeros(n,1)))
 linear=s.linear_eq_to_matrix(equations,v)[0]
 rank=len(linear.to_DM().rref()[1]);return n*n-rank

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
 local=local_catalog();local_matrices=[quotient_matrices(a) for a in local]
 cases=list(partitions(local,5));assert len(cases)==20
 R=s.Matrix([[1,2,3,5],[2,5,7,11],[3,7,13,17],[5,11,17,29]])
 assert R.det()!=0
 adapter=S61();results=[];variables=s.Matrix(X)
 vertical,pivots=vertical_rref()
 for case in cases:
  names=[local[k][0] for k in case];ns=[local[k][1] for k in case]
  raw=[];unit=s.zeros(5,1);offset=0
  for k in case:
   n=local[k][1];unit[offset]=1
   for mat in local_matrices[k]:
    block=s.zeros(5);block[offset:offset+n,offset:offset+n]=mat;raw.append(block)
   offset+=n
  traces=[s.trace(a) for a in raw]
  cols=[unit]+[s.eye(5)[:,i]-traces[i]*unit/5 for i in range(1,5)]
  T=s.Matrix.hstack(*cols);assert T.det()!=0
  T=T*s.diag(1,R)
  mats=[T.inv()*sum((T[a,i]*raw[a] for a in range(5)),s.zeros(5))*T for i in range(5)]
  assert mats[0]==s.eye(5) and all(s.trace(mats[i])==0 for i in range(1,5))
  regular=mats[1:];m=[a[1:5,1:5] for a in regular]
  q=s.Matrix.vstack(*(a[0,1:5] for a in regular));assert q==q.T
  assert q==s.Matrix(4,4,lambda i,j:s.trace(m[i]*m[j])/3)
  fields=s.Matrix.hstack(*(a*variables for a in regular))
  P=(fields*J*fields.T).applyfunc(s.expand)
  normalized=adapter.normalize(P);u=normalized['parameters']
  # Public adapter verifies full equations, reconstruction, and independent Schouten.
  check=adapter.verify_point(normalized)
  affine=P[1:5,1:5]-s.Matrix(X[1:5])*P[0,1:5]+P[0,1:5].T*s.Matrix(X[1:5]).T
  affine=affine.subs(X[0],1).applyfunc(s.expand)
  pf=s.factor(affine[0,1]*affine[2,3]-affine[0,2]*affine[1,3]+affine[0,3]*affine[1,2])
  norm=s.factor((s.eye(5)+sum((X[i+1]*regular[i] for i in range(4)),s.zeros(5))).det(method='domain-ge'))
  assert s.expand(pf-norm)==0
  beta=s.Matrix([(J*m[k]+m[k].T*J)[i,j] for (i,j),k in COORDS])
  d=derivation_dimension(raw,unit)
  tangent=[]
  for i,j in product(range(4),repeat=2):
   E=s.zeros(4);E[i,j]=1
   dm=[-E*m[k]+m[k]*E+sum((E[a,k]*m[a] for a in range(4)),s.zeros(4)) for k in range(4)]
   dq=s.Matrix(4,4,lambda a,b:s.trace(dm[a]*m[b]+m[a]*dm[b])/3)
   dregular=[s.Matrix.vstack(s.Matrix([[0]]).row_join(dq[k,:]),s.zeros(4,1).row_join(dm[k])) for k in range(4)]
   df=s.Matrix.hstack(*(a*variables for a in dregular))
   dP=(df*J*fields.T+fields*J*df.T).applyfunc(s.expand)
   tangent.append(s.Matrix(quotient_coefficients(dP)))
  differential_rank=len(s.Matrix.hstack(*tangent).to_DM().rref()[1])
  assert differential_rank<=16-d
  results.append({'algebra_type':' x '.join(names),'local_lengths':ns,'GL4_algebra_orbit_dimension':16-d,'automorphism_dimension':d,'sample_Poisson_map_differential_rank':differential_rank,'image_dimension_proved_equal_orbit_dimension':differential_rank==16-d,'sample_first_jet_coefficient_rank':len(coefficient_map(beta).to_DM().rref()[1]),'sample_trace_form_rank':len(q.to_DM().rref()[1]),'trace_zero_regular_multiplication':[[[str(a) for a in row] for row in mat.tolist()] for mat in regular],'sample_S61_coordinates':[str(a) for a in u],'sample_norm_factorization':str(norm),'sample_verification':check})
  print(names,'orbit',16-d,'image rank',differential_rank,flush=True)
 report={'id':'TORIC_ALGEBRA_BOUNDARY_20261101','all_passed':True,'field':'QQ for certificates, C for classification','local_classes_rank1_through5':[1,1,2,4,9],'total_rank5_algebra_types':20,'all_cases':results,'sample_change_of_trace_zero_basis':[[str(a) for a in row] for row in R.tolist()],'classification_primary_source':'https://math.mit.edu/~poonen/papers/dimension6.pdf, Table1 and section1; corrected2025version','geometric_source':'https://math.mit.edu/~poonen/papers/moduli.pdf, Remark1.2 and section6; Mazzola1980 Corollary4','proof_note':'research/notes/TORIC_ALGEBRA_BOUNDARY_20261101.md','limitations':['Twenty algebra types are not twenty new Poisson components.','Algebra-model boundary images lie in T01, but surjectivity onto the whole T01 boundary is not asserted.','A sample differential rank lower than the orbit dimension is only a lower bound on image dimension.','This does not classify PGL5 orbits of the resulting Poisson structures.'],'input_sha256':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/weight1_structure.py']}}
 args.output.parent.mkdir(parents=True,exist_ok=True)
 with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
if __name__=='__main__':main()

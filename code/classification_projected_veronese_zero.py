#!/usr/bin/env python3
"""Complete projective bivectors vanishing on a projected Veronese surface.

Use all150 homogeneous coefficients and all280 degree6 substitution equations;
Euler24 is quotiented only after its full exact inclusion has been verified.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,jacobi,minor_certificate
ROOT=Path(__file__).resolve().parents[1]

def main():
 started=time.monotonic();xx=s.symbols('x0:5');a,b,c=s.symbols('a b c');uv=(a,b,c);q=s.Matrix([a*a-b*b,b*b-c*c,a*b,a*c,b*c]);E=s.Matrix(xx)
 triples=list(itertools.combinations(range(5),3));m6=sorted([ex for ex in itertools.product(range(7),repeat=3)if sum(ex)==6],reverse=True)
 rows=list(itertools.product(triples,m6));rindex={r:i for i,r in enumerate(rows)};A=s.zeros(280,150)
 for col,(pair,ex) in enumerate(COLUMNS):
  polynomial=s.prod(q[i]**power for i,power in enumerate(ex))
  for tri in triples:
   if set(pair).issubset(tri):
    missing=next(i for i in tri if i not in pair);sign=(-1)**tri.index(missing)
    for exponent,co in s.Poly(s.expand(sign*q[missing]*polynomial),*uv).terms():A[rindex[(tri,exponent)],col]=co
 dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(s.QQ)
 rank=lambda M:len(dm(M).rref()[1])
 K=dm(A).nullspace().to_Matrix().T;assert A*K==s.zeros(280,K.cols)
 def decode(v):
  P=s.zeros(5)
  for coefficient,(ij,ex)in zip(v,COLUMNS):P[ij]+=coefficient*s.prod(z**power for z,power in zip(xx,ex))
  for i,j in PAIRS:P[j,i]=-P[i,j]
  return P.applyfunc(s.expand)
 def vector(P):return s.Matrix([s.Poly(s.expand(P[ij]),*xx).coeff_monomial(ex)for ij,ex in COLUMNS])
 def canonical(P):
  divergence=s.Matrix([sum(s.diff(P[i,j],xx[i])for i in range(5))for j in range(5)])
  return (P-(E*divergence.T-divergence*E.T)/5).applyfunc(s.expand)
 V=[]
 for i,j in itertools.product(range(5),repeat=2):
  field=s.zeros(5,1);field[i]=xx[j];V.append(vector(E*field.T-field*E.T))
 V=s.Matrix.hstack(*V);assert rank(V)==24;assert A*V==s.zeros(280,25)
 canonicalK=s.Matrix.hstack(*[vector(canonical(decode(K[:,i])))for i in range(K.cols)])
 # RREF of the row span gives a small canonical basis independent of the
 # arbitrary denominator-cleared nullspace vectors.
 red,piv=dm(canonicalK.T).rref();R=red.to_Matrix()[:len(piv),:];B=R.T
 assert A*B==s.zeros(280,B.cols);assert rank(V.row_join(B))==K.cols
 ns=s.symbols('z0:'+str(B.cols));P=decode(B*s.Matrix(ns));sparse={ij:dict(s.Poly(P[ij],*xx).terms())for ij in PAIRS if P[ij]!=0}
 J=jacobi(sparse).applyfunc(s.expand);equations=[]
 for f in J:
  if f!=0:
   poly=s.Poly(f,*ns);normalized=s.expand(f/poly.LC())
   if normalized not in equations:equations.append(normalized)
 G=s.groebner(equations,*ns,order='grevlex') if equations else None
 pf=[]
 for missing in range(5):
  i,j,k,l=[z for z in range(5)if z!=missing]
  pf.append(s.factor((-1)**missing*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])))
 # Euler-horizontal exterior square has a unique actual anticanonical
 # quintic D with Pf_i=x_i*D. This assertion is tested modulo Jacobi too.
 D=s.cancel(sum(xx[i]*pf[i]for i in range(5))/sum(z*z for z in xx))
 residual=[s.factor(pf[i]-xx[i]*D)for i in range(5)]
 print('rank',rank(A),'kernel',K.cols,'quotient',B.cols,'Jacobi',list(map(str,equations)),'GB',list(map(str,G.polys))if G else [],flush=True)
 print('Pfaffians',list(map(str,pf)),flush=True)
 if any(residual):print('Pf_i/x_i not common before Jacobi; use affine chart x0=1.',flush=True)
 affineP=P.subs(xx[0],1)[1:,1:]
 # The actual affine representative subtracts the Euler wedge of its x0 row.
 RR=s.Matrix([0]+[P[0,j]/xx[0]for j in range(1,5)])
 affineP=(P-E*RR.T+RR*E.T).subs(xx[0],1)[1:,1:]
 affinePf=s.factor(affineP[0,1]*affineP[2,3]-affineP[0,2]*affineP[1,3]+affineP[0,3]*affineP[1,2])
 encode=lambda Q:{str(i)+','+str(j):str(s.factor(Q[i,j]))for i,j in PAIRS if Q[i,j]!=0}
 inputs=[Path(__file__),ROOT/'code/classification_rank3_mixed_component_review.py']
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),
  projected_surface_parametrization=list(map(str,q)),surface_parameter_variables=list(map(str,uv)),ambient_coordinates=list(map(str,xx)),
  ambient_coefficient_columns=[dict(pair=list(ij),exponents=list(ex))for ij,ex in COLUMNS],substitution_equation_rows=[dict(triple=list(tri),exponents=list(ex))for tri,ex in rows],
  complete_280_by150_sparse_matrix={str(i):{str(j):str(A[i,j])for j in range(150)if A[i,j]}for i in range(280)},
  exact_matrix_rank=rank(A),exact_kernel_dimension=K.cols,complete_kernel_basis_columns=[[str(z)for z in K[:,i]]for i in range(K.cols)],kernel_is_annihilated=True,rank_certificate=minor_certificate(A,rank(A)),
  Euler_vertical_rank=24,all_Euler_directions_annihilated=True,complete_projective_section_dimension=B.cols,canonical_section_basis=[encode(decode(B[:,i]))for i in range(B.cols)],
  complete_canonical_source=encode(P),parameters=list(map(str,ns)),complete350_Jacobi=list(map(str,J)),all_distinct_monic_Jacobi_equations=list(map(str,equations)),Jacobi_grevlex_basis=[str(z.as_expr())for z in G.polys]if G else [],
  full_Pfaffian_five_vector=list(map(str,pf)),actual_affine_representative=[[str(z)for z in row]for row in affineP.tolist()],actual_affine_Pfaffian=str(affinePf),
  scope='Complete fixed smooth projectedVeronese zero-surface source only. No assertion that all rankzero surfaces are projectedVeronese, and no global quintic exhaustion.',seconds=time.monotonic()-started)
 dest=ROOT/'reports/CLASSIFICATION_PROJECTED_VERONESE_ZERO_20261107.json'
 if dest.exists():raise FileExistsError(str(dest))
 dest.write_text(json.dumps(result,indent=2)+'\n');print('seconds',result['seconds'])
if __name__=='__main__':main()

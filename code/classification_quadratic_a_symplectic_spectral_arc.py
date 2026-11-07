#!/usr/bin/env python3
"""The complete second-order spectral coefficient forced by the pole-two frame.

Uses the three actual whole sources from the preceding five-parameter theorem;
no higher Jacobi jet or arbitrary product coefficient is introduced.
"""
from pathlib import Path
import hashlib,itertools as it,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_quadratic_a_symplectic_source import anti,jac,encode
from classification_rank3_mixed_component_review import jacobi as projective_jacobi,PAIRS
ROOT=Path(__file__).resolve().parents[1]
def main():
 begun=time.monotonic();u1,u2,v1,v2=xx=s.symbols('u1 u2 v1 v2');ep,lam,r0,r1,c0,c1=s.symbols('epsilon lambda r0 r1 c0 c1');a,f,p,q,r=s.symbols('a f p q r');b=q+r;e=p+q;J=anti({(0,2):1,(1,3):1});Q=anti({(0,1):p*u1*u2,(0,2):a*u1*v1+b*u2*v2,(0,3):q*u1*v2,(1,2):(p+q+r)*u2*v1,(1,3):e*u1*v1+f*u2*v2,(2,3):r*v1*v2});alpha=lam-s.Rational(1,2);beta=-lam-s.Rational(1,2);target=J+anti({(0,1):alpha*u2*u2,(0,2):beta*u2*v1,(0,3):-u2*v1-beta*u1*v1-alpha*u2*v2,(1,3):alpha*u2*v1,(2,3):beta*v1*v1});G=s.Matrix([[1,1/ep],[0,1]]);T=s.diag(G,G.inv().T);back=T.inv()*s.Matrix(xx);spectral=ep*r0+ep*ep*r1;period=c0+ep*c1
 sources=[('separated',{a:-spectral,f:-spectral-ep,p:-period,q:period,r:-period}),('rank3_rank5',{a:-spectral,f:-(spectral+ep)/2,p:-period,q:period,r:-(spectral+ep)/2-period}),('rank5_rank3',{a:-spectral/2,f:-spectral-ep,p:-spectral/2-period,q:period,r:-period})];records=[];dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ)
 for name,binding in sources:
  source=(J+Q).subs(binding,simultaneous=True);assert jac(source,xx)==[0]*4;family=(T*source*T.T).subs(dict(zip(xx,back)),simultaneous=True).applyfunc(s.cancel);equations=[];labels=[]
  for i,j in it.combinations(range(4),2):
   for(ex,*powers),co in s.Poly(s.expand(ep**2*(family-target)[i,j]),ep,*xx).terms():
    if ex<=2:equations.append(co);labels.append(dict(pair=[i,j],epsilon_power=ex-2,exponents=powers))
  A,v=s.linear_eq_to_matrix(equations,[r0,r1,c0,c1]);pc=list(dm(A).rref()[1]);pr=list(dm(A[:,pc].T).rref()[1]);inv=dm(A[pr,pc]).inv().to_Matrix();res=(v-A[:,pc]*inv*v[pr,:]).applyfunc(s.cancel);bad=next((i for i,z in enumerate(res)if z!=0),None);row=dict(name=name,actual_source_binding={str(k):str(z)for k,z in binding.items()},linear_shape=list(A.shape),rank=len(pc),matrix=[[str(z)for z in row]for row in A.tolist()],rhs=list(map(str,v)),row_labels=labels)
  if bad is not None:
   L=s.zeros(1,A.rows);L[bad]=1
   for j,k in enumerate(pr):L[k]-=(A[bad,pc]*inv)[j]
   assert L*A==s.zeros(1,4)and(L*v)[0]!=0;row.update(status='POLE_BOUNDED_SPECTRAL_ARC_INCONSISTENT',rank_augmented=len(pc)+1,left_functional={str(i):str(z)for i,z in enumerate(L)if z!=0},left_rhs=str((L*v)[0]))
  else:
   values=s.zeros(4,1)
   for j,k in enumerate(pc):values[k]=(inv*v[pr,:])[j]
   sol=dict(zip([r0,r1,c0,c1],values));final=family.subs(sol).applyfunc(s.cancel);assert all(s.denom(z).subs(ep,0)!=0 for z in final);assert(final.subs(ep,0)-target).applyfunc(s.expand)==s.zeros(4);assert jac(final,xx)==[0]*4;h=s.Symbol('h');hh=(h,*xx);Ph=s.zeros(5)
   for i,j in it.combinations(range(4),2):Ph[i+1,j+1]=sum(c*h**(2-sum(ex))*s.prod(z**k for z,k in zip(xx,ex))for ex,c in s.Poly(final[i,j],*xx).terms());Ph[j+1,i+1]=-Ph[i+1,j+1]
   sparse={ij:dict(s.Poly(s.expand(Ph[ij]),*hh).terms())for ij in PAIRS if Ph[ij]!=0};assert projective_jacobi(sparse).applyfunc(s.expand)==s.zeros(350,1);pf=s.factor(final[0,1]*final[2,3]-final[0,2]*final[1,3]+final[0,3]*final[1,2]);row.update(status='ALL_LAMBDA_ACTUAL_POISSON_CURVE_FOUND',spectral_period_solution={str(k):str(z)for k,z in sol.items()},actual_regular_affine_family=encode(final),actual_regular_homogeneous_family=encode(Ph),actual_affine_Pfaffian=str(pf),all350_Jacobi_exact_zero=True,endpoint_all_lambda_exact=True,true_source_periods={str(k):str(s.factor(z.subs(sol)))for k,z in binding.items()},independent_review_required=True)
  records.append(row)
 paths=[Path(__file__),ROOT/'code/classification_quadratic_a_symplectic_source.py',ROOT/'reports/CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_SOURCE_20261107.json',ROOT/'code/classification_rank3_mixed_component_review.py'];out=dict(task_id='E_multiple_hyperplane',all_passed=True,status='WHOLE_QUADRATIC_A_CONTAINMENT_CANDIDATE'if any(row['status']=='ALL_LAMBDA_ACTUAL_POISSON_CURVE_FOUND'for row in records)else'COMPLETE_SPECTRAL_ARC_NEGATIVE_STOP_ROUTE',whole_three_parent_sources_dependency='Prior five-coefficient exact Jacobi source theorem; each parent belongs to known H+2Q closed union by generic irreducible quadratic factors and closure.',apriori_degree_reason='Block symplectic frame has inverse pole1. After exact source identities, its linear spectral coefficient has pole2 and the true free period tensor pole1. The eigenvalue second coefficient r1 affects the endpoint without changing the nonzero Jordan modular endpoint. No higher spectral or period coefficient can affect poles or endpoint.',generic_source_magnetic_complements='Prior universal all6mag normalizer on dense alpha*beta*I open and irreducibility of the entire accepted7source. Positive all-lambda membership returns all complements by closedness.',curves=records,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_bytes=Path(__file__).read_text(),seconds=time.monotonic()-begun);path=ROOT/'reports/CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_SPECTRAL_ARC_20261107.json';assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status=out['status'],curves=[dict(name=r['name'],status=r['status'],solution=r.get('spectral_period_solution'),failure=r.get('left_rhs'))for r in records],seconds=out['seconds'])))
if __name__=='__main__':main()

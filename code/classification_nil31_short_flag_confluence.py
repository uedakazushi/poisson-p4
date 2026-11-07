#!/usr/bin/env python3
"""One actual independent-three-H+quadric flag confluence for nil31 short 2/3.

The four source coefficients are genuine independent parameters of a universal
Poisson cone. No coefficient of a product of parameters is freed artificially.
Failure concerns only this fixed geometric flag, not any component membership.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2))
def wedge(v,w):return v*w.T-w*v.T
def canonical(P,xx):
 d=s.Matrix([sum(s.diff(P[j,i],xx[j])for j in range(5))for i in range(5)])
 return (P-wedge(s.Matrix(xx),d)/5).applyfunc(s.expand)
def main():
 begun=time.monotonic();xx=s.symbols('h x y z w');h,x,y,z,w=xx;e=s.Symbol('epsilon');ss=s.symbols('hs X Y Z W');hs,X,Y,Z,W=ss
 D1=s.Matrix([0,X,0,0,W]);D2=s.Matrix([0,0,Y,Z,0]);N=s.Matrix([0,0,0,Y,-X]);vZ=s.eye(5)[:,3];vW=s.eye(5)[:,4]
 bases=[wedge(D1,D2),(X*Z+Y*W)*wedge(vZ,vW),wedge(N,D2),wedge(D1,N)]
 L=s.Matrix([h+e*y,h,z,y+e*x,-w]).jacobian(s.Matrix(xx));Li=L.inv();det=s.factor(L.det());assert det!=0
 pulled=[canonical((Li*P.subs(dict(zip(ss,L*s.Matrix(xx))),simultaneous=True)*Li.T).applyfunc(s.expand),xx)for P in bases]
 target=s.zeros(5);br=[h*z-h*y+2*z*w-y*y/3,h*z-2*y*z/3,h*h+y*w/3,-2*z*z/3,(h*y+z*w)/3,2*h*z/3]
 for (i,j),a in zip(itertools.combinations(range(1,5),2),br):target[i,j]=a;target[j,i]=-a
 target=canonical(target,xx);mons=list(itertools.combinations_with_replacement(range(5),2));cs=s.symbols('c0:20');periods=[sum(cs[5*i+k]*e**k for k in range(5))for i in range(4)]
 candidate=sum((a*P for a,P in zip(periods,pulled)),s.zeros(5));diff=(candidate-target).applyfunc(s.expand);rows=[];rhs=[];labels=[]
 # Canonical Euler quotient has 126 effective coordinates; retain all150.
 for i,j in PAIRS:
  poly=s.Poly(s.expand(e**4*diff[i,j]),e,*xx)
  assert all(sum(ex[1:])==2 for ex,c in poly.terms())
  for (k,*ex),coef in poly.terms():
   if k<=4:
    M,b=s.linear_eq_to_matrix([coef],cs);rows.append(list(M));rhs.append(b[0]);labels.append(dict(pair=[i,j],epsilon_power=k-4,exponents=ex))
 A=s.Matrix(rows);b=s.Matrix(rhs);dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(s.QQ)
 piv=list(dm(A).rref()[1]);rr=list(dm(A[:,piv].T).rref()[1]);inv=dm(A[rr,piv]).inv().to_Matrix();res=(b-A[:,piv]*inv*b[rr,:]).applyfunc(s.cancel);bad=next((i for i,a in enumerate(res)if a),None)
 out=dict(task_id='E_multiple_hyperplane',all_passed=True,actual_source_configuration='Q=XZ+YW, independent hs,X,Y',source_four_actual_bases=[{str(i)+','+str(j):str(P[i,j])for i,j in PAIRS if P[i,j]!=0}for P in bases],source_coordinate_forms=list(map(str,L*s.Matrix(xx))),source_frame_determinant=str(det),inverse_frame=[[str(a)for a in row]for row in Li.tolist()],source_periods=list(map(str,periods)),source_period_degree_bound=4,bound_reason='This chosen polynomial actual subfamily has degree4; inverse frame pole2 gives bivector pole4. No completeness claim for arbitrary flags or source arcs.',full150_projective_Euler_coordinates=True,linear_system_shape=list(A.shape),rank=len(piv),row_labels=labels,full_matrix=[[str(a)for a in row]for row in A.tolist()],full_rhs=list(map(str,b)),minor_rows=rr,minor_columns=piv,minor_determinant=str(dm(A[rr,piv]).det()),target_anchor='t=d=1 of the accepted whole-source canonical nil31 short 2/3 family',target_canonical_tensor={str(i)+','+str(j):str(target[i,j])for i,j in PAIRS if target[i,j]!=0})
 if bad is not None:
  lam=s.zeros(1,A.rows);lam[bad]=1
  for j,k in enumerate(rr):lam[k]-=(A[bad,piv]*inv)[j]
  assert lam*A==s.zeros(1,20) and (lam*b)[0]!=0
  out.update(status='THIS_ACTUAL_FIXED_FLAG_CLASS_HAS_NO_REGULAR_TARGET_LIMIT',rank_augmented=len(piv)+1,left_inconsistency_certificate={str(j):str(a)for j,a in enumerate(lam)if a},left_rhs=str((lam*b)[0]),left_M_zero=True,whole_source_membership='UNKNOWN; this fixed flag failure proves no component exclusion',next_route='Do not increase Taylor degree for this flag. Diagnose flag geometry and change actual family or classification route.')
 else:
  sol=s.zeros(20,1)
  for j,k in enumerate(piv):sol[k]=(inv*b[rr,:])[j]
  bind=dict(zip(cs,sol));period_values=[s.expand(p.subs(bind))for p in periods];family=candidate.subs(bind).applyfunc(s.cancel)
  assert all(s.denom(a).subs(e,0)!=0 for a in family)
  assert(family.subs(e,0)-target).applyfunc(s.cancel)==s.zeros(5)
  triples={inds:s.expand(sum(family[i,l]*s.diff(family[j,k],xx[l])+family[j,l]*s.diff(family[k,i],xx[l])+family[k,l]*s.diff(family[i,j],xx[l])for l in range(5)))for inds in itertools.combinations(range(5),3)for i,j,k in [inds]}
  J=[s.expand(sum((-1)**k*xx[inds[k]]*triples[inds[:k]+inds[k+1:]]for k in range(4)))for inds in itertools.combinations(range(5),4)]
  assert J==[0]*5
  out.update(status='ACTUAL_POISSON_CONFLUENCE_FOUND',chosen_periods=list(map(str,period_values)),whole_actual_homogeneous_family={str(i)+','+str(j):str(family[i,j])for i,j in PAIRS if family[i,j]!=0},all350_projective_Jacobi_zero=True,target_endpoint_exact=True,independent_review_required=True)
 paths=[Path(__file__),ROOT/'reports/CLASSIFICATION_NIL31_SHORT_TWO_THIRDS_CANONICAL_20261107.json',ROOT/'research/classification/THREE_HYPERPLANES_QUADRIC_INDEPENDENT.md',ROOT/'reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_REVIEW_20261107.json'];out['input_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths};out['producer_bytes']=Path(__file__).read_text();out['seconds']=time.monotonic()-begun
 path=ROOT/'reports/CLASSIFICATION_NIL31_SHORT_FLAG_CONFLUENCE_20261107.json';assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in ['status','linear_system_shape','rank','seconds']}))
if __name__=='__main__':main()

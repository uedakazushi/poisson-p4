#!/usr/bin/env python3
"""Complete five-parameter semisimple symplectic source and genuine contractions.

No local Mora, unit elimination or jet extension. Every generic parent is one
of three actual whole Poisson sources, with factorized H+2Q divisor.
"""
from pathlib import Path
import hashlib,itertools as it,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
def anti(entries,n=4):
 P=s.zeros(n)
 for(i,j),f in entries.items():P[i,j]=f;P[j,i]=-f
 return P
def jac(P,xx):return[s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l])for l in range(len(xx))))for i,j,k in it.combinations(range(len(xx)),3)]
def encode(P):return{str(i)+','+str(j):str(s.factor(P[i,j]))for i,j in it.combinations(range(P.rows),2)if P[i,j]!=0}
def main():
 started=time.monotonic();u1,u2,v1,v2=xx=s.symbols('u1 u2 v1 v2');a,f,p,q,r=s.symbols('a f p q r');b=q+r;e=p+q;J=anti({(0,2):1,(1,3):1});Q=anti({(0,1):p*u1*u2,(0,2):a*u1*v1+b*u2*v2,(0,3):q*u1*v2,(1,2):(p+q+r)*u2*v1,(1,3):e*u1*v1+f*u2*v2,(2,3):r*v1*v2});P=J+Q
 assert[j1-j2 for j1,j2 in zip(jac(P,xx),jac(Q,xx))]==[0]*4
 coefficients=[];labels=[]
 for tri,expr in zip(it.combinations(range(4),3),jac(P,xx)):
  for ex,c in s.Poly(expr,*xx).terms():
   if c!=0:coefficients.append(c);labels.append(dict(triple=list(tri),exponents=list(ex)))
 generators=[b*e,b*(b-f),e*(e-a)];params=(a,f,p,q,r);mons=sorted(set(ex for z in coefficients+generators for ex,c in s.Poly(z,*params).terms()),reverse=True);vec=lambda z:s.Matrix([[s.Poly(z,*params).coeff_monomial(m)for m in mons]])
 C=s.Matrix.vstack(*map(vec,coefficients));G=s.Matrix.vstack(*map(vec,generators));dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(s.QQ);assert len(dm(C).rref()[1])==len(dm(G).rref()[1])==len(dm(C.col_join(G)).rref()[1])==3
 # Explicit two-direction constant coefficient certificates for the ideals.
 cols=list(dm(G).rref()[1]);GC=G[:,cols]
 to_g=(C[:,cols]*GC.inv());assert to_g*G==C
 independent_rows=list(dm(C.T).rref()[1]);CB=C[independent_rows,cols];to_c_small=G[:,cols]*CB.inv();to_c=s.zeros(3,C.rows)
 for j,k in enumerate(independent_rows):to_c[:,k]=to_c_small[:,j]
 assert to_c*C==G
 pf=s.factor(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2]);assert s.expand(pf+1+(a+e)*u1*v1+(b+f)*u2*v2+a*e*(u1*v1)**2+a*f*u1*u2*v1*v2+b*f*(u2*v2)**2)==0
 # The entire seven-parameter target, all magnetic coefficients retained.
 lam=s.Symbol('lambda');alpha=lam-s.Rational(1,2);beta=-lam-s.Rational(1,2);q1,q2,x,y=old=s.symbols('q1 q2 x y');mag=s.symbols('m0:6');R=sum(c*m for c,m in zip(mag,[1,x,y,x*x,x*y,y*y]));target_old=anti({(0,1):R+beta*y*q1-alpha*x*q2,(0,2):alpha*x*x,(0,3):1+beta*x*y,(1,2):1+alpha*x*y,(1,3):beta*y*y})
 source_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json';original=json.loads(source_path.read_text())['cases'][7];assert original['name']=='nil22_image_quadratic_A';env=dict(zip(['h','x1','x2','x3','x4'],[s.Integer(1),q1,q2,x,y]));env.update(lam=lam);env.update({str(m):m for m in mag})
 for i,j in it.combinations(range(4),2):assert s.expand(s.sympify(original['full_six_magnetic_source'].get(str(i+1)+','+str(j+1),'0'),locals=env)-target_old[i,j])==0
 I=mag[0]+mag[4];normal=s.Matrix([(q1+mag[2]/beta+mag[5]*y/(2*beta))/I,(q2-mag[1]/alpha-mag[3]*x/(2*alpha)+mag[4]*y)/I,x,y]);N=normal.jacobian(s.Matrix(old));shift=normal-N*s.Matrix(old);back=N.inv()*(s.Matrix(old)-shift);normalized=(I*N*target_old*N.T).subs(dict(zip(old,back)),simultaneous=True).applyfunc(s.cancel);expected=target_old.subs(dict(zip(mag,[1,0,0,0,0,0])));assert(normalized-expected).applyfunc(s.cancel)==s.zeros(4)
 Cmap=s.Matrix([q1,x,y,y-q2]);Cmat=Cmap.jacobian(s.Matrix(old));target=(Cmat*expected*Cmat.T).subs(dict(zip(old,Cmat.inv()*s.Matrix(xx))),simultaneous=True).applyfunc(s.expand);assert(target-J-anti({(0,1):alpha*u2*u2,(0,2):beta*u2*v1,(0,3):-u2*v1-beta*u1*v1-alpha*u2*v2,(1,3):alpha*u2*v1,(2,3):beta*v1*v1})).applyfunc(s.expand)==s.zeros(4)
 assert jac(target,xx)==[0]*4
 ep,rho,c0,c1=s.symbols('epsilon rho c0 c1');period=c0+ep*c1;G2=s.Matrix([[1,1/ep],[0,1]]);T=s.diag(G2,G2.inv().T);assert T*J*T.T==J;inverse=T.inv()*s.Matrix(xx)
 sources=[('separated',dict(a=-ep*rho,f=-ep*(rho+1),p=-period,q=period,r=-period)),('rank3_rank5',dict(a=-ep*rho,f=-ep*(rho+1)/2,p=-period,q=period,r=-ep*(rho+1)/2-period)),('rank5_rank3',dict(a=-ep*rho/2,f=-ep*(rho+1),p=-ep*rho/2-period,q=period,r=-period))];records=[]
 for name,values in sources:
  bind={str(k):v for k,v in values.items()};sub={z:bind[str(z)]for z in params};source=P.subs(sub,simultaneous=True);assert jac(source,xx)==[0]*4;factor=s.factor((P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2]).subs(sub,simultaneous=True));family=(T*source*T.T).subs(dict(zip(xx,inverse)),simultaneous=True).applyfunc(s.cancel);difference=(family-target).applyfunc(s.cancel);equations=[];rowlabels=[]
  for i,j in it.combinations(range(4),2):
   for(ex,*powers),coef in s.Poly(s.expand(ep**4*difference[i,j]),ep,*xx).terms():
    if ex<=4:equations.append(coef);rowlabels.append(dict(pair=[i,j],epsilon_power=ex-4,exponents=powers))
  A,v=s.linear_eq_to_matrix(equations,[rho,c0,c1]);piv=list(dm(A).rref()[1]);rows=list(dm(A[:,piv].T).rref()[1]);inv=dm(A[rows,piv]).inv().to_Matrix();res=(v-A[:,piv]*inv*v[rows,:]).applyfunc(s.cancel);bad=next((k for k,z in enumerate(res)if z!=0),None);record=dict(name=name,source_binding={str(k):str(v)for k,v in sub.items()},actual_source_Pfaffian=str(factor),all_source_Jacobi_zero=True,linear_shape=list(A.shape),rank=len(piv),linear_matrix=[[str(z)for z in row]for row in A.tolist()],linear_rhs=list(map(str,v)),row_labels=rowlabels)
  if bad is not None:
   L=s.zeros(1,A.rows);L[bad]=1
   for j,k in enumerate(rows):L[k]-=(A[bad,piv]*inv)[j]
   assert L*A==s.zeros(1,3)and(L*v)[0]!=0;record.update(status='FIXED_ACTUAL_SOURCE_CONTRACTION_NOT_MATCHED',left_functional={str(i):str(z)for i,z in enumerate(L)if z!=0},left_rhs=str((L*v)[0]))
  else:
   sol=s.zeros(3,1)
   for j,k in enumerate(piv):sol[k]=(inv*v[rows,:])[j]
   values=dict(zip([rho,c0,c1],sol));final=family.subs(values).applyfunc(s.cancel);assert all(s.denom(z).subs(ep,0)!=0 for z in final);assert(final.subs(ep,0)-target).applyfunc(s.cancel)==s.zeros(4);assert jac(final,xx)==[0]*4
   record.update(status='ACTUAL_FULL_LAMBDA_POISSON_CURVE_FOUND',parameter_solution={str(k):str(z)for k,z in values.items()},whole_actual_regular_family=encode(final),all4_Jacobi_zero=True,full_lambda_endpoint_exact=True,generic_known_component_union='Accepted entire H+two distinct irreducible quadric locus, extended to these irreducible polynomial parent sources by closure.',independent_review_required=True)
  records.append(record)
 input_paths=[Path(__file__),source_path,ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_MAGNETIC_NORMALFORMS_20261107.json',ROOT/'research/classification/NIL22_IMAGE_SOURCES_AND_COMPONENTS.md',ROOT/'research/classification/TWO_QUADRICS_SUPPORT.md'];out=dict(task_id='E_multiple_hyperplane',all_passed=True,status='ACTUAL_WHOLE_SOURCE_CONTAINMENT_CANDIDATE'if any(r['status']=='ACTUAL_FULL_LAMBDA_POISSON_CURVE_FOUND'for r in records)else'COMPLETE_DISTINCT_EIGEN_SOURCE_CLASSIFIED_FIXED_CONTRACTIONS_NEGATIVE',symplectic_J=encode(J),five_parameter_Q=encode(Q),complete_Jacobi_coefficients=list(map(str,coefficients)),Jacobi_coefficient_labels=labels,three_ideal_generators=list(map(str,generators)),Jacobi_to_generators=[[str(z)for z in row]for row in to_g.tolist()],generators_to_Jacobi=[[str(z)for z in row]for row in to_c.tolist()],whole_semisimple_source_support_union=['b=e=0','e=0,f=b','b=0,a=e'],nonregular_walls='All distinct eigenvalues including ratio-1 andzero eigenvalue allowed in the fiveparameter source. Repeated eigenvalue scalar is a closed block condition and cannot approach nonzero Jordan target via block symplectic conjugacy. Nonsemisimple Jordan parents remain an explicit separate stratum, not classified here.',all_six_magnetic_original_source_exact=True,all_six_magnetic_normalization_exact=True,normalization_localization='alpha*beta*I !=0; dense open of the accepted irreducible entire7source, so positive entirelambda membership returns every complement by source closure.',full_normalization_map=list(map(str,normal)),full_normalization_inverse=list(map(str,back)),normalization_determinant=str(s.factor(N.det())),symplectic_standardization=list(map(str,Cmap)),target_standard=encode(target),actual_symplectic_contraction_matrix=[[str(z)for z in row]for row in T.tolist()],actual_parent_period='q=c0+epsilon*c1, actual single source period; no products freed',curves=records,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in input_paths},producer_bytes=Path(__file__).read_text(),seconds=time.monotonic()-started);path=ROOT/'reports/CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_SOURCE_20261107.json';assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status=out['status'],curves=[dict(name=r['name'],status=r['status'],solution=r.get('parameter_solution'),failure=r.get('left_rhs'))for r in records],seconds=out['seconds'])))
if __name__=='__main__':main()

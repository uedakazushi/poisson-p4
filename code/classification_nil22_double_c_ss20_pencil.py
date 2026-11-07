#!/usr/bin/env python3
"""Full90 affinequadratic SS20-pencil entrance for wholegeneric doubleC source.
No restrictive NwedgeY ansatz and no arbitrary magneticpoint.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();q1,q2,x,y=s.symbols('q1 q2 x y');xx=(q1,q2,x,y);a,b=s.symbols('a b');zz=s.symbols('z0:90');pairs=list(itertools.combinations(range(4),2));mons=[s.Integer(1),*xx,*[xx[i]*xx[j]for i,j in itertools.combinations_with_replacement(range(4),2)]];triples=list(itertools.combinations(range(4),3));mon3=[s.prod(xx[i]for i in factors)for d in range(4)for factors in itertools.combinations_with_replacement(range(4),d)]
 source_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json';source=json.loads(source_path.read_text())['cases'][2];assert source['name']=='nil22_image_double_C';env=dict(zip(['h','x1','x2','x3','x4'],[1,*xx]));env.update(dict(zip(['m'+str(i)for i in range(6)],[a,b,0,1,0,0])));P=s.zeros(4)
 for ij,text in source['full_six_magnetic_source'].items():
  i,j=map(int,ij.split(','));P[i-1,j-1]=s.sympify(text,locals=env);P[j-1,i-1]=-P[i-1,j-1]
 def jac(Q):return [s.expand(sum(Q[i,l]*s.diff(Q[j,k],xx[l])+Q[j,l]*s.diff(Q[k,i],xx[l])+Q[k,l]*s.diff(Q[i,j],xx[l])for l in range(4)))for i,j,k in triples]
 def modular(Q):return s.Matrix([sum(s.diff(Q[i,j],xx[i])for i in range(4))for j in range(4)])
 assert jac(P)==[0]*4;N=modular(P);A0=N.jacobian(s.Matrix(xx));assert N==s.Matrix([x,y,0,0]);columns=[];anti=[];basis=[]
 for ij,mon in itertools.product(pairs,mons):
  Q=s.zeros(4);Q[ij]=mon;Q[ij[1],ij[0]]=-mon;basis.append(Q);JJ=jac(P+Q);columns.append(s.Matrix([s.Poly(f,*xx).coeff_monomial(m)for f in JJ for m in mon3]));B=modular(Q).jacobian(s.Matrix(xx));anti.append(s.Matrix(list(A0*B+B*A0)))
 D=s.Matrix.hstack(*columns);H=s.Matrix.hstack(*anti);A=D.col_join(H);rhs=s.zeros(D.rows,1).col_join(s.Matrix(list(s.eye(4))));K=s.QQ.frac_field(a,b);dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(K);Am=dm(A);Aug=dm(A.row_join(rhs));r=len(Am.rref()[1]);reduced,pivots=Aug.rref();ar=len(pivots);print('full90linear',A.shape,r,ar,round(time.monotonic()-started,3),flush=True)
 result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),source_path,ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_MAGNETIC_NORMALFORMS_20261107.json',ROOT/'research/classification/HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_SUPPORT.md']},producer_bytes=Path(__file__).read_text(),canonical_generic_magnetic_F0='a+b*x+x^2',canonical_generic_normalization='Originalfull3residualmagneticvector(I0,I1,I2), onI2!=0 use q_new=q/I2,tensorscaleI2, giving(a,b,1). Wholeirreducsource boundary followsclosedness ONLYafterpositiveknownmembership.',all90_affine_quadratic_unknowns=True,linear_shape=list(A.shape),linear_rank=r,augmented_rank=ar,complete_linear_matrix=[[str(z)for z in row]for row in A.tolist()],linear_rhs=list(map(str,rhs)),basis_pairs=[list(ij)for ij in pairs],basis_monomials=list(map(str,mons)),scope='FullcrossJacobi andN0B+BN0=I necessaryforactualsemisimpleoppositepair pencil. Ifconsistent, thecanonicalparticular is testedexactly. Failureofoneparticular isnotfullpencilnonexistence. No membership withoutallselfJacobi/Bsquare/Pfzero hypotheses.',seconds=time.monotonic()-started)
 if ar>r:
  left=Am.transpose().nullspace(divide_last=True).to_Matrix();row=next(left.row(i)for i in range(left.rows)if s.cancel((left.row(i)*rhs)[0])!=0);assert(row*A).applyfunc(s.cancel)==s.zeros(1,90);result.update(no_canonical_SS20_transverse_pencil=True,left_functional={str(i):str(z)for i,z in enumerate(row)if z!=0},nonzero_pairing=str(s.cancel((row*rhs)[0])))
 else:
  RREF=reduced.to_Matrix();part=s.zeros(90,1)
  for i,pivot in enumerate(pivots):part[pivot]=RREF[i,90]
  assert(A*part-rhs).applyfunc(s.cancel)==s.zeros(A.rows,1);Q=sum([co*B for co,B in zip(part,basis)],s.zeros(4)).applyfunc(s.cancel);J=[s.factor(s.cancel(f))for f in jac(Q)];Z=modular(Q).applyfunc(s.cancel);B=Z.jacobian(s.Matrix(xx));mu=s.factor(s.trace(B*B)/4);SS=(B*B-mu*s.eye(4)).applyfunc(s.cancel)==s.zeros(4)
  result.update(complete_free_unknown_indices=[i for i in range(90)if i not in pivots],canonical_particular_Q={str(i)+','+str(j):str(Q[i,j])for i,j in pairs if Q[i,j]!=0},canonical_particular_all4_selfJacobi=list(map(str,J)),canonical_particular_Poisson=not any(J),canonical_particular_modular=list(map(str,Z)),canonical_particular_modular_square_scalar=SS,modular_square_scalar=str(mu),generic_SS20_actual_family_found=not any(J)and SS)
  if not any(J)and SS:
   e=s.Symbol('epsilon');family=P+e*Q;zero=(-(A0+e*B).inv()*(e*Z.subs(dict.fromkeys(xx,0)))).applyfunc(s.cancel);C=family.subs(dict(zip(xx,zero)),simultaneous=True).applyfunc(s.cancel);pf=s.factor(s.cancel(C[0,1]*C[2,3]-C[0,2]*C[1,3]+C[0,3]*C[1,2]));assert pf!=0;result.update(modular_zero=list(map(str,zero)),constant_Pf_at_modular_zero=str(pf),generic_rank2_SS20_allhypotheses=True,whole_double_C_MPS_containment='Genericcanonicalsource openonly, extendtoentireoriginalirreduc6source byexplicitmagneticnormalizer/closedfiniteMPSunion afterindependentreview.')
  else:result['next_step']='No membershipfromcanonicalparticular. Completepencilremainingfreefamily needsseparatefinitejustification; no arbitraryparticularpointfarm.'
 result['seconds']=time.monotonic()-started;out=ROOT/'reports/CLASSIFICATION_NIL22_DOUBLE_C_SS20_PENCIL_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in ['linear_rank','augmented_rank','seconds']}))
if __name__=='__main__':main()

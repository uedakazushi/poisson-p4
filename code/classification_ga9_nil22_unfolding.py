#!/usr/bin/env python3
"""Complete literal cotangent-pencil test for one prescribed Ga9 unfolding.

Negative scope: N_e=y*d_x+(1+e*w)*d_z, f={y,w}=0 and P_e=P0+eR.
This is not a higher-jet calculation or a Ga9/MPS nonmembership theorem.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import minor_certificate
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();x,y,z,w=s.symbols('x y z w');xx=(x,y,z,w);base=(y,w);mons=[1,y,w,y*y,y*w,w*w];aa=s.symbols('a0:24');mag=s.symbols('b0:6')
 M=s.Matrix(2,2,[sum(aa[6*k+j]*m for j,m in enumerate(mons))for k in range(4)]);div=[sum(s.diff(M[i,j],base[j])for j in range(2))for i in range(2)]
 R=s.zeros(4);R[0,1]=M[0,0];R[0,3]=M[0,1];R[1,2]=-M[1,0];R[2,3]=M[1,1];R[0,2]=(w+div[1])*x-div[0]*z+sum(b*m for b,m in zip(mag,mons))
 for i,j in itertools.combinations(range(4),2):R[j,i]=-R[i,j]
 triples=list(itertools.combinations(range(4),3))
 def jac(P):return [s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l])for l in range(4)))for i,j,k in triples]
 JR=jac(R);assert all(s.diff(f,b)==0 for f in JR for b in mag)
 dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ)
 rank=lambda A:len(dm(A).rref()[1])
 records=[];inputs=[Path(__file__),ROOT/'code/classification_rank3_mixed_component_review.py']
 for label in ['QUADRATIC_W_POWER1','QUADRATIC_W_POWER3','LINEAR_W']:
  source_path=ROOT/('reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL211_ZERO_C_'+label+'_NORMALFORM_20261107.json');inputs.append(source_path)
  packet=json.loads(source_path.read_text())['source_magnetic_normalform'];P=s.zeros(4);den=s.Rational(packet['canonical_integer_scale'])
  for ij,f in packet['canonical_brackets'].items():
   i,j=map(int,ij.split(','));P[i-1,j-1]=s.sympify(f,locals={'h':1,'x':x,'y':y,'z':z,'w':w})/den;P[j-1,i-1]=-P[i-1,j-1]
  assert jac(P)==[0]*4
  cross=[s.expand(a-b)for a,b in zip(jac(P+R),JR)];assert all(s.diff(f,b)==0 for f in cross for b in mag)
  linear=[];labels=[]
  for tri,f in zip(triples,cross):
   for ex,co in s.Poly(f,*xx).terms():
    if co!=0:linear.append(co);labels.append(dict(type='cross_Jacobi',triple=list(tri),exponents=list(ex)))
  for ex,co in s.Poly(s.expand(-w*div[0]-M[0,1]),*base).terms():
   if co!=0:linear.append(co);labels.append(dict(type='necessary_Lie_N1_R',base_exponents=list(ex)))
  A,rhs=s.linear_eq_to_matrix(linear,aa);r,ar=rank(A),rank(A.row_join(rhs));record=dict(name=label,canonical_target_brackets=packet['canonical_brackets'],canonical_integer_scale=str(den),linear_system_shape=list(A.shape),linear_rank=r,augmented_rank=ar,
   complete24_anchor_unknowns=True,all_six_magnetic_corrections_Jacobi_independent=True,linear_equation_labels=labels,coefficient_matrix=[[str(z)for z in row]for row in A.tolist()],rhs=list(map(str,rhs)),rank_certificate=minor_certificate(A,r),augmented_rank_certificate=minor_certificate(A.row_join(rhs),ar))
  if ar>r:
   left=dm(A.T).nullspace().to_Matrix();which=next(i for i in range(left.rows)if(left*rhs)[i]!=0);row=left.row(which);assert row*A==s.zeros(1,24)
   record.update(no_actual_pencil=True,obstruction_stage='cross_Jacobi_linear',left_functional={str(i):str(a)for i,a in enumerate(row)if a},nonzero_pairing=str((row*rhs)[0]))
  else:
   sol=list(s.linsolve((A,rhs),aa))[0];free=sorted(set().union(*(f.free_symbols for f in sol)),key=str);rr=R.subs(dict(zip(aa,sol)),simultaneous=True).subs(dict.fromkeys(mag,0));equations=[]
   for tri,f in zip(triples,jac(rr)):
    for ex,co in s.Poly(f,*xx).terms():
     if co!=0:equations.append(dict(triple=list(tri),exponents=list(ex),coefficient=str(s.factor(co))))
   constant=next(q for q in equations if not s.sympify(q['coefficient']).free_symbols);assert s.Rational(constant['coefficient'])!=0
   record.update(no_actual_pencil=True,obstruction_stage='self_Jacobi_nonzero_constant',complete_linear_solution={str(z):str(s.factor(f))for z,f in zip(aa,sol)},remaining_anchor_parameters=list(map(str,free)),all_self_Jacobi_coefficients=equations,constant_obstruction=constant)
  records.append(record)
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),cases=records,
  scope='No actual affine-linear-in-e cotangent Poisson pencil at any of the three canonical Ga9 c0 orbits with the prescribed N_e=y*d_x+(1+e*w)*d_z. This does NOT exclude other modular deformations, nonlinear actualfamilies, MPS membership or componenthood. No source-membership conclusion is obtained.',
  linear_completeness='R has all24 degree<=2 coefficients in the two base anchor rows. divR=w*d_z fixes the two fiberlinear coefficients; all6 magnetic corrections are included and verified irrelevant to both cross/selfJacobi. Every PoissonR must satisfy [divR,R]=0, giving the necessary linear condition -w*divU1-U1_w=0. The complete crossJacobi system and all selfJacobi are retained.',
  next_route='Stop this prescribedcofactor unfolding after the compatible nonPoisson linear-module limits and the complete actualpencil obstruction. Do not add higher jets, larger Laurentcutoffs or arbitrary points.',seconds=time.monotonic()-started)
 out=ROOT/'reports/CLASSIFICATION_GA9_NIL22_UNFOLDING_DIAGNOSTIC_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,three_actual_pencil_obstructions=True,seconds=result['seconds'])))
if __name__=='__main__':main()

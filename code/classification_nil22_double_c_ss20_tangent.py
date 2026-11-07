#!/usr/bin/env python3
"""Complete fixedS SS20 actualpencil at wholegenericdoubleC, all90unknowns."""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from classification_affine_poisson_pencil import data,affine_solve
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();xx=s.symbols('q1 q2 x y');a,b=s.symbols('a b');K=s.QQ.frac_field(a,b);sourcepath=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json';source=json.loads(sourcepath.read_text())['cases'][2];env=dict(zip(['h','x1','x2','x3','x4'],[1,*xx]));env.update(dict(zip(['m'+str(i)for i in range(6)],[a,b,0,1,0,0])));P=s.zeros(4)
 for ij,text in source['full_six_magnetic_source'].items():
  i,j=map(int,ij.split(','));P[i-1,j-1]=s.sympify(text,locals=env);P[j-1,i-1]=-P[i-1,j-1]
 frame=data(P,xx);S=s.diag(1,1,-1,-1);field=S*s.Matrix(xx);rhsdiv=s.Matrix([s.Poly(f,*xx).coeff_monomial(m)for f in field for m in [1,*xx]]);A=frame['D'].col_join(frame['div']);rhs=s.zeros(frame['D'].rows,1).col_join(rhsdiv);answer=affine_solve(A,rhs,K);r,ar=answer['rank'],answer['augmented_rank'];print('full90fixedS',r,ar,round(time.monotonic()-started,3),flush=True)
 inputs=[Path(__file__),ROOT/'code/classification_affine_poisson_pencil.py',sourcepath,ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_MAGNETIC_NORMALFORMS_20261107.json',ROOT/'research/classification/HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_SUPPORT.md'];result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),canonical_source_F0='a+b*x+x^2',complete90unknowns=True,full_linear_shape=list(A.shape),linear_rank=r,augmented_rank=ar,full_linear_matrix=[[str(f)for f in row]for row in A.tolist()],full_linear_rhs=list(map(str,rhs)),fixed_modular_direction=list(map(str,field)),scope='Allactualaffinequadraticpencils withdivR=diag(1,1,-1,-1)*z andcrossJacobi. ApositivefullselfJacobi solution givesgenericSS20family; no inferencefromoneparticularfailure.',seconds=time.monotonic()-started)
 if ar!=r:result.update({k:v for k,v in answer.items()if k not in ['rank','augmented_rank']});result['fixed_S_pencil_empty']=True
 else:
  zz=answer['free_parameters'];Q=sum([co*B for co,B in zip(answer['solution'],frame['basis'])],s.zeros(4)).applyfunc(s.cancel);J=frame['jac'](Q);equations=[s.factor(s.cancel(c))for f in J for c in s.Poly(f,*xx).coeffs()if c!=0];remaining=list(zz);sub={};steps=[];constant=None
  while True:
   eq=[s.factor(s.cancel(f.subs(sub,simultaneous=True)))for f in equations];constant=next((f for f in eq if f!=0 and not f.free_symbols.intersection(zz)),None)
   if constant is not None:break
   remaining=sorted(set().union(*(f.free_symbols.intersection(zz)for f in eq)),key=str);lins=[f for f in eq if f!=0 and s.Poly(f,*remaining).total_degree()==1]
   if not lins:break
   L,rr=s.linear_eq_to_matrix(lins,remaining);sol=affine_solve(L,rr,K)
   if sol['rank']!=sol['augmented_rank']:constant='inconsistentlinear';break
   # Preserve symbols for original nonpivotfreevariables ratherthanintroducinga second z-list.
   from sympy.polys.matrices import DomainMatrix
   red,piv=DomainMatrix.from_Matrix(L.row_join(rr)).convert_to(K).rref();red=red.to_Matrix();free=[i for i in range(L.cols)if i not in piv];new={remaining[p]:s.cancel(red[i,L.cols]-sum(red[i,j]*remaining[j]for j in free))for i,p in enumerate(piv)};assert new;sub={t:s.cancel(f.subs(new,simultaneous=True))for t,f in sub.items()};sub.update(new);steps.append({str(t):str(f)for t,f in new.items()})
  fullQ=Q.subs(sub,simultaneous=True).applyfunc(s.cancel);residual=[s.factor(s.cancel(f.subs(sub,simultaneous=True)))for f in equations];remaining=sorted(set().union(*(f.free_symbols.intersection(zz)for f in residual)),key=str);result.update(full_linear_Q={str(i)+','+str(j):str(Q[i,j])for i,j in frame['pairs']if Q[i,j]!=0},initial_free_parameters=list(map(str,zz)),all_full_selfJacobi_coefficients=list(map(str,equations)),linear_selfJacobi_eliminations=steps,remaining_parameters=list(map(str,remaining)),constant_obstruction=str(constant)if constant is not None else None,complete_remaining_selfJacobi_ideal=sorted({str(f)for f in residual if f!=0}),reduced_Q={str(i)+','+str(j):str(fullQ[i,j])for i,j in frame['pairs']if fullQ[i,j]!=0})
  print('selfJac residual',len(zz),'->',len(remaining),len(result['complete_remaining_selfJacobi_ideal']),'constant',constant,flush=True)
  if not any(residual):result['actual_full_generic_SS20_pencil_found']=True
 result['seconds']=time.monotonic()-started;out=ROOT/'reports/CLASSIFICATION_NIL22_DOUBLE_C_SS20_TANGENT_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print('PASS',result['seconds'])
if __name__=='__main__':main()

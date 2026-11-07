#!/usr/bin/env python3
"""Exact obstruction to lifting a particular integrable 3D nil31 deformation.

This is a route diagnostic, not a nonexistence theorem for all 4D deformations.
The first-order inconsistency rules out every regular formal extension of this
fixed 3D family by x-independent quadratic Poisson derivations, even when the
prescribed modular s-component is omitted.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as S
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]

def main():
 begun=time.monotonic();y,z,u,e=S.symbols('y z u e');v=(y,z,u)
 q=y-z*u+e*(y*y-u*u)/2
 pi=S.Matrix([[0,-2*z*(z+e*u)/3,(y+z*u+e*(y*y-u*u)/2)/3],[2*z*(z+e*u)/3,0,2*z*(1+e*y)/3],[-(y+z*u+e*(y*y-u*u)/2)/3,-2*z*(1+e*y)/3,0]])
 V0=S.Matrix([z-y+2*z*u-y*y/3,z-2*y*z/3,1+y*u/3])
 mons=[S.Integer(1),y,z,u,y*y,y*z,y*u,z*z,z*u,u*u]
 cs=S.symbols('a0:30');V1=S.Matrix([sum(cs[10*i+j]*mons[j] for j in range(10)) for i in range(3)])
 def lie(P,V):
  return [S.expand(sum(V[k]*S.diff(P[i,j],v[k])-P[k,j]*S.diff(V[i],v[k])-P[i,k]*S.diff(V[j],v[k]) for k in range(3))) for i,j in itertools.combinations(range(3),2)]
 assert lie(pi.subs(e,0),V0)==[0]*3
 jac=S.expand(sum(pi[0,k]*S.diff(pi[1,2],v[k])+pi[1,k]*S.diff(pi[2,0],v[k])+pi[2,k]*S.diff(pi[0,1],v[k]) for k in range(3)))
 assert jac==0
 modular=S.Matrix([S.expand(sum(S.diff(pi[i,j],v[i]) for i in range(3))) for j in range(3)])
 assert modular==S.Matrix([z+e*u,0,1+e*y])
 equations=[];labels=[]
 for ij,f in zip(itertools.combinations(range(3),2),lie(pi,V0+e*V1)):
  for powers,coef in S.Poly(f.coeff(e,1),*v).terms():
   equations.append(coef);labels.append(dict(pair=list(ij),exponents=list(powers)))
 M,b=S.linear_eq_to_matrix(equations,cs)
 dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(S.QQ)
 cols=list(dm(M).rref()[1]);rows=list(dm(M[:,cols].T).rref()[1]);block=M[rows,cols]
 inv=dm(block).inv().to_Matrix();residual=b-M[:,cols]*inv*b[rows,:]
 i=next(i for i,a in enumerate(residual) if a)
 lam=S.zeros(1,M.rows);lam[i]=1
 selected=M[i,cols]*inv
 for j,row in enumerate(rows):lam[row]-=selected[j]
 assert lam*M==S.zeros(1,30)
 value=(lam*b)[0];assert value!=0
 A=S.Matrix([[0,1,0,0],[e*e,0,1,e],[0,0,0,0],[0,0,0,0]])
 assert A*(A*A-e*e*S.eye(4))==S.zeros(4)
 here=Path(__file__);inputs=[here,ROOT/'reports/CLASSIFICATION_NIL31_SHORT_TWO_THIRDS_CANONICAL_20261107.json']
 out=dict(task_id='E_multiple_hyperplane',status='FIXED_THREE_VARIABLE_FAMILY_HAS_NO_FIRST_ORDER_DERIVATION_LIFT',
  all_passed=True,coefficient_ring='QQ[y,z,u,e]',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},producer_bytes=here.read_text(),
  adapted_coordinates='s=x,u=w+e*x',q=str(q),pi_brackets=[str(pi[i,j]) for i,j in itertools.combinations(range(3),2)],
  whole_three_variable_Jacobi_zero=True,three_variable_modular=list(map(str,modular)),
  target_V0=list(map(str,V0)),target_lie_zero=True,
  first_order_unknown_monomials=list(map(str,mons)),first_order_unknown_count=30,
  first_order_system_shape=list(M.shape),first_order_matrix=M.tolist(),first_order_rhs=list(b),row_labels=labels,
  rank_M=len(cols),rank_augmented=len(cols)+1,
  left_inconsistency_certificate={str(j):str(a) for j,a in enumerate(lam) if a},left_M_zero=True,left_rhs=str(value),
  retained_minor_rows=rows,retained_minor_columns=cols,retained_minor_det=str(dm(block).det()),
  preliminary_exact_degree_one_system=[104,30],preliminary_exact_degree_two_system=[167,60],
  no_divergence_condition_used=True,
  rigorous_conclusion='No regular formal family V(e)=V0+e*V1+... of quadratic x-independent Poisson derivations can extend this fixed pi(e), because already the full first-order Lie equation is inconsistent.',
  limitation='This does not exclude x-dependent 4D Poisson deformations, a different 3D pi(e), a different epsilon-dependent translation direction, or membership in any known component. The semisimple modular FIELD candidate alone proves no E membership.',
  next_route='Do not add higher jets to this fixed candidate. Change the actual geometric construction or target a different finite source.',seconds=time.monotonic()-begun)
 target=ROOT/'reports/CLASSIFICATION_NIL31_SHORT_SEMISIMPLE_CANDIDATE_20261107.json'
 if target.exists():raise FileExistsError(str(target))
 target.write_text(json.dumps(out,indent=2,default=str)+'\n')
 print(json.dumps(dict(status=out['status'],shape=list(M.shape),rank=len(cols),left_rhs=str(value),seconds=out['seconds'])))
if __name__=='__main__':main()

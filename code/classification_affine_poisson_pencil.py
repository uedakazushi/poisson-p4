"""Reusable exact full90 affinequadratic Poisson-pencil linear algebra."""
import itertools
import sympy as s
from sympy.polys.matrices import DomainMatrix

def data(P,coords):
 pairs=list(itertools.combinations(range(4),2));mons=[s.Integer(1),*coords,*[coords[i]*coords[j]for i,j in itertools.combinations_with_replacement(range(4),2)]];triples=list(itertools.combinations(range(4),3));mon3=[s.prod(coords[i]for i in factors)for d in range(4)for factors in itertools.combinations_with_replacement(range(4),d)]
 def jac(Q):return [s.expand(sum(Q[i,l]*s.diff(Q[j,k],coords[l])+Q[j,l]*s.diff(Q[k,i],coords[l])+Q[k,l]*s.diff(Q[i,j],coords[l])for l in range(4)))for i,j,k in triples]
 def modular(Q):return s.Matrix([sum(s.diff(Q[i,j],coords[i])for i in range(4))for j in range(4)])
 assert jac(P)==[0]*4;basis=[];columns=[];divcolumns=[]
 for ij,mon in itertools.product(pairs,mons):
  Q=s.zeros(4);Q[ij]=mon;Q[ij[1],ij[0]]=-mon;basis.append(Q);columns.append(s.Matrix([s.Poly(f,*coords).coeff_monomial(m)for f in jac(P+Q)for m in mon3]));divcolumns.append(s.Matrix([s.Poly(f,*coords).coeff_monomial(m)for f in modular(Q)for m in [1,*coords]]))
 return dict(pairs=pairs,mons=mons,mon3=mon3,jac=jac,modular=modular,basis=basis,D=s.Matrix.hstack(*columns),div=s.Matrix.hstack(*divcolumns))

def affine_solve(A,rhs,K):
 dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(K);reduced,pivots=dm(A.row_join(rhs)).rref();r=len(dm(A).rref()[1]);ar=len(pivots)
 if r!=ar:
  left=dm(A.T).nullspace(divide_last=True).to_Matrix();row=next(left.row(i)for i in range(left.rows)if s.cancel((left.row(i)*rhs)[0])!=0);assert(row*A).applyfunc(s.cancel)==s.zeros(1,A.cols)
  return dict(rank=r,augmented_rank=ar,left_functional={str(i):str(a)for i,a in enumerate(row)if a},nonzero_pairing=str(s.cancel((row*rhs)[0])))
 free=[i for i in range(A.cols)if i not in pivots];zz=s.symbols('z0:'+str(len(free)));reduced=reduced.to_Matrix();solution=s.zeros(A.cols,1)
 for i,z in zip(free,zz):solution[i]=z
 for i,pivot in enumerate(pivots):solution[pivot]=s.cancel(reduced[i,A.cols]-sum(reduced[i,j]*solution[j]for j in free))
 assert(A*solution-rhs).applyfunc(s.cancel)==s.zeros(A.rows,1)
 return dict(rank=r,augmented_rank=ar,solution=solution,free_parameters=zz)

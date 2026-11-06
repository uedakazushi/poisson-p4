#!/usr/bin/env python3
"""The complete normalized fiber-cubic subcase of the rank1-isotropic source."""
import hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
start=time.monotonic();path=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_AZERO_CZERO_RANK1_ISOTROPIC_20261107.json');src=json.loads(path.read_text());xx=s.symbols('x0:5');ns=s.symbols('n0:22');P=s.zeros(5)
fixed={ns[7]:-2,ns[18]:-3,ns[11]:0,ns[6]:0,ns[15]:0};pars=[n for n in ns if n not in fixed]
for ij in itertools.combinations(range(5),2):P[ij]=s.expand(s.sympify(src['brackets'][''.join(map(str,ij))]).subs(fixed));P[ij[::-1]]=-P[ij]
a,b,c,r,z=xx;q=2*a*c+b*b;C=a*a*z+b**3+3*a*b*c;X=s.Matrix([0,0,0,a,0]);equations=[]
for j in range(5):equations.extend(s.Poly(s.expand(sum(s.diff(C,xx[i])*P[i,j] for i in range(5))+C*X[j]),*xx).coeffs())
M,rhs=s.linear_eq_to_matrix(equations,pars);dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ)
rank=len(dm(M).rref()[1]);aug=len(dm(M.row_join(rhs)).rref()[1]);result={'task_id':'E_reduced_nonnormal','source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'canonical_fixed_parameters':{str(k):str(v) for k,v in fixed.items()},'canonical_cubic':str(C),'linear_casimir_matrix_shape':list(M.shape),'linear_rank_QQ':rank,'augmented_rank_QQ':aug,'inconsistent':aug>rank}
if aug==rank:
 sol=next(iter(s.linsolve((M,rhs),pars)));sub=dict(zip(pars,sol));new=P.subs(sub,simultaneous=True).applyfunc(s.expand);free=sorted(set().union(*(f.free_symbols for f in sol))&set(pars),key=str)
 jac={ijk:s.expand(sum(new[i,l]*s.diff(new[j,k],xx[l])+new[j,l]*s.diff(new[k,i],xx[l])+new[k,l]*s.diff(new[i,j],xx[l]) for l in range(5))) for ijk in itertools.combinations(range(5),3) for i,j,k in [ijk]}
 coeff=[]
 for f in jac.values():
  for cc in s.Poly(f,*xx).coeffs():
   cc=s.factor(cc)
   if cc and cc not in coeff:coeff.append(cc)
 pf=[]
 for missing in range(5):
  i,j,k,l=[v for v in range(5) if v!=missing];pf.append((-1)**missing*(new[i,j]*new[k,l]-new[i,k]*new[j,l]+new[i,l]*new[j,k]))
 D=s.factor(sum(x*f for x,f in zip(xx,pf))/5)
 result.update(full_linear_substitution={str(k):str(v) for k,v in sub.items()},free_parameters=list(map(str,free)),brackets={''.join(map(str,ij)):str(new[ij]) for ij in itertools.combinations(range(5),2)},remaining_Jacobi_coefficients=list(map(str,coeff)),Euler_contracted_Pfaffian=str(D),divisor_ratio=s.cancel(D/(q*C)).__str__())
result.update(producer_source=Path(__file__).read_text(),seconds=time.monotonic()-start)
out=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_AZERO_FIBER_20261107.json')
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['producer_source','brackets']},indent=2))

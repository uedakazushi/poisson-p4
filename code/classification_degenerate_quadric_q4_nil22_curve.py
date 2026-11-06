#!/usr/bin/env python3
"""Actual one-parameter curve from the finite full equivariant I-slice."""
import hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
start=time.monotonic();path=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL22_WEIGHTED_20261107.json');src=json.loads(path.read_text());xx=s.symbols('h x1 x2 x3 x4');tau=s.symbols('tau0:3');t=s.Symbol('t');P=s.zeros(5)
for key,val in src['full_curve_brackets'].items():
    ij=tuple(map(int,key.split(',')));P[ij]=s.expand(s.sympify(val).subs({tau[0]:0,tau[1]:t,tau[2]:0}));P[ij[::-1]]=-P[ij]
J={ijk:s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l]) for l in range(5))) for ijk in itertools.combinations(range(5),3) for i,j,k in [ijk]}
full={ijkl:s.expand(sum((-1)**at*xx[ell]*J[ijkl[:at]+ijkl[at+1:]] for at,ell in enumerate(ijkl))) for ijkl in itertools.combinations(range(5),4)};assert all(v==0 for v in full.values())
pf=[]
for missing in range(5):
    i,j,k,l=[v for v in range(5) if v!=missing];pf.append(s.expand((-1)**missing*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])))
D=s.factor(sum(x*f for x,f in zip(xx,pf))/5);factors=s.factor_list(D)[1]
quads=[f for f,e in factors if s.Poly(f,*xx).total_degree()==2]
q=quads[0] if len(quads)==1 else None
record={'task_id':'E_reduced_nonnormal','all_parameter_projective_Jacobi_zero':True,'family':{','.join(map(str,ij)):str(P[ij]) for ij in itertools.combinations(range(5),2) if P[ij]},'Euler_contracted_Pfaffian':str(D),'factor_degrees_multiplicities':[(s.Poly(f,*xx).total_degree(),e) for f,e in factors],'quadratic_factor':None if q is None else str(q),'quadratic_hessian_determinant':None if q is None else str(s.factor(s.hessian(q,xx).det())),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-start}
out=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL22_CURVE_20261107.json')
with out.open('x') as f:json.dump(record,f,indent=2);f.write('\n')
print(json.dumps({k:record[k] for k in ['Euler_contracted_Pfaffian','factor_degrees_multiplicities','quadratic_factor','quadratic_hessian_determinant','seconds']},indent=2))

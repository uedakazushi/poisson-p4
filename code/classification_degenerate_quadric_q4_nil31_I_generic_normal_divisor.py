#!/usr/bin/env python3
"""One finite Macaulay certificate that the actual quintic is normal."""
import hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
start=time.monotonic();p=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL31_I_WEIGHTED_20261107.json');src=json.loads(p.read_text());h,b,c,d,z=s.symbols('h x1 x2 x3 x4');tau=s.symbols('tau0:3');xx=(h,b,c,d,z)
P=s.zeros(5)
for key,val in src['full_curve_brackets'].items():
    ij=tuple(map(int,key.split(',')));P[ij]=s.expand(s.sympify(val).subs(dict(zip(tau,[1,2,3]))));P[ij[::-1]]=-P[ij]
pf=[]
for missing in range(5):
    i,j,k,l=[v for v in range(5) if v!=missing];pf.append((-1)**missing*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k]))
D=s.expand(sum(x*f for x,f in zip(xx,pf)));subs={h:-2*z-3*b,c:-5*b-7*d-11*z};yy=(b,d,z)
polys=[s.expand(s.diff(D,v).subs(subs,simultaneous=True)) for v in xx]
def mons(degree):return [(i,j,degree-i-j) for i in range(degree,-1,-1) for j in range(degree-i,-1,-1)]
N=8;rows=mons(N);multipliers=mons(N-4);columns=[(i,e) for i in range(5) for e in multipliers]
M=s.zeros(len(rows),len(columns));rowindex={e:i for i,e in enumerate(rows)}
for j,(at,e) in enumerate(columns):
    for ex,a in s.Poly(polys[at],*yy).terms():M[rowindex[tuple(x+y for x,y in zip(ex,e))],j]=a
mat=DomainMatrix.from_Matrix(M).convert_to(s.QQ);rref,pivots=mat.rref();rank=len(pivots);det=None
if rank==len(rows):det=str(mat.extract(range(len(rows)),pivots).det())
report={'task_id':'E_reduced_nonnormal','source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'witness_parameter':'(tau0,tau1,tau2)=(1,2,3)','actual_quintic':str(D),'linear_sections':['h+3*x1+2*x4','x2+5*x1+7*x3+11*x4'],'restricted_five_derivatives':list(map(str,polys)),'Macaulay_total_degree':N,'row_monomial_exponents':rows,'column_derivative_and_multiplier_exponents':columns,'matrix_shape':list(M.shape),'rank_QQ':rank,'full_row_rank':rank==len(rows),'minor_columns':list(pivots),'minor_determinant':det,'projective_singular_locus_dimension_at_most_one':rank==len(rows),'implication':'If full rank: the chosen codimension-two linear section has empty projective intersection with Sing(D). Projective dimension theorem forces dim Sing(D)<=1; the quintic hypersurface in P4 is reduced and normal by the hypersurface Serre criterion.','producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-start}
out=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL31_I_GENERIC_NORMAL_DIVISOR_20261107.json')
with out.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
print(json.dumps({k:report[k] for k in ['matrix_shape','rank_QQ','full_row_rank','seconds']},indent=2))

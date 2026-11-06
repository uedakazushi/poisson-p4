#!/usr/bin/env python3
"""One exact genuine transverse quadratic test for all six complete source orbits."""
import hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_three_hyperplanes_quadric_independent import ambient_matrices
from classification_rank3_mixed_component_review import COLUMNS,jacobi,sparse
resource.setrlimit(resource.RLIMIT_CPU,(85,86));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
start=time.monotonic();path=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_CANONICAL_20261107.json');src=json.loads(path.read_text());bounds=json.loads(Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_COMPONENTS_20261107.json').read_text());xx=s.symbols('a b c r z');dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ);records=[]
def enc(P):return {str(ij):str(s.expand(P[ij])) for ij in itertools.combinations(range(5),2) if P[ij]}
def tensor(v):
 R=s.zeros(5)
 for ((i,j),ex),a in zip(COLUMNS,v):R[i,j]+=a*s.prod(x**e for x,e in zip(xx,ex))
 return R-R.T
for case in src['cases']:
 begun=time.monotonic();P=s.zeros(5)
 for ij in itertools.combinations(range(5),2):P[ij]=s.sympify(case['canonical_point_brackets'].get(str(ij),'0'));P[ij[::-1]]=-P[ij]
 D,V,G,F=ambient_matrices(P,xx,[s.zeros(5)]);rref,pivots=dm(D).rref();rref=rref.to_Matrix();drank=len(pivots);grank=len(dm(G).rref()[1]);expected=next(c for c in bounds['cases'] if c['case']==case['case']);assert drank==expected['ranks_QQ']['full_D'] and grank==48
 K=[]
 for j in range(150):
  if j in pivots:continue
  v=s.zeros(150,1);v[j]=1
  for i,k in enumerate(pivots):v[k]=-rref[i,j]
  v*=s.ilcm(1,*[s.denom(a) for a in v]);K.append(v)
 K=s.Matrix.hstack(*K);chosen=[j-G.cols for j in dm(G.row_join(K)).rref()[1] if j>=G.cols];assert len(chosen)==150-drank-grank
 normals=[tensor(K[:,j]) for j in chosen];assert D*K[:,chosen]==s.zeros(350,len(chosen));taus=s.symbols('n0:'+str(len(normals)));monpairs=list(itertools.combinations_with_replacement(range(len(normals)),2));Q=[];selfQ=[jacobi(sparse(R,xx)) for R in normals]
 for i,j in monpairs:Q.append(selfQ[i] if i==j else jacobi(sparse(normals[i]+normals[j],xx))-selfQ[i]-selfQ[j])
 Q=s.Matrix.hstack(*Q);ind=D[:,list(pivots)];rows=list(dm(ind.T).rref()[1]);block=dm(D.extract(rows,list(pivots)));inverse=block.inv();coeff=inverse.matmul(dm(Q[rows,:])).to_Matrix();res=Q-ind*coeff;witness_rows=list(dm(res.T).rref()[1]);polys=[s.factor(sum(res[i,j]*taus[a]*taus[b] for j,(a,b) in enumerate(monpairs))) for i in witness_rows];rpolys=dm(res.extract(witness_rows,range(res.cols))).rref()[0].to_Matrix();basis=[s.factor(sum(rpolys[i,j]*taus[a]*taus[b] for j,(a,b) in enumerate(monpairs))) for i in range(rpolys.rows) if any(rpolys[i,j] for j in range(rpolys.cols))]
 funcs=[]
 for row in witness_rows:
  lam=s.zeros(1,350);lam[row]=1;weights=dm(ind[row,:]).matmul(inverse).to_Matrix()
  for j,k in enumerate(rows):lam[k]-=weights[j]
  assert lam*D==s.zeros(1,150);funcs.append({'row':row,'functional':{str(i):str(a) for i,a in enumerate(lam) if a},'pairings':list(map(str,lam*Q))})
 if not basis:zero=False;supportdim=len(normals)
 elif len(normals)==1:zero=True;supportdim=0
 else:
  gcd=basis[0]
  for f in basis[1:]:gcd=s.gcd(gcd,f)
  zero=s.Poly(gcd,*taus).total_degree()==0;supportdim=0 if zero else 1
 rec={'case':case['case'],'full_D_rank_QQ':drank,'actual_group_rank_QQ':grank,'full_kernel_dimension':150-drank,'true_normal_count':len(normals),'literal_normals':[enc(R) for R in normals],'normal_variable_order':list(map(str,taus)),'quadratic_monomial_pairs':monpairs,'full_D_times_normals_zero':True,'left_null_witnesses':funcs,'true_quadratic_rowspace_basis':list(map(str,basis)),'quadratic_initial_support_dimension':supportdim,'quadratic_initial_support_is_origin':zero,'component23_if_support_is_origin':zero,'all_normal_linear_combinations_are_actual_Poisson':Q==s.zeros(350,Q.cols),'implicit_D_minor_rows':rows,'implicit_D_minor_columns':list(pivots),'implicit_D_minor_determinant':str(block.det()),'seconds':time.monotonic()-begun}
 records.append(rec);print(json.dumps({k:rec[k] for k in ['case','true_normal_count','true_quadratic_rowspace_basis','quadratic_initial_support_dimension','all_normal_linear_combinations_are_actual_Poisson','seconds']}),flush=True)
result={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-start,'scope':'One true-normal quadratic computation retaining the full150 derivative and actual orbit directions for each complete source. Positive-dimensional initial support is not equated with an actual deformation.'}
out=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_NORMALS_20261107.json')
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')

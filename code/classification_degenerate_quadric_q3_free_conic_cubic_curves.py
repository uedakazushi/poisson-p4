#!/usr/bin/env python3
"""Check one specified quadratic actual-curve candidate per unresolved normal line."""
import hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_three_hyperplanes_quadric_independent import ambient_matrices
from classification_rank3_mixed_component_review import COLUMNS,jacobi,sparse
resource.setrlimit(resource.RLIMIT_CPU,(85,86));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
start=time.monotonic();sourcepath=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_CANONICAL_20261107.json');normalpath=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_NORMALS_20261107.json');sources=json.loads(sourcepath.read_text());normaldata=json.loads(normalpath.read_text());xx=s.symbols('a b c r z');t=s.Symbol('t');dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ);records=[]
def dec(data):
 P=s.zeros(5)
 for ij in itertools.combinations(range(5),2):P[ij]=s.sympify(data.get(str(ij),'0'));P[ij[::-1]]=-P[ij]
 return P
def enc(P):return {str(ij):str(s.expand(P[ij])) for ij in itertools.combinations(range(5),2) if P[ij]}
def tensor(v):
 P=s.zeros(5)
 for ((i,j),e),a in zip(COLUMNS,v):P[i,j]+=a*s.prod(x**n for x,n in zip(xx,e))
 return P-P.T
for name,index in [('cusp_outside_3',1),('cusp_at_cusp_2',0),('cusp_at_cusp_3',0)]:
 begun=time.monotonic();case=next(c for c in sources['cases'] if c['case']==name);normal=next(c for c in normaldata['cases'] if c['case']==name);P=dec(case['canonical_point_brackets']);R=dec(normal['literal_normals'][index]);D,*_=ambient_matrices(P,xx,[s.zeros(5)]);rows=normal['implicit_D_minor_rows'];cols=normal['implicit_D_minor_columns'];Q=jacobi(sparse(R,xx));block=dm(D.extract(rows,cols));coeff=block.inv().matmul(dm(-Q[rows,:])).to_Matrix();v=s.zeros(150,1)
 for i,j in enumerate(cols):v[j]=coeff[i]
 assert D*v+Q==s.zeros(350,1);T=tensor(v);curve=P+t*R+t*t*T
 raw={ijk:s.expand(sum(curve[i,l]*s.diff(curve[j,k],xx[l])+curve[j,l]*s.diff(curve[k,i],xx[l])+curve[k,l]*s.diff(curve[i,j],xx[l]) for l in range(5))) for ijk in itertools.combinations(range(5),3) for i,j,k in [ijk]}
 full={ijkl:s.expand(sum((-1)**at*xx[ell]*raw[ijkl[:at]+ijkl[at+1:]] for at,ell in enumerate(ijkl))) for ijkl in itertools.combinations(range(5),4)};good=all(v==0 for v in full.values());rec={'case':name,'normal_index':index,'second_order_correction':enc(T),'exact_D_correction_plus_quadratic_zero':True,'specified_quadratic_curve_is_actual_projective_Poisson':good,'curve':enc(curve)}
 if good:
  pf=[]
  for missing in range(5):
   i,j,k,l=[v for v in range(5) if v!=missing];pf.append((-1)**missing*(curve[i,j]*curve[k,l]-curve[i,k]*curve[j,l]+curve[i,l]*curve[j,k]))
  potential=s.factor(sum(x*f for x,f in zip(xx,pf))/5);fac=s.factor_list(potential)[1];rec.update(actual_Euler_Pfaffian=str(potential),factor_degrees=[(s.Poly(f,*xx).total_degree(),int(e)) for f,e in fac],quadric_factors=[{'factor':str(f),'Hessian_rank_generic':s.hessian(f,xx).rank(),'Hessian_determinant':str(s.factor(s.hessian(f,xx).det()))} for f,e in fac if s.Poly(f,*xx).total_degree()==2])
 else:
  rec['remaining_projective_Jacobi']={str(ijkl):str(s.factor(f)) for ijkl,f in full.items() if f};rec['stop']='This specified finite curve fails. No higher jet or actual deformation is inferred.'
 rec['seconds']=time.monotonic()-begun;records.append(rec);print(json.dumps({k:rec[k] for k in ['case','specified_quadratic_curve_is_actual_projective_Poisson','seconds']}),flush=True)
result={'task_id':'E_reduced_nonnormal','specified_finite_candidates_checked':True,'cases':records,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [sourcepath,normalpath]},'producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-start}
out=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_CURVES_20261107.json')
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')

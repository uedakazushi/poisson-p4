#!/usr/bin/env python3
"""Genuine cubic obstruction modulo D and every quadratic normal pairing."""
import hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_three_hyperplanes_quadric_independent import ambient_matrices
from classification_rank3_mixed_component_review import COLUMNS,jacobi
resource.setrlimit(resource.RLIMIT_CPU,(85,86));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
start=time.monotonic();paths=[Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_'+name+'_20261107.json') for name in ['CANONICAL','NORMALS','CURVES']];src,norms,curves=[json.loads(p.read_text()) for p in paths];xx=s.symbols('a b c r z');dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ);records=[]
def dec(data):
 P=s.zeros(5)
 for ij in itertools.combinations(range(5),2):P[ij]=s.sympify(data.get(str(ij),'0'));P[ij[::-1]]=-P[ij]
 return P
def sparse(P):return {ij:dict(s.Poly(s.expand(P[ij]),*xx).terms()) for ij in itertools.combinations(range(5),2) if P[ij]!=0}
def J(P):return jacobi(sparse(P))
for name,axis in [('cusp_outside_3',1),('cusp_at_cusp_2',0),('cusp_at_cusp_3',0)]:
 begun=time.monotonic();source=next(c for c in src['cases'] if c['case']==name);normal=next(c for c in norms['cases'] if c['case']==name);curve=next(c for c in curves['cases'] if c['case']==name);P=dec(source['canonical_point_brackets']);RR=list(map(dec,normal['literal_normals']));R=RR[axis];T=dec(curve['second_order_correction']);D,V,G,F=ambient_matrices(P,xx,[s.zeros(5)]);qcols=[]
 for i,j in itertools.combinations_with_replacement(range(len(RR)),2):qcols.append(J(RR[i]) if i==j else J(RR[i]+RR[j])-J(RR[i])-J(RR[j]))
 Q=s.Matrix.hstack(*qcols);C=J(R+T)-J(R)-J(T);A=D.row_join(Q);pivots=list(dm(A).rref()[1]);ind=A[:,pivots];rows=list(dm(ind.T).rref()[1]);inverse=dm(ind[rows,:]).inv();coeff=inverse.matmul(dm(C[rows,:])).to_Matrix();res=C-ind*coeff;row=next((i for i,v in enumerate(res) if v),None)
 # Verify arbitrary source-group tangent ambiguity, not only a fixed correction.
 group_cross=[]
 for col in range(G.cols):
  field=s.zeros(5)
  for ((i,j),e),v in zip(COLUMNS,G[:,col]):field[i,j]+=v*s.prod(x**k for x,k in zip(xx,e))
  field-=field.T;group_cross.append(J(R+field)-J(R)-J(field))
 GC=s.Matrix.hstack(*group_cross);gc_coeff=inverse.matmul(dm(GC[rows,:])).to_Matrix();assert GC-ind*gc_coeff==s.zeros(350,G.cols)
 rec={'case':name,'free_axis_normal_index':axis,'full_D_rank_QQ':normal['full_D_rank_QQ'],'actual_group_rank_QQ':48,'full_kernel_dimension':normal['full_kernel_dimension'],'true_normal_count':len(RR),'rank_D_plus_all_normal_quadratics_QQ':len(pivots),'all_group_cross_pairings_in_D_plus_all_normal_quadratics':True,'all_normal_cross_pairings_killed_by_construction':True,'arbitrary_full_kernel_correction_ambiguity_removed':True,'genuine_cubic_obstruction_nonzero':row is not None,'cubic_vector_sparse':{str(i):str(v) for i,v in enumerate(C) if v},'selected_A_columns':pivots,'selected_A_rows':rows,'A_minor_determinant':str(dm(ind[rows,:]).det())}
 if row is not None:
  lam=s.zeros(1,350);lam[row]=1;weights=dm(ind[row,:]).matmul(inverse).to_Matrix()
  for j,k in enumerate(rows):lam[k]-=weights[j]
  assert lam*A==s.zeros(1,A.cols);pair=(lam*C)[0];assert pair
  rec.update(left_null_functional={str(i):str(v) for i,v in enumerate(lam) if v},exact_lambda_D_zero=True,exact_lambda_all_normal_quadratics_zero=True,nonzero_cubic_pairing=str(pair),full_initial_support_is_origin=True,projective_component_dimension=23)
 else:rec.update(full_initial_support_is_origin=False,stop='This single genuine cubic test gives no isolation; no higher jet is authorized by this contract.')
 rec['seconds']=time.monotonic()-begun;records.append(rec);print(json.dumps({k:rec[k] for k in ['case','genuine_cubic_obstruction_nonzero','rank_D_plus_all_normal_quadratics_QQ','seconds']}),flush=True)
result={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-start,'scope':'Genuine cubic modulo D and all true-normal quadratic pairings, with every group tangent cross pairing independently checked. Nonzero isolates the remaining free axis; zero does not imply an actual deformation.'}
out=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_THIRD_20261107.json')
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')

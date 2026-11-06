#!/usr/bin/env python3
"""The single true transverse normal at the nil31-long-II orbit."""
import argparse,hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_three_hyperplanes_quadric_independent import ambient_matrices
from classification_rank3_mixed_component_review import COLUMNS,jacobi,sparse,vector

ROOT=Path(__file__).resolve().parents[1]
SOURCE='reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_COMPONENTS_20261107.json'
def dm(A):return DomainMatrix.from_Matrix(A).convert_to(s.QQ)
def decode(data):
    A=s.zeros(5)
    for ij in itertools.combinations(range(5),2):A[ij]=s.sympify(data.get(str(ij),'0'));A[ij[::-1]]=-A[ij]
    return A
def encode(A):return {str(ij):str(s.expand(A[ij])) for ij in itertools.combinations(range(5),2) if A[ij]}
def digest(A):return hashlib.sha256(json.dumps(A.tolist(),default=str).encode()).hexdigest()

def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();source=json.loads((ROOT/SOURCE).read_text());saved=next(r for r in source['cases'] if r['source']=='nil31_long_II')
    xx=s.symbols('x0:5');P=decode(saved['point_brackets']);directions=[decode(d) for d in saved['all_actual_parameter_directions']]
    D,V,M,F=ambient_matrices(P,xx,directions)
    for key,A in [('full350_D',D),('vertical',V),('orbit_vertical_scale',M),('actual_full_source',F)]:assert digest(A)==saved['matrix_sha256'][key]
    rr,pivots=dm(D).rref();rr=rr.to_Matrix();assert len(pivots)==102
    kernel=[]
    for j in range(150):
        if j in pivots:continue
        v=s.zeros(150,1);v[j]=1
        for i,k in enumerate(pivots):v[k]=-rr[i,j]
        v*=s.ilcm(1,*[s.denom(c) for c in v]);kernel.append(v)
    K=s.Matrix.hstack(*kernel);assert D*K==s.zeros(350,48)
    selected=dm(F.row_join(K)).rref()[1]
    normal_indices=[i-F.cols for i in selected if i>=F.cols]
    assert len([i for i in selected if i<F.cols])==47 and len(normal_indices)==1
    nv=K[:,normal_indices[0]];R=s.zeros(5)
    for (ij,exponents),a in zip(COLUMNS,nv):R[ij]+=a*s.prod(x**e for x,e in zip(xx,exponents))
    for ij in itertools.combinations(range(5),2):R[ij[::-1]]=-R[ij]
    assert D*nv==s.zeros(350,1)
    Q=jacobi(sparse(R,xx));independent=D[:,list(pivots)]
    rows=list(dm(independent.T).rref()[1]);block=dm(D.extract(rows,list(pivots)));inverse=block.inv()
    coefficients=inverse.matmul(dm(Q.extract(rows,[0]))).to_Matrix()
    obstruction=Q-independent*coefficients
    assert obstruction.extract(rows,[0])==s.zeros(102,1)
    nonzero=next((i for i,a in enumerate(obstruction) if a),None)
    result={'task_id':'E_reduced_nonnormal','all_passed':True,'source':'nil31_long_II',
        'full_D_rank_QQ':102,'actual_family_rank_QQ':47,'full_kernel_dimension_QQ':48,
        'true_normal_count':1,'literal_normal':encode(R),'normal_kernel_index':normal_indices[0],
        'D_times_normal_zero':True,'actual_F_plus_normal_rank':48,
        'quadratic_Jacobi_vector_sparse':{str(i):str(a) for i,a in enumerate(Q) if a},
        'eliminated_D_rows':rows,'eliminated_D_columns':list(pivots),'D_minor_determinant':str(block.det()),
        'quadratic_obstruction_nonzero':nonzero is not None,'normal_itself_is_projectively_Poisson':Q==s.zeros(350,1),
        'quadratic_support_zero_dimensional':nonzero is not None,
        'projective_component_dimension_if_support_is_zero_dimensional':22 if nonzero is not None else None}
    if nonzero is not None:
        weights=dm(independent[nonzero,:]).matmul(inverse).to_Matrix()
        functional=s.zeros(1,350);functional[nonzero]=1
        for j,row in enumerate(rows):functional[row]-=weights[j]
        assert functional*D==s.zeros(1,150)
        pairing=(functional*Q)[0];assert pairing==obstruction[nonzero] and pairing
        result.update(left_null_functional={str(i):str(v) for i,v in enumerate(functional) if v},
            exact_functional_times_all_D_columns_zero=True,quadratic_pairing=str(pairing),
            true_quadratic_initial_support='(normal)=0')
    else:
        correction=s.zeros(150,1)
        for i,j in enumerate(pivots):correction[j]=-coefficients[i]
        assert D*correction+Q==s.zeros(350,1)
        result.update(second_order_correction_sparse={str(i):str(v) for i,v in enumerate(correction) if v},
            exact_D_correction_plus_quadratic_zero=True,
            implication='Quadratic obstruction does not isolate the orbit; no componenthood or actual deformation is inferred from this second-order correction.')
    paths=[SOURCE,'code/classification_three_hyperplanes_quadric_independent.py','code/classification_rank3_mixed_component_review.py']
    result.update(input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun,global_exhaustion_proved=False)
    with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','quadratic_obstruction_nonzero','normal_itself_is_projectively_Poisson','quadratic_support_zero_dimensional','seconds']},indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output)

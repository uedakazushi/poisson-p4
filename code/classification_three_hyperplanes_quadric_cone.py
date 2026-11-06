#!/usr/bin/env python3
"""Test the complete three-hyperplane logarithmic ansatz modulo true source.

This is distinct from the cone-algebroid ansatz: Pxy may be nonzero, mixed
entries may involve either fiber coordinate, and Pzw is unrestricted quadratic.
"""
import argparse,hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,jacobi,sparse,vector
from classification_three_hyperplanes_quadric_independent import ambient_matrices

ROOT=Path(__file__).resolve().parents[1]
SOURCE='reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_20261107.json'
def dm(A):return DomainMatrix.from_Matrix(A).convert_to(s.QQ)
def decode(data):
    A=s.zeros(5)
    for ij in itertools.combinations(range(5),2):
        A[ij]=s.sympify(data.get(str(ij),'0'));A[ij[::-1]]=-A[ij]
    return A
def from_vector(v,x):
    A=s.zeros(5)
    for (ij,e),value in zip(COLUMNS,v):A[ij]+=value*s.prod(xx**power for xx,power in zip(x,e))
    for ij in itertools.combinations(range(5),2):A[ij[::-1]]=-A[ij]
    return A
def encode(A):return {str(ij):str(A[ij])for ij in itertools.combinations(range(5),2)if A[ij]}
def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();h,x,y,z,w=xx=s.symbols('h x y z w')
    record=json.loads((ROOT/SOURCE).read_text())['cases'][1];assert record['case']=='m1_t2'
    P=decode(record['point_brackets']);D,V,M,F=ambient_matrices(P,xx,[decode(a)for a in record['all6_actual_parameter_directions']])
    assert hashlib.sha256(json.dumps(D.tolist(),default=str).encode()).hexdigest()==record['matrix_sha256']['full350_D']
    base2=[h*h,h*x,h*y,x*x,x*y,y*y];basis=[]
    A=s.zeros(5);A[1,2]=x*y;A[2,1]=-x*y;basis.append(A)
    for ij,factor in [((1,3),x),((1,4),x),((2,3),y),((2,4),y)]:
        for variable in xx:
            A=s.zeros(5);A[ij]=factor*variable;A[ij[::-1]]=-A[ij];basis.append(A)
    for i,j in itertools.combinations_with_replacement(range(5),2):
        mon=xx[i]*xx[j]
        A=s.zeros(5);A[3,4]=mon;A[4,3]=-mon;basis.append(A)
    C=s.Matrix.hstack(*[vector(sparse(A,xx))for A in basis]);assert C.cols==36
    rr,pivots=dm(D*C).rref();rr=rr.to_Matrix();kernel=[]
    for j in range(36):
        if j in pivots:continue
        v=s.zeros(36,1);v[j]=1
        for i,p in enumerate(pivots):v[p]=-rr[i,j]
        v*=s.ilcm(1,*[s.denom(a)for a in v]);kernel.append(v)
    K=C*s.Matrix.hstack(*kernel);assert D*K==s.zeros(350,K.cols)
    combined=list(dm(F.row_join(K)).rref()[1]);indices=[j-F.cols for j in combined if j>=F.cols]
    assert len(indices)<=1;directions=[]
    for index in indices:
        A=from_vector(K[:,index],xx);J=jacobi(sparse(A,xx))
        directions.append(dict(brackets=encode(A),projective_self_Jacobi_zero=J==s.zeros(350,1),
                               projective_Jacobi_nonzero_coefficients=sum(bool(v)for v in J)))
    inputs=[SOURCE,'code/classification_three_hyperplanes_quadric_cone.py',
            'code/classification_three_hyperplanes_quadric_independent.py','code/classification_rank3_mixed_component_review.py']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,case='m1_t2',
        complete_three_hyperplane_ansatz_dimension=36,logarithmic_linearized_rank=len(pivots),logarithmic_tangent_kernel_dimension=K.cols,
        actual_source_rank=48,source_plus_cone_tangent_rank=len(combined),
        genuine_normal_has_three_hyperplane_representative=bool(indices),logarithmic_normal_representatives=directions,
        independent_matrix_hashes=dict(cone=hashlib.sha256(json.dumps(C.tolist(),default=str).encode()).hexdigest(),
                                       cone_jacobian=hashlib.sha256(json.dumps((D*C).tolist(),default=str).encode()).hexdigest()),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in inputs},
        producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps(result|{'producer_source':'[stored in report]'}))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)

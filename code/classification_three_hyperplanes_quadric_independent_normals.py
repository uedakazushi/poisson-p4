#!/usr/bin/env python3
"""One true-normal obstruction computation for each nonsmooth conic source."""
import argparse, hashlib, itertools, json, resource, time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,jacobi,sparse
from classification_three_hyperplanes_quadric_independent import ambient_matrices

ROOT=Path(__file__).resolve().parents[1]
SOURCE='reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_20261107.json'

def dm(A):return DomainMatrix.from_Matrix(A).convert_to(s.QQ)

def decode(encoded):
    A=s.zeros(5)
    for ij in itertools.combinations(range(5),2):
        A[ij]=s.sympify(encoded.get(str(ij),'0'));A[ij[::-1]]=-A[ij]
    return A

def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(90,91))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();saved=json.loads((ROOT/SOURCE).read_text())
    for snapshot in saved['source_snapshots']:
        assert hashlib.sha256(snapshot['source'].encode()).hexdigest()==snapshot['sha256']
        assert hashlib.sha256((ROOT/snapshot['path']).read_bytes()).hexdigest()==snapshot['sha256']
    xx=s.symbols('h x y z w');records=[]
    for record in saved['cases'][1:]:
        started=time.monotonic();P=decode(record['point_brackets'])
        directions=[decode(data)for data in record['all6_actual_parameter_directions']]
        D,V,M,F=ambient_matrices(P,xx,directions)
        for name,A in [('full350_D',D),('vertical',V),('orbit_vertical_scale',M),('actual_full_source',F)]:
            assert hashlib.sha256(json.dumps(A.tolist(),default=str).encode()).hexdigest()==record['matrix_sha256'][name]
        reduced,pivots=dm(D).rref();reduced=reduced.to_Matrix();rank=len(pivots)
        kernel=[]
        for j in range(150):
            if j in pivots:continue
            col=s.zeros(150,1);col[j]=1
            for i,pivot in enumerate(pivots):col[pivot]=-reduced[i,j]
            kernel.append(col*s.ilcm(1,*[s.denom(v)for v in col]))
        K=s.Matrix.hstack(*kernel);assert D*K==s.zeros(350,150-rank)
        combined=list(dm(F.row_join(K)).rref()[1])
        fpivots=[j for j in combined if j<F.cols]
        indices=[j-F.cols for j in combined if j>=F.cols]
        count=len(indices);assert len(fpivots)==48 and count==102-rank
        normals=[]
        for index in indices:
            A=s.zeros(5)
            for (ij,e),value in zip(COLUMNS,K[:,index]):
                A[ij]+=value*s.prod(v**power for v,power in zip(xx,e))
            for ij in itertools.combinations(range(5),2):A[ij[::-1]]=-A[ij]
            normals.append(A)
        squares=[jacobi(sparse(A,xx))for A in normals]
        pairs=list(itertools.combinations_with_replacement(range(count),2))
        Q=s.Matrix.hstack(*[squares[i]if i==j else jacobi(sparse(normals[i]+normals[j],xx))-squares[i]-squares[j]for i,j in pairs])
        independent=D[:,list(pivots)]
        rows=list(dm(independent.T).rref()[1]);block=dm(D.extract(rows,list(pivots)))
        inverse=block.inv();correction=inverse.matmul(dm(Q.extract(rows,range(Q.cols)))).to_Matrix()
        obstruction=Q-independent*correction
        assert obstruction.extract(rows,range(Q.cols))==s.zeros(rank,Q.cols)
        qr,qpivots=dm(obstruction).rref();qr=qr.to_Matrix()
        n=s.symbols('n0:'+str(count));mons=[n[i]*n[j]for i,j in pairs]
        polynomials=[s.expand(sum(qr[i,j]*mons[j]for j in range(Q.cols)))for i in range(len(qpivots))]
        gb=s.groebner(polynomials,*n,order='grevlex')if polynomials else None
        isolated=bool(gb and gb.is_zero_dimensional)
        functional_rows=list(dm(obstruction.T).rref()[1]);functionals=[]
        for row in functional_rows:
            weights=dm(independent[row,:]).matmul(inverse).to_Matrix()
            functional=s.zeros(1,350);functional[row]=1
            for j,selected in enumerate(rows):functional[selected]-=weights[j]
            assert functional*D==s.zeros(1,150)
            coeff=list(functional*Q);assert coeff==list(obstruction[row,:])
            functionals.append(dict(row=row,nonzero_coefficients={str(j):str(v)for j,v in enumerate(functional)if v},quadratic_coefficients=list(map(str,coeff))))
        result=dict(case=record['case'],frozen_all_four_matrix_hashes_reproduced=True,
            full_D_rank=rank,actual_source_tangent_rank=48,full_kernel_dimension=150-rank,
            normal_count=count,normal_kernel_indices=indices,source_independent_column_indices=fpivots,
            normals=[{str(ij):str(A[ij])for ij in itertools.combinations(range(5),2)if A[ij]}for A in normals],
            normal_monomial_order=list(map(str,mons)),D_times_all_normals_zero=True,
            eliminated_D_rows=rows,eliminated_D_columns=list(pivots),invertible_D_minor_determinant=str(block.det()),
            quadratic_obstruction_rank=len(qpivots),quadratic_initial_generators=list(map(str,polynomials)),
            left_cokernel_functionals=functionals,grevlex_basis=[str(poly.as_expr())for poly in gb.polys]if gb else [],
            quadratic_support_zero_dimensional=isolated,
            local_projective_dimension_upper_if_zero_dimensional=23 if isolated else None,
            source_componenthood_if_formal_slice_argument_accepted=isolated,
            local_uniqueness_and_full_initial_ideal_not_claimed=True,
            obstruction_matrix_sha256=hashlib.sha256(json.dumps(obstruction.tolist(),default=str).encode()).hexdigest(),
            seconds=time.monotonic()-started)
        records.append(result)
        print(json.dumps({k:result[k]for k in ['case','normal_count','quadratic_initial_generators','quadratic_support_zero_dimensional','seconds']}),flush=True)
    inputs=[SOURCE,'code/classification_three_hyperplanes_quadric_independent_normals.py',
            'code/classification_three_hyperplanes_quadric_independent.py','code/classification_rank3_mixed_component_review.py']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,cases=records,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in inputs},
        producer_source=Path(__file__).read_text(),global_exhaustion_proved=False,seconds=time.monotonic()-begun)
    with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)

#!/usr/bin/env python3
"""One cubic obstruction after the unique one-normal quadratic correction."""
import argparse,hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,jacobi
from classification_three_hyperplanes_quadric_independent import ambient_matrices

ROOT=Path(__file__).resolve().parents[1]
SOURCE='reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_20261107.json'
NORMAL='reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_NORMALS_20261107.json'

def dm(A):return DomainMatrix.from_Matrix(A).convert_to(s.QQ)
def decode(data):
    A=s.zeros(5)
    for ij in itertools.combinations(range(5),2):
        A[ij]=s.sympify(data.get(str(ij),'0'));A[ij[::-1]]=-A[ij]
    return A
def sp(A,x):
    return {ij:{e:v for e,v in s.Poly(A[ij],*x).terms()if v}for ij in itertools.combinations(range(5),2)if A[ij]}
def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();x=s.symbols('h x y z w')
    saved=json.loads((ROOT/SOURCE).read_text())['cases'][1]
    normal=json.loads((ROOT/NORMAL).read_text())['cases'][0]
    assert saved['case']==normal['case']=='m1_t2'
    P=decode(saved['point_brackets']);R=decode(normal['normals'][0])
    D,V,M,F=ambient_matrices(P,x,[decode(a)for a in saved['all6_actual_parameter_directions']])
    assert hashlib.sha256(json.dumps(D.tolist(),default=str).encode()).hexdigest()==saved['matrix_sha256']['full350_D']
    cols=normal['eliminated_D_columns'];rows=normal['eliminated_D_rows'];independent=D[:,cols]
    inverse=dm(D.extract(rows,cols)).inv()
    q2=jacobi(sp(R,x));correction=-inverse.matmul(dm(q2.extract(rows,[0]))).to_Matrix()
    assert q2+independent*correction==s.zeros(350,1)
    S=s.zeros(5)
    for index,value in zip(cols,correction):
        ij,e=COLUMNS[index];S[ij]+=value*s.prod(v**power for v,power in zip(x,e))
    for ij in itertools.combinations(range(5),2):S[ij[::-1]]=-S[ij]
    q4=jacobi(sp(S,x));q3=jacobi(sp(R+S,x))-q2-q4
    cubic=q3-independent*inverse.matmul(dm(q3.extract(rows,[0]))).to_Matrix()
    candidates=[i for i,value in enumerate(cubic)if value]
    functionals=[]
    if candidates:
        # Prefer the sparsest single exact annihilator; no extra algebraic search.
        for row in candidates:
            weights=dm(independent[row,:]).matmul(inverse).to_Matrix()
            functional=s.zeros(1,350);functional[row]=1
            for j,selected in enumerate(rows):functional[selected]-=weights[j]
            data={str(j):str(value)for j,value in enumerate(functional)if value}
            functionals.append((len(data),row,functional,data))
        _,row,functional,data=min(functionals,key=lambda item:item[0])
        assert functional*D==s.zeros(1,150)and(functional*q3)[0]==cubic[row]
        witness=dict(row=row,nonzero_coefficients=data,cubic_pairing=str(cubic[row]),annihilates_full_D=True)
    else:witness=None
    inputs=[SOURCE,NORMAL,'code/classification_three_hyperplanes_quadric_m1_cubic.py',
            'code/classification_three_hyperplanes_quadric_independent.py','code/classification_rank3_mixed_component_review.py']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,case='m1_t2',
        all_source_and_normal_inputs_frozen=True,full_D_rank=101,actual_source_rank=48,normal_count=1,
        second_order_correction={str(ij):str(S[ij])for ij in itertools.combinations(range(5),2)if S[ij]},
        full_second_order_equations_zero=True,cubic_obstruction_nonzero=bool(candidates),
        cubic_left_cokernel_witness=witness,
        all_t_curve_if_order3_and_order4_zero=(q3==s.zeros(350,1)and q4==s.zeros(350,1)),
        raw_order3_nonzero_coefficients=sum(bool(v)for v in q3),raw_order4_nonzero_coefficients=sum(bool(v)for v in q4),
        source_componenthood_if_formal_slice_argument_accepted=bool(candidates),
        local_projective_dimension_upper_if_nonzero_cubic=23 if candidates else None,
        exact_nil_thickness_and_local_uniqueness_not_claimed=True,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in inputs},
        producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:result[k]for k in ['cubic_obstruction_nonzero','cubic_left_cokernel_witness','all_t_curve_if_order3_and_order4_zero','seconds']}))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)

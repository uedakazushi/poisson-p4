#!/usr/bin/env python3
"""Genuine cubic initials after cancelling all quadratic initial rows.

Each of the two remaining multi-normal sources receives one bounded third-order
test.  Cubic terms of equations with nonzero quadratic part are never themselves
treated as initial forms.
"""
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
def encode(A):return {str(ij):str(A[ij])for ij in itertools.combinations(range(5),2)if A[ij]}
def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(90,91));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();xx=s.symbols('h x y z w');records=[]
    saved=json.loads((ROOT/SOURCE).read_text())['cases']
    normal_saved=json.loads((ROOT/NORMAL).read_text())['cases']
    for record,normal in zip(saved[2:],normal_saved[1:]):
        started=time.monotonic();assert record['case']==normal['case'];count=normal['normal_count']
        P=decode(record['point_brackets']);R=[decode(a)for a in normal['normals']]
        D,V,M,F=ambient_matrices(P,xx,[decode(a)for a in record['all6_actual_parameter_directions']])
        assert hashlib.sha256(json.dumps(D.tolist(),default=str).encode()).hexdigest()==record['matrix_sha256']['full350_D']
        cols=normal['eliminated_D_columns'];rows=normal['eliminated_D_rows'];independent=D[:,cols]
        inverse=dm(D.extract(rows,cols)).inv()
        pairs=list(itertools.combinations_with_replacement(range(count),2))
        triples=list(itertools.combinations_with_replacement(range(count),3))
        squares=[jacobi(sp(A,xx))for A in R]
        Q2=s.Matrix.hstack(*[squares[i]if i==j else jacobi(sp(R[i]+R[j],xx))-squares[i]-squares[j]for i,j in pairs])
        correction=-inverse.matmul(dm(Q2.extract(rows,range(Q2.cols)))).to_Matrix()
        O2=Q2+independent*correction
        assert O2.extract(rows,range(Q2.cols))==s.zeros(len(rows),Q2.cols)
        S=[]
        for j in range(Q2.cols):
            A=s.zeros(5)
            for index,value in zip(cols,correction[:,j]):
                ij,e=COLUMNS[index];A[ij]+=value*s.prod(v**power for v,power in zip(xx,e))
            for ij in itertools.combinations(range(5),2):A[ij[::-1]]=-A[ij]
            S.append(A)
        ssquares=[jacobi(sp(A,xx))for A in S]
        Q3=s.zeros(350,len(triples))
        for i,A in enumerate(R):
            for j,B in enumerate(S):
                index=triples.index(tuple(sorted((i,*pairs[j]))))
                Q3[:,index]+=jacobi(sp(A+B,xx))-squares[i]-ssquares[j]
        O3=Q3-independent*inverse.matmul(dm(Q3.extract(rows,range(Q3.cols)))).to_Matrix()
        qrows=list(dm(O2.T).rref()[1]);assert len(qrows)==1
        qrow=qrows[0];qcolumn=next(i for i,v in enumerate(O2[qrow,:])if v)
        multipliers=O2[:,qcolumn]/O2[qrow,qcolumn]
        assert O2-multipliers*O2[qrow,:]==s.zeros(350,Q2.cols)
        G3=O3-multipliers*O3[qrow,:]
        reduced,pivots=dm(G3).rref();reduced=reduced.to_Matrix()
        n=s.symbols('n0:'+str(count));mons2=[n[i]*n[j]for i,j in pairs]
        mons3=[n[i]*n[j]*n[k]for i,j,k in triples]
        quadrics=[s.sympify(q)for q in normal['quadratic_initial_generators']]
        cubics=[s.expand(sum(reduced[i,j]*mons3[j]for j in range(len(mons3))))for i in range(len(pivots))]
        gb=s.groebner(quadrics+cubics,*n,order='grevlex')
        isolated=bool(gb.is_zero_dimensional)
        def functional(row):
            weights=dm(independent[row,:]).matmul(inverse).to_Matrix()
            lam=s.zeros(1,350);lam[row]=1
            for j,selected in enumerate(rows):lam[selected]-=weights[j]
            return lam
        qlam=functional(qrow);witnesses=[]
        for row in dm(G3.T).rref()[1]:
            lam=functional(row)-multipliers[row]*qlam
            assert lam*D==s.zeros(1,150)
            assert lam*Q2==s.zeros(1,Q2.cols)
            coeff=list(lam*Q3);assert coeff==list(G3[row,:])
            witnesses.append(dict(row=row,nonzero_coefficients={str(j):str(v)for j,v in enumerate(lam)if v},
                                  annihilates_full_D=True,annihilates_all_quadratic_normal_coefficients=True,
                                  cubic_coefficients=list(map(str,coeff))))
        result=dict(case=record['case'],normal_count=count,
            fixed_point_and_full_D_hash_reproduced=True,
            quadratic_normal_monomials=list(map(str,mons2)),cubic_normal_monomials=list(map(str,mons3)),
            second_order_correction_coefficients=[encode(A)for A in S],
            full_selected_linear_equations_solved_through_order2=True,
            full_quadratic_rowspace_rank=1,quadratic_pivot_row=qrow,
            all_cubic_forms_are_genuine_initials_after_constant_quadratic_cancellation=True,
            quadratic_initial_generators=list(map(str,quadrics)),cubic_initial_generators=list(map(str,cubics)),
            cubic_initial_rank=len(pivots),left_cokernel_functionals=witnesses,
            grevlex_basis=[str(poly.as_expr())for poly in gb.polys],
            quadratic_and_cubic_support_zero_dimensional=isolated,
            source_componenthood_if_formal_slice_argument_accepted=isolated,
            local_projective_dimension_upper_if_isolated=23 if isolated else None,
            longer_jets_not_authorized_by_this_result=True,seconds=time.monotonic()-started)
        records.append(result)
        print(json.dumps({k:result[k]for k in ['case','cubic_initial_rank','cubic_initial_generators','quadratic_and_cubic_support_zero_dimensional','seconds']}),flush=True)
    inputs=[SOURCE,NORMAL,'code/classification_three_hyperplanes_quadric_independent_cubics.py',
            'code/classification_three_hyperplanes_quadric_independent.py','code/classification_rank3_mixed_component_review.py']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,cases=records,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in inputs},
        producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)

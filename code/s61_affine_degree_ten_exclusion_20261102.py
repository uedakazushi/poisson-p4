#!/usr/bin/env python3
"""Exclude the five possible conjugate quadric cones by a universal pencil calculation."""
import argparse
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path
import time

import sympy as sp

from s61_model import S61, X, J

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_AFFINE_DEGREE_TEN_EXCLUSION_20261102'
ALLOWED=[0,9,11,16,25,26,27,28,30,32,33,35,37,42]


def recover_rows(matrix,targets):
    columns=list(matrix.to_DM().rref()[1])
    rows=list(matrix[:,columns].T.to_DM().rref()[1])
    coefficients=targets[:,columns]*matrix.extract(rows,columns).inv()
    assert coefficients*matrix[rows,:]==targets
    return rows,coefficients


def sparse(matrix):
    return [[[j,str(v)] for j,v in enumerate(row) if v] for row in matrix.tolist()]


def verify():
    start=time.monotonic()
    model=S61()
    x,_,z,_=X[1:]
    weights=model.reduced['parameter_weights']
    one=[i for i,w in enumerate(weights) if w==1]
    two=[i for i,w in enumerate(weights) if w==2]
    def quadratic(polynomial):
        representative=X[0]**2*J
        representative[2,4]-=polynomial
        representative[4,2]+=polynomial
        return sp.Matrix([model.normalize(representative)['parameters'][i] for i in two])
    A,C=quadratic(x*x),quadratic(z*z)
    assert A==-sp.eye(45)[:,25] and C==-sp.eye(45)[:,32]
    beta={i:sp.Symbol('u'+str(i+1))for i in one}
    beta[3]=beta[4]=0
    beta[13]=sp.Rational(4,5)*beta[48]
    h=beta[50]-5*beta[25]
    variables=[beta[i] for i in one if i not in (3,4,13)]
    assert len(variables)==13
    H=[];linear=[]
    for row,equation in enumerate(model.reduced['equations']):
        hh=sp.zeros(45);linear_row=sp.zeros(45,1)
        for mon,num,den in equation:
            assert sum(weights[i] for i in mon)==(4 if row<70 else 5)
            coordinates=[two.index(i) for i in mon if i in two]
            coeff=sp.Rational(num,den)*sp.prod(beta[i] for i in mon if i in beta)
            if len(coordinates)==2:
                i,j=coordinates
                hh[i,j]+=coeff/2
                hh[j,i]+=coeff/2
            if len(coordinates)==1:
                linear_row[coordinates[0]]+=coeff
        H.append(hh);linear.append(linear_row)
    jac=lambda q:sp.Matrix.vstack(*[2*q.T*h0 for h0 in H])
    da,dc=jac(A),jac(C)
    first=da[:70,:].col_join(dc[:70,:])
    kernel=sp.Matrix.hstack(*first.nullspace())
    assert first.rank()==17 and kernel.shape==(45,28)
    remainder=da[70:,:].col_join(dc[70:,:])*kernel
    # For each original row, concatenate all thirteen coefficient rows in the same order.
    stack=sp.Matrix([[expr.diff(v) for v in variables for expr in row]
                     for row in remainder.tolist()])
    assert stack.shape==(128,364)
    excluded=[i for i in range(45) if i not in ALLOWED and any(kernel.row(i))]
    targets=sp.Matrix([[h.diff(v)*entry for v in variables for entry in kernel.row(i)]
                       for i in excluded])
    selected,multipliers=recover_rows(stack,targets)
    assert len(selected)==30 and len(excluded)==25
    # Verify again over QQ[u], without the coefficient-stack encoding.
    assert (multipliers*remainder[selected,:]-h*kernel[excluded,:]).applyfunc(sp.expand)==sp.zeros(25,28)
    for i in range(45):
        if i not in ALLOWED and i not in excluded:
            assert not any(kernel.row(i))
    # For these rows the B^2 term vanishes identically on the certified coordinate subspace.
    keep=[i for i in range(70) if not any(H[i].extract(ALLOWED,ALLOWED))]
    d=da[:70,:].col_join(dc[keep,:])
    shift=sp.Matrix([(A.T*linear[i])[0] for i in range(70)]+[(C.T*linear[i])[0]for i in keep])
    pc=list(d.to_DM().rref()[1])
    pr=list(d[:,pc].T.to_DM().rref()[1])
    lift=d[:,pc]*d.extract(pr,pc).inv()
    assert lift*d[pr,:]==d
    residual=(shift-lift*shift[pr,:]).applyfunc(sp.expand)
    polynomials=[sp.Poly(f,*variables)for f in residual]
    target=sp.Poly(h*h,*variables)
    monomials=sorted(set().union(*(set(f.monoms())for f in polynomials+[target])))
    matrix=sp.Matrix([[f.coeff_monomial(mon)for mon in monomials]for f in polynomials])
    indices,square=recover_rows(matrix,sp.Matrix([[target.coeff_monomial(mon)for mon in monomials]]))
    assert sp.expand((square*residual[indices,:])[0]-h*h)==0
    # Express h^2 directly by the original selected affine-linear coefficient equations.
    row=sp.zeros(1,len(shift))
    for j,k in enumerate(indices):
        row[k]+=square[j]
        for m,n in enumerate(pr):
            row[n]-=square[j]*lift[k,m]
    assert row*d==sp.zeros(1,45) and sp.expand((row*shift)[0]-h*h)==0
    deps=['S61_AFFINE_SURFACE_DEGREE_OBSTRUCTION_20261102','S61_GENERIC_BOUNDARY_INDEX_20261102',
          'S61_AFFINE_SURFACE_BOUND_20261102','S61_BINARY_H11_RANK_DROP_RESIDUAL_FINITE_20261102']
    files=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),ROOT/'code/s61_model.py']
    for name in deps:
        f=ROOT/'reports'/(name+'.json')
        assert json.loads(f.read_text())['all_passed']
        files.extend([f,ROOT/'research/notes'/(name+'.md')])
    files.extend(ROOT/'legacy/v02/data'/n for n in ['symplectic_slice_certificate.json',
                'symplectic_slice_equations.json','symplectic_slice_reduced_equations.json'])
    return dict(id=STEM,all_passed=True,field='QQ identities, generic base K=C(W16)',
                quadratic_cone_normal_form='q=E+s*A+t*B+t^2*C',
                directions_A_C=[sparse(A.T)[0],sparse(C.T)[0]],
                normalized_source=['u4=u5=0','5*u14=4*u49'],remaining_variables=list(map(str,variables)),
                invertible_marked_factor=str(h),first_kernel_dimension=28,first_kernel=sparse(kernel),
                lower_stack_shape=[128,364],constant_row_span_rank=len(selected),
                selected_lower_rows=selected,excluded_actual_B_coordinates=excluded,
                horizontal_multipliers=sparse(multipliers),remaining_B_coordinates=ALLOWED,
                zero_B2_rows_at_C=keep,stacked_constant_rows_order='s coefficient rows 0..69; t^2 coefficient rows in zero_B2_rows_at_C',
                direct_square_witness=sparse(row)[0],square_identity=str(h*h),
                actual_degree_ten_surface_over_K_excluded=True,
                nonzero_surface_degree_over_K_is_multiple_of_five_and_at_least=15,
                two_variable_projection_rank_if_some_surface_at_least=15,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files},
                seconds=time.monotonic()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=verify()
    with args.output.open('x') as out:
        json.dump(report,out,indent=2)
        out.write('\n')
    print('PASS',len(report['direct_square_witness']),round(report['seconds'],3))

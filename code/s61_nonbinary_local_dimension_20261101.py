#!/usr/bin/env python3
"""Local dimension one at a nonbinary constant-factor point at infinity.

Use the original 134 quadratic homogenizations as primary input; the 135th
consequence supplies an independent smaller tangent-space cross-check.
The geometric lower bound is supplied by the smooth source/submersion and
unique constant-factor line in S61_R0_INFINITY_20261101.
"""
from __future__ import annotations
import argparse,hashlib,json
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize,EXTRA
ROOT=Path(__file__).resolve().parents[1]


def calculate(model,extra):
    point=s.zeros(61,1)
    for index,value in {14:1,49:3,26:2,35:3,50:4,51:5,28:s.Rational(-5,7),38:s.Rational(19,2),41:-1}.items():
        point[index-1]=value
    beta={i:point[i] for i,w in enumerate(model.reduced['parameter_weights']) if w==1}
    two,polynomials=specialize(model,beta,extra)
    q0=s.Matrix([point[i] for i in two]);fixed=two.index(27)
    columns=[i for i in range(45) if i!=fixed];count=len(polynomials)
    derivative=s.zeros(count,45);inhomogeneous=s.zeros(count,1)
    for row,f in enumerate(polynomials):
        value=0
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon;value+=c*q0[i]*q0[j]
                derivative[row,i]+=c*q0[j];derivative[row,j]+=c*q0[i]
            elif len(mon)==1:inhomogeneous[row,0]+=c*q0[mon[0]]
        assert value==0
    jacobian=derivative[:,columns].row_join(inhomogeneous)
    tangent=s.Matrix.hstack(*jacobian.nullspace())
    cokernel=s.Matrix.hstack(*jacobian.T.nullspace()).T
    n=tangent.cols
    assert tangent[-1,:]==s.zeros(1,n)
    linear_forms={i:list(tangent[columns.index(i),:]) if i!=fixed else [0]*n for i in range(45)}
    t=list(tangent[-1,:]);monomials=list(combinations_with_replacement(range(n),2));coefficients=[]
    for polynomial in polynomials:
        out={}
        for mon,c in polynomial.items():
            forms=[linear_forms[i] for i in mon]+[t]*(2-len(mon))
            for i in range(n):
                for j in range(n):
                    key=tuple(sorted((i,j)));out[key]=out.get(key,0)+c*forms[0][i]*forms[1][j]
        coefficients.append([out.get(mon,0) for mon in monomials])
    projected=(cokernel.to_DM()*s.Matrix(coefficients).to_DM()).to_Matrix()
    _,pivots=projected.T.to_DM().rref();selected=list(pivots);basis=projected[selected,:]
    hessians=[]
    for row in basis.tolist():
        f=dict(zip(monomials,row))
        hessians.extend([[f[tuple(sorted((i,j)))]*(2 if i==j else 1) for j in range(n)] for i in range(n)])
    hessian=s.Matrix(hessians);radical=s.Matrix.hstack(*hessian.nullspace())
    expected=(42,3,3,2) if not extra else (43,2,1,1)
    assert (jacobian.rank(),n,basis.rows,hessian.rank())==expected
    assert radical.cols==1 and basis.rows==(n-1)*n//2
    assert jacobian*tangent==s.zeros(count,n)
    assert cokernel*jacobian==s.zeros(cokernel.rows,45)
    def rows(matrix):return [[str(v) for v in row] for row in matrix.tolist()]
    return dict(equation_count=count,point61=list(map(str,point)),q_to_u_indices_zero_based=two,
        chart_fixed_u=28,chart_fixed_value='-5/7',chart_columns=columns,
        leading_jacobian_rank=derivative.rank(),homogenized_jacobian_rank=jacobian.rank(),
        jacobian=rows(jacobian),tangent_basis=rows(tangent),homogenizer_zero_on_tangent=True,
        cokernel_rows=rows(cokernel[selected,:]),quadratic_monomials=monomials,
        quadratic_obstruction_basis=rows(basis),quadratic_common_radical=rows(radical),
        quadratic_span_rank=basis.rows,common_hessian_rank=hessian.rank(),local_dimension_upper_bound=1)


def verify():
    model=S61();primary=calculate(model,False);augmented=calculate(model,True)
    source=ROOT/'reports/S61_R0_INFINITY_20261101.json'
    lower=json.loads(source.read_text())
    incidence=lower['nonbinary']['beta_line_incidence_minor']
    assert incidence['rank']==17
    beta=lower['nonbinary']['dominance_minor']
    assert beta['rank']==16
    inputs=[Path(__file__),ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
        ROOT/'code/s61_model.py',EXTRA,source,ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json']
    return dict(id='S61_NONBINARY_LOCAL_DIMENSION_20261101',all_passed=True,field='QQ',
        primary134=primary,augmented135=augmented,local_dimension=1,
        lower_bound='The smooth source of the R0 family maps submersively to W at the displayed point (rank16), and jointly to (beta,[v]) with rank17. Its fixed-beta fiber therefore maps nontrivially through the point. The constant factor line is uniquely and regularly recovered from Q on this open, so this gives a curve in the q fiber.',
        upper_bound='Quadratic cokernel obstructions span all quadratic products in the normal space to one tangent line; the tangent cone support has dimension at most1.',
        limitations=['Local dimension of the unsaturated homogenization only.','The reduced local ring is not determined; additional tangent branches have not been excluded.','No generic affine finiteness, degree, or global exhaustion conclusion.'],
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    report=verify();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,local_dimension=1,primary_tangent_dimension=3,augmented_tangent_dimension=2)))
if __name__=='__main__':main()

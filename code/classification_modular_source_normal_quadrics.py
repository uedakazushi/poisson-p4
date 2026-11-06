#!/usr/bin/env python3
"""True normal quadratic support for one previously verified complete source."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as S
from sympy.polys.matrices import DomainMatrix

import classification_rank3_mixed_component_review as arithmetic

ROOT=Path(__file__).resolve().parents[1]


def sparse(p,coordinates):
    return {ij:dict(S.Poly(p[ij],*coordinates).terms()) for ij in arithmetic.PAIRS if p[ij]!=0}


def dm(mat):
    return DomainMatrix.from_Matrix(mat).convert_to(S.QQ)


def main(source_path,branch_name,output):
    began=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    signal.alarm(90)
    report=json.loads(source_path.read_text())
    branch=next(b for b in report['branches'] if b['name']==branch_name)
    coords=S.symbols('h x1 x2 x3 x4')
    spec=branch['full_source_specification']
    parameters=tuple(S.Symbol(v) for v in spec['free_parameters'])
    env={str(v):v for v in (*coords,*parameters)}
    family=S.zeros(5)
    for ij in arithmetic.PAIRS:
        family[ij]=S.sympify(branch['universal_brackets'].get(f'{ij[0]},{ij[1]}','0'),locals=env)
        family[ij[::-1]]=-family[ij]
    point={env[k]:S.Rational(v) for k,v in spec['witness'].items()}
    P=family.subs(point).applyfunc(S.expand)
    ps=sparse(P,coords)
    assert arithmetic.jacobi(ps)==S.zeros(350,1)
    cols=[]
    for ij,e in arithmetic.COLUMNS:
        perturbed={pair:dict(coeffs) for pair,coeffs in ps.items()}
        perturbed.setdefault(ij,{})[e]=perturbed.get(ij,{}).get(e,0)+1
        cols.append(arithmetic.jacobi(perturbed))
    D=S.Matrix.hstack(*cols)
    X=S.Matrix(coords)
    actual=[]
    for a,b in itertools.product(range(5),repeat=2):
        linear=S.zeros(5)
        linear[a,b]=1
        field=linear*X
        actual.append(arithmetic.vector(sparse(X*field.T-field*X.T,coords)))
    for a,b in itertools.product(range(5),repeat=2):
        linear=S.zeros(5)
        linear[a,b]=1
        action=linear*P+P*linear.T-coords[b]*P.diff(coords[a])
        actual.append(arithmetic.vector(sparse(action,coords)))
    actual.append(arithmetic.vector(ps))
    actual.extend(arithmetic.vector(sparse(family.diff(v).subs(point),coords)) for v in parameters)
    F=S.Matrix.hstack(*actual)
    assert D*F==S.zeros(350,F.cols)
    rr,pivots=dm(D).rref()
    rr=rr.to_Matrix()
    rankd=len(pivots)
    free=[j for j in range(150) if j not in pivots]
    kernel=[]
    for j in free:
        col=S.zeros(150,1)
        col[j]=1
        for row,pivot in enumerate(pivots):
            col[pivot]=-rr[row,j]
        denominator=S.ilcm(1,*[S.denom(v) for v in col])
        kernel.append(col*denominator)
    K=S.Matrix.hstack(*kernel)
    assert D*K==S.zeros(350,K.cols)
    full_pivots=list(dm(F.row_join(K)).rref()[1])
    source_pivots=[j for j in full_pivots if j<F.cols]
    normal_indices=[j-F.cols for j in full_pivots if j>=F.cols]
    normal_count=len(normal_indices)
    assert 0<normal_count<=5
    normals=[]
    normal_sparse=[]
    for j in normal_indices:
        data={}
        for (ij,e),v in zip(arithmetic.COLUMNS,K[:,j]):
            if v:
                data.setdefault(ij,{})[e]=v
        normal_sparse.append(data)
        R=S.zeros(5)
        for ij,poly in data.items():
            R[ij]=sum(c*S.prod(x**a for x,a in zip(coords,e)) for e,c in poly.items())
            R[ij[::-1]]=-R[ij]
        normals.append(R)
    pairs=list(itertools.combinations_with_replacement(range(normal_count),2))
    squares=[arithmetic.jacobi(r) for r in normal_sparse]
    qcols=[squares[i] if i==j else arithmetic.jacobi(sparse(normals[i]+normals[j],coords))-squares[i]-squares[j]
           for i,j in pairs]
    Q=S.Matrix.hstack(*qcols)
    independent_D=D[:,list(pivots)]
    rows=list(dm(independent_D.T).rref()[1])
    block=dm(D.extract(rows,list(pivots)))
    correction=block.inv().matmul(dm(Q.extract(rows,range(len(pairs))))).to_Matrix()
    obstruction=Q-independent_D*correction
    assert obstruction.extract(rows,range(len(pairs)))==S.zeros(rankd,len(pairs))
    reduced,q_pivots=dm(obstruction).rref()
    reduced=reduced.to_Matrix()
    qrank=len(q_pivots)
    variables=S.symbols(f't0:{normal_count}')
    monomials=[variables[i]*variables[j] for i,j in pairs]
    quadrics=[S.expand(sum(reduced[i,j]*monomials[j] for j in range(len(pairs)))) for i in range(qrank)]
    gb=S.groebner(quadrics,*variables,order='grevlex') if quadrics else None
    zero_dimensional=bool(gb and gb.is_zero_dimensional)
    files=[Path(__file__).resolve(),Path(arithmetic.__file__).resolve(),source_path.resolve()]
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,branch=branch_name,
        exact_D_rank=rankd,actual_source_tangent_rank=len(source_pivots),
        full_kernel_dimension=K.cols,normal_count=normal_count,
        normal_indices=normal_indices,source_independent_column_indices=source_pivots,
        all_actual_source_and_normal_directions_killed=True,
        normal_monomial_order=list(map(str,monomials)),
        normal_bivectors=[{f'{i},{j}':str(R[i,j]) for i,j in arithmetic.PAIRS if R[i,j]!=0} for R in normals],
        genuine_quadratic_initial_generators=list(map(str,quadrics)),
        quadratic_obstruction_rank=qrank,
        groebner_basis=[str(p.as_expr()) for p in gb.polys] if gb else [],
        quadratic_support_zero_dimensional=zero_dimensional,
        projective_component_dimension_if_zero_dimensional=len(source_pivots)-25 if zero_dimensional else None,
        all_positive_dimensional_branches_retained=True,
        higher_initial_ideal_or_exact_nilpotence_claimed=False,
        eliminated_D_rows=rows,eliminated_D_columns=list(pivots),
        matrix_sha256={label:hashlib.sha256(json.dumps(mat.tolist(),default=str).encode()).hexdigest()
                       for label,mat in [('D',D),('source',F),('obstruction',obstruction)]},
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        seconds=time.monotonic()-began)
    with output.open('x') as handle:
        json.dump(result,handle,indent=2)
        handle.write('\n')
    signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['branch','exact_D_rank','actual_source_tangent_rank','normal_count',
        'genuine_quadratic_initial_generators','quadratic_support_zero_dimensional','seconds']}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_report',type=Path)
    parser.add_argument('branch')
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    main(args.source_report,args.branch,args.output)

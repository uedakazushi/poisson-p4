#!/usr/bin/env python3
"""Exact finite equivariant implicit slice, not a precision-dependent jet.

The stable linear slice is constructed in all150 ambient coefficients. A
negative-weight one-dimensional normal forces all eliminated coordinates to
be finite monomials; the whole formal support therefore has a finite bound.
"""
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
import classification_rank3_mixed_component_review as engine


def dm(M):
    return DomainMatrix.from_Matrix(M).convert_to(s.QQ)


def sparse(P,xx):
    return {ij:dict(s.Poly(s.expand(P[ij]),*xx).terms())
            for ij in engine.PAIRS if P[ij]!=0}


def tensor(v,xx):
    M=s.zeros(5)
    for ((i,j),exponent),coefficient in zip(engine.COLUMNS,v):
        M[i,j]+=coefficient*s.prod(x**e for x,e in zip(xx,exponent))
    return M-M.T


def brackets(P):
    return {f'{i},{j}':str(s.factor(P[i,j])) for i,j in engine.PAIRS if P[i,j]!=0}


def coefficient(J,tau,k):
    return s.Matrix([s.Poly(v,tau).coeff_monomial(tau**k) for v in J])


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(90,91))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    started=time.monotonic()
    path=Path('reports/CLASSIFICATION_NIL22_EQUIVARIANT_20261107.json')
    source=json.loads(path.read_text())
    xx=s.symbols('h x1 x2 x3 x4');tau=s.Symbol('tau')
    env={str(v):v for v in xx};env['n1']=s.Integer(1)
    P=s.zeros(5)
    for key,value in source['canonical_brackets'].items():
        i,j=map(int,key.split(','));P[i,j]=s.sympify(value,locals=env)
    P-=P.T
    ps=sparse(P,xx);assert engine.jacobi(ps)==s.zeros(350,1)
    weights=source['weight'];pweight=source['P_weight']
    absolute=[weights[i]+weights[j]-sum(a*b for a,b in zip(weights,e)) for (i,j),e in engine.COLUMNS]
    relative=[v-pweight for v in absolute]
    assert min(relative)==-8 and max(relative)==6
    assert all(v==0 or a==pweight for v,a in zip(engine.vector(ps),absolute))
    Dcols=[]
    for ij,mon in engine.COLUMNS:
        v={key:dict(value) for key,value in ps.items()}
        v.setdefault(ij,{})[mon]=v.get(ij,{}).get(mon,0)+1
        Dcols.append(engine.jacobi(v))
    D=s.Matrix.hstack(*Dcols)
    V=[];O=[];E=s.Matrix(xx)
    for a,b in itertools.product(range(5),repeat=2):
        A=s.zeros(5);A[a,b]=1;field=A*E
        V.append(engine.vector(sparse(E*field.T-field*E.T,xx)))
        O.append(engine.vector(sparse(A*P+P*A.T-xx[b]*P.diff(xx[a]),xx)))
    G=s.Matrix.hstack(*V,*O,engine.vector(ps))
    assert D*G==s.zeros(350,G.cols)
    # Each group direction is homogeneous, so the actual tangent is stable.
    for col in range(G.cols):
        assert len({relative[i] for i in range(150) if G[i,col]})<=1
    reduced,pivot_rows=dm(G.T).rref();reduced=reduced.to_Matrix()
    assert len(pivot_rows)==47
    for row in range(47):
        assert len({relative[i] for i in range(150) if reduced[row,i]})==1
    slice_indices=[i for i in range(150) if i not in pivot_rows]
    assert len(slice_indices)==103
    DS=D[:,slice_indices]
    implicit_pivots=list(dm(DS).rref()[1]);assert len(implicit_pivots)==102
    implicit_indices=[slice_indices[i] for i in implicit_pivots]
    missing=[i for i in range(103) if i not in implicit_pivots];assert len(missing)==1
    K=dm(DS).nullspace().to_Matrix().T;assert K.shape==(103,1)
    K/=K[missing[0]]
    rv=s.zeros(150,1)
    for i,index in enumerate(slice_indices):rv[index]=K[i]
    rweights={relative[i] for i in range(150) if rv[i]};assert rweights=={-1}
    assert D*rv==s.zeros(350,1)
    R=tensor(rv,xx)
    J=D[:,implicit_indices]
    rows=list(dm(J.T).rref()[1]);assert len(rows)==102
    inverse=dm(J[rows,:]).inv().to_Matrix()
    determinant=dm(J[rows,:]).det();assert determinant
    curve=P+tau*R
    corrections=[]
    obstruction=None
    bound=-min(relative)

    def certify_obstruction(F,k):
        nz=next(i for i,v in enumerate(F) if v)
        functional=s.zeros(1,350);functional[0,nz]=1
        selected=J[nz,:]*inverse
        for i,row in enumerate(rows):functional[0,row]-=selected[0,i]
        assert functional*D==s.zeros(1,150)
        value=(functional*F)[0];assert value
        return dict(first_nonzero_degree=k,
                    left_full_D_zero=True,
                    left_functional={str(i):str(functional[0,i]) for i in range(350) if functional[0,i]},
                    nonzero_pairing=str(value),
                    full_residual_coefficient={str(i):str(F[i]) for i in range(350) if F[i]})

    # These are all possible monomials of the UNIQUE equivariant IFT map.
    # The bound comes from the full ambient representation, not a jet budget.
    for k in range(2,bound+1):
        actual=engine.jacobi(sparse(curve,xx)).applyfunc(s.expand)
        F=coefficient(actual,tau,k)
        values=-inverse*F[rows,:]
        qv=s.zeros(150,1)
        for i,index in enumerate(implicit_indices):
            assert values[i]==0 or relative[index]==-k
            qv[index]=values[i]
        residual=F+D*qv
        Q=tensor(qv,xx)
        curve+=tau**k*Q
        if Q!=s.zeros(5):corrections.append(dict(degree=k,brackets=brackets(Q)))
        assert residual[rows,:]==s.zeros(102,1)
        if residual!=s.zeros(350,1):
            obstruction=certify_obstruction(residual,k)
            break
    actual=engine.jacobi(sparse(curve,xx)).applyfunc(s.expand)
    if obstruction is None:
        assert actual[rows,:]==s.zeros(102,1)
        nonzero_degrees=[mon[0] for v in actual if v for mon,c in s.Poly(v,tau).terms() if c]
        if nonzero_degrees:
            k=min(nonzero_degrees);assert k>bound
            obstruction=certify_obstruction(coefficient(actual,tau,k),k)
    if obstruction:
        k=obstruction['first_nonzero_degree']
        assert all(coefficient(actual,tau,j)==s.zeros(350,1) for j in range(k))
        status='FULL_REDUCED_SLICE_IS_POINT'
        component=22
    else:
        assert actual==s.zeros(350,1)
        status='EXACT_TRANSVERSE_POLYNOMIAL_CURVE'
        component=23
    result=dict(task_id='E_multiple_hyperplane',status=status,
                source_identity_report=str(path),source_identity_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                canonical_brackets=brackets(P),coordinate_weights=weights,
                relative_action='lambda^-1 times diagonal pushforward',
                relative_ambient_weight_range=[min(relative),max(relative)],normal_relative_weight=-1,
                all_implicit_functions_polynomial_degree_bound=bound,
                full_D_rank=102,actual_orbit_vertical_scale_rank=47,
                stable_slice_dimension=103,stable_slice_coordinate_indices=slice_indices,
                actual_group_tangent_pivot_rows=list(pivot_rows),
                implicit_coordinate_indices=implicit_indices,selected_equation_rows=rows,
                implicit_minor_determinant=str(determinant),
                full_kernel_dimension=48,normal_brackets=brackets(R),
                implicit_corrections=corrections,full_curve_brackets=brackets(curve),
                obstruction=obstruction,projective_full_component_dimension=component,
                exact_nilpotence_thickness_not_claimed=True,
                producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                producer_bytes=Path(__file__).read_text(),
                engine_sha256=hashlib.sha256(Path(engine.__file__).read_bytes()).hexdigest(),
                seconds=round(time.monotonic()-started,4))
    Path('reports/CLASSIFICATION_MODULAR_WEIGHTED_SLICE_20261107.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=status,component_dimension=component,obstruction=obstruction,seconds=result['seconds'])))


if __name__=='__main__':main()

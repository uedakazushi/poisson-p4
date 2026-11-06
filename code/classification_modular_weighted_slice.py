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
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args()
    if args.output.exists():raise FileExistsError('Preserve existing certificates: '+str(args.output))
    resource.setrlimit(resource.RLIMIT_CPU,(85,86))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    started=time.monotonic()
    source=json.loads(args.input.read_text())
    xx=s.symbols('h x1 x2 x3 x4')
    env={str(v):v for v in xx};env['n1']=s.Integer(1)
    env.update({k:s.Rational(v) for k,v in source.get('canonical_parameter_values',{}).items()})
    P=s.zeros(5)
    for key,value in source['canonical_brackets'].items():
        i,j=map(int,key.split(','));P[i,j]=s.sympify(value,locals=env)
    P-=P.T;ps=sparse(P,xx)
    assert engine.jacobi(ps)==s.zeros(350,1)
    weights=source['weight'];pweight=source['P_weight']
    absolute=[weights[i]+weights[j]-sum(a*b for a,b in zip(weights,e)) for (i,j),e in engine.COLUMNS]
    relative=[v-pweight for v in absolute]
    assert all(not v or a==pweight for v,a in zip(engine.vector(ps),absolute))
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
    for col in range(G.cols):assert len({relative[i] for i in range(150) if G[i,col]})<=1
    reduced,pivot_rows=dm(G.T).rref();reduced=reduced.to_Matrix();gdim=len(pivot_rows)
    for row in range(gdim):assert len({relative[i] for i in range(150) if reduced[row,i]})==1
    slice_indices=[i for i in range(150) if i not in pivot_rows]
    DS=D[:,slice_indices]
    implicit_pivots=list(dm(DS).rref()[1]);drank=len(implicit_pivots)
    implicit_indices=[slice_indices[i] for i in implicit_pivots]
    missing=[i for i in range(len(slice_indices)) if i not in implicit_pivots]
    K=dm(DS).nullspace().to_Matrix().T
    K=K*K[missing,:].inv()
    normals=[];normal_weights=[]
    for j in range(K.cols):
        rv=s.zeros(150,1)
        for i,index in enumerate(slice_indices):rv[index]=K[i,j]
        rw={relative[i] for i in range(150) if rv[i]};assert len(rw)==1
        assert D*rv==s.zeros(350,1)
        normals.append(tensor(rv,xx));normal_weights.append(next(iter(rw)))
    if 'expected_full_D_rank' in source:assert drank==source['expected_full_D_rank']
    if 'expected_actual_group_rank' in source:assert gdim==source['expected_actual_group_rank']
    base=dict(task_id='E_multiple_hyperplane',branch=source.get('branch','unspecified'),
              source_identity_report=str(args.input),source_identity_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
              canonical_brackets=brackets(P),coordinate_weights=weights,P_weight=pweight,
              relative_action=f'lambda^(-{pweight}) times diagonal pushforward',
              relative_ambient_weight_range=[min(relative),max(relative)],normal_relative_weights=normal_weights,
              full_D_rank=drank,actual_orbit_vertical_scale_rank=gdim,
              stable_slice_dimension=len(slice_indices),stable_slice_coordinate_indices=slice_indices,
              actual_group_tangent_pivot_rows=list(pivot_rows),implicit_coordinate_indices=implicit_indices,
              full_kernel_dimension=150-drank,normal_brackets=[brackets(R) for R in normals],
              producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),producer_bytes=Path(__file__).read_text(),
              engine_sha256=hashlib.sha256(Path(engine.__file__).read_bytes()).hexdigest(),
              exact_nilpotence_thickness_not_claimed=True)
    def save(extra):
        result={**base,**extra,'seconds':round(time.monotonic()-started,4)}
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({k:result[k] for k in ['status','normal_relative_weights','seconds']}))
    if not normal_weights:
        save(dict(status='SMOOTH_ORBIT_COMPONENT',projective_local_dimension=gdim-25));return
    if not (all(w<0 for w in normal_weights) or all(w>0 for w in normal_weights)):
        save(dict(status='NORMAL_WEIGHTS_NOT_STRICTLY_ONE_SIGN',full_slice_support_unresolved=True));return
    sign=-1 if normal_weights[0]<0 else 1
    costs=[sign*w for w in normal_weights]
    bound=max(sign*w for w in relative)
    tau=s.symbols('tau0:'+str(len(normals)))
    J=D[:,implicit_indices]
    rows=list(dm(J.T).rref()[1]);assert len(rows)==drank
    inverse=dm(J[rows,:]).inv().to_Matrix()
    determinant=dm(J[rows,:]).det();assert determinant
    base.update(selected_equation_rows=rows,implicit_minor_determinant=str(determinant),
                normal_variables=list(map(str,tau)),positive_normal_costs=costs,
                all_implicit_functions_weighted_degree_bound=bound)
    curve=P+sum((t*R for t,R in zip(tau,normals)),s.zeros(5))
    corrections=[]
    candidates=list(itertools.product(*[range(bound//c+1) for c in costs]))
    candidates=[e for e in candidates if sum(e)>=2 and sum(a*b for a,b in zip(e,costs))<=bound]
    candidates.sort(key=lambda e:(sum(a*b for a,b in zip(e,costs)),e))
    for degree in sorted({sum(a*b for a,b in zip(e,costs)) for e in candidates}):
        actual=engine.jacobi(sparse(curve,xx)).applyfunc(s.expand)
        polynomials=[s.Poly(f,*tau) for f in actual]
        for exponent in [e for e in candidates if sum(a*b for a,b in zip(e,costs))==degree]:
            F=s.Matrix([poly.coeff_monomial(exponent) for poly in polynomials])
            values=-inverse*F[rows,:]
            qv=s.zeros(150,1)
            for i,index in enumerate(implicit_indices):
                assert values[i]==0 or relative[index]==sign*degree
                qv[index]=values[i]
            assert (F+D*qv)[rows,:]==s.zeros(drank,1)
            Q=tensor(qv,xx)
            if Q!=s.zeros(5):
                curve+=s.prod(t**e for t,e in zip(tau,exponent))*Q
                corrections.append(dict(normal_exponent=exponent,weighted_degree=degree,brackets=brackets(Q)))
    actual=engine.jacobi(sparse(curve,xx)).applyfunc(s.expand)
    assert actual[rows,:]==s.zeros(drank,1)
    nonzero=[]
    for value in actual:
        if value:
            poly=s.Poly(value,*tau,domain=s.QQ).monic().as_expr()
            if poly not in nonzero:nonzero.append(poly)
    basis=s.groebner(nonzero,*tau,order='grevlex',domain=s.QQ) if nonzero else None
    if basis is None:dim=len(normals)
    elif basis.is_zero_dimensional:dim=0
    elif len(normals)<=2:dim=len(normals)-1
    else:dim=None
    assert basis is None or not any(p==1 for p in basis.polys)
    support=dict(status='FINITE_FULL_TRANSVERSE_IDEAL',implicit_corrections=corrections,
                 full_curve_brackets=brackets(curve),all_selected_IFT_equations_identically_zero=True,
                 full350_residual_polynomials={str(i):str(s.factor(v)) for i,v in enumerate(actual) if v},
                 residual_ideal_generators=list(map(str,nonzero)),
                 residual_groebner_basis=[] if basis is None else list(map(str,basis.exprs)),
                 groebner_order='grevlex',reduced_slice_support_dimension=dim,
                 reduced_slice_supported_only_at_origin=dim==0,
                 projective_local_dimension=None if dim is None else gdim+dim-25,
                 full_finite_IFT_not_a_truncated_jet=True)
    if dim==0:support['status']='FULL_REDUCED_SLICE_IS_POINT';support['projective_orbit_closure_component_dimension']=gdim-25
    save(support)


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Canonicalize one complete modular source and construct equivariant curves.

This transports the existing certified point and normal; no new point is probed.
The maximum correction degree is two, exactly as in the earlier finite test.
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


def sparse(P, xx):
    return {ij: dict(s.Poly(s.expand(P[ij]), *xx).terms())
            for ij in engine.PAIRS if P[ij] != 0}


def brackets(P):
    return {f'{i},{j}': str(s.factor(P[i, j])) for i, j in engine.PAIRS if P[i,j] != 0}


def tensor(v, xx):
    M = s.zeros(5)
    for ((i,j), exponent), coefficient in zip(engine.COLUMNS, v):
        M[i,j] += coefficient*s.prod(x**e for x,e in zip(xx,exponent))
    return M-M.T


def primitive(v):
    den = s.ilcm(*[a.q for a in v])
    v = v*den
    return v/s.igcd(*[int(a) for a in v])


def main():
    resource.setrlimit(resource.RLIMIT_CPU, (90,91))
    resource.setrlimit(resource.RLIMIT_AS, (1536*1024**2,1536*1024**2))
    started=time.monotonic()
    source_path=Path('reports/CLASSIFICATION_NIL22_TOP_COMPONENTS_20261107.json')
    prior_path=Path('reports/CLASSIFICATION_MODULAR_ONE_NORMAL_20261107.json')
    source=next(v for v in json.loads(source_path.read_text())['branches'] if v['name']=='A_minus_one')
    prior=next(v for v in json.loads(prior_path.read_text())['cases'] if v['branch']=='A_minus_one')
    h,x,y,z,w=xx=s.symbols('h x1 x2 x3 x4')
    nn=s.symbols('n0:15')
    env={str(v):v for v in xx+nn}

    def parse(data):
        P=s.zeros(5)
        for key,value in data.items():
            i,j=map(int,key.split(','));P[i,j]=s.sympify(value,locals=env)
        return P-P.T

    family=parse(source['universal_brackets'])
    n1,n2,n6,n7,n8,n9,u=[nn[k] for k in (1,2,6,7,8,9,10)]
    W=w-u*h;Z=z-u*y
    d=n2+n6*u;b=n6;c=n9;e=n7+n8*u+n9*u*u;f=n8+2*n9*u
    delta=-3*c/4;j=3*f;g=-3*e/2
    e0=-3*(b-delta*f)/5;c0=-3*(d-delta*e)/2
    coordinates=s.Matrix([h,x+delta*Z+e0*W+c0*h,y+delta*W,Z+j*W+g*h,W])
    A=coordinates.jacobian(xx)
    assert s.expand(A.det())==1
    canonical=s.zeros(5)
    for pair,value in zip(itertools.combinations(range(1,5),2),
                          [y*y-2*h*x,3*n1*h*h-2*x*w-2*y*z,-3*y*w,
                           -2*h*z-y*w,-3*h*w,-3*w*w]):
        canonical[pair]=value
    canonical-=canonical.T
    substitution=dict(zip(xx,coordinates))
    # This is the universal all-seven-parameter coordinate identity.
    assert all(s.expand(v)==0 for v in A*family*A.T-canonical.subs(substitution,simultaneous=True))
    at={env[k]:s.Rational(v) for k,v in source['full_source_specification']['witness'].items()}
    Ap=A.subs(at)
    inverse_sub=dict(zip(xx,Ap.inv()*s.Matrix(xx)))
    P=canonical.subs(at)
    Rold=parse(prior['normal_brackets'])
    R=(Ap*Rold*Ap.T).subs(inverse_sub,simultaneous=True).applyfunc(s.expand)
    ps=sparse(P,xx)
    assert engine.jacobi(ps)==s.zeros(350,1)
    Dcols=[]
    for ij,mon in engine.COLUMNS:
        v={key:dict(value) for key,value in ps.items()}
        v.setdefault(ij,{})[mon]=v.get(ij,{}).get(mon,0)+1
        Dcols.append(engine.jacobi(v))
    D=s.Matrix.hstack(*Dcols)
    V=[];O=[];E=s.Matrix(xx)
    for a,b in itertools.product(range(5),repeat=2):
        F=s.zeros(5);F[a,b]=1;field=F*E
        V.append(engine.vector(sparse(E*field.T-field*E.T,xx)))
        O.append(engine.vector(sparse(F*P+P*F.T-xx[b]*P.diff(xx[a]),xx)))
    G=s.Matrix.hstack(*V,*O,engine.vector(ps))
    assert D*G==s.zeros(350,G.cols)
    gp=dm(G).rref()[1]
    assert len(gp)==47
    assert len(dm(D).rref()[1])==102
    rv=engine.vector(sparse(R,xx))
    assert D*rv==s.zeros(350,1)
    weights=(0,2,1,-1,-2)
    column_weights=[weights[i]+weights[j]-sum(a*b for a,b in zip(weights,mon)) for (i,j),mon in engine.COLUMNS]
    assert all(v==0 or wt==1 for v,wt in zip(engine.vector(ps),column_weights))
    # Remove orbit directions within each weight using exact row reduction.
    quotient=[]
    for weight in sorted(set(column_weights)):
        inds=[i for i,wt in enumerate(column_weights) if wt==weight]
        piece=s.zeros(150,1)
        for i in inds:piece[i]=rv[i]
        if piece==s.zeros(150,1):continue
        # Orbit columns need not be homogeneous, but projecting is legitimate
        # because the orbit tangent is preserved by the diagonal stabilizer.
        Gw=G.copy()
        for i in range(150):
            if i not in inds:Gw[i,:]=s.zeros(1,G.cols)
        assert D*Gw==s.zeros(350,G.cols)
        rr,pivots=dm(Gw.T).rref()
        rr=rr.to_Matrix()
        for row,pivot in enumerate(pivots):piece-=piece[pivot]*rr[row,:].T
        if piece!=s.zeros(150,1):quotient.append((weight,primitive(piece)))
    assert len(quotient)==1
    weight,rv=quotient[0];R=tensor(rv,xx)
    assert D*rv==s.zeros(350,1)
    assert len(dm(G.row_join(rv)).rref()[1])==48
    J2=engine.jacobi(sparse(R,xx))
    result=dict(universal_source_identity=True,affine_determinant=1,
                coordinates=[str(v) for v in coordinates],canonical_brackets=brackets(canonical),
                same_prior_point_and_normal=True,weight=list(weights),P_weight=1,
                true_normal_weight=weight,normal_brackets=brackets(R),
                full_D_rank=102,orbit_vertical_scale_rank=47,plus_normal_rank=48)
    if J2==s.zeros(350,1):
        result.update(status='EXACT_POISSON_PENCIL',curve_brackets=brackets(P+s.Symbol('tau')*R))
    else:
        target_weight=2*weight-1
        inds=[i for i,wgt in enumerate(column_weights) if wgt==target_weight]
        Dw=D[:,inds]
        pivots=list(dm(Dw).rref()[1]);cols=[inds[i] for i in pivots]
        J=D[:,cols];rows=list(dm(J.T).rref()[1])
        correction=dm(J[rows,:]).inv().matmul(dm(J2[rows,:])).to_Matrix()
        residual=J2-J*correction
        result.update(correction_weight=target_weight,quadratic_in_homogeneous_image=residual==s.zeros(350,1))
        if residual!=s.zeros(350,1):
            result['status']='QUADRATIC_OBSTRUCTION_NO_CURVE'
        else:
            qv=s.zeros(150,1)
            for i,col in enumerate(cols):qv[col]=-correction[i]
            Q=tensor(qv,xx);tau=s.Symbol('tau')
            actual=engine.jacobi(sparse(P+tau*R+tau*tau*Q,xx))
            result.update(status='EXACT_POISSON_QUADRATIC_CURVE' if actual==s.zeros(350,1) else 'EQUIVARIANT_QUADRATIC_CURVE_UNRESOLVED',
                          second_order_brackets=brackets(Q),
                          curve_brackets=brackets(P+tau*R+tau*tau*Q),
                          residual_coefficients=sorted(set(str(s.factor(v)) for v in actual if v)))
    result.update(producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  producer_bytes=Path(__file__).read_text(),
                  exact_inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source_path,prior_path,Path(engine.__file__))},
                  seconds=round(time.monotonic()-started,4))
    out=Path('reports/CLASSIFICATION_NIL22_EQUIVARIANT_20261107.json')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','true_normal_weight','normal_brackets','seconds')}))


if __name__=='__main__':main()

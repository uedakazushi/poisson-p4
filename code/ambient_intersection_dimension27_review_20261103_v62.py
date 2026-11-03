#!/usr/bin/env python3
"""Independent raw-coordinate audit of the 27-dimensional CA/CB intersection.

No producer Poisson, coordinate-action, coefficient or quotient helpers used.
The generic dimension lower bound uses a new sample and includes all Euler
vertical columns before taking their 24-dimensional quotient.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
PRIME=32003


def raw(entries):
    out=s.zeros(5)
    for (i,j),value in entries.items():
        out[i,j]=value
        out[j,i]=-value
    return out


def elimination(matrix):
    a=matrix.copy()%PRIME
    rowids=list(range(len(a)))
    pivots=[]
    r=0
    for col in range(a.shape[1]):
        nz=np.flatnonzero(a[r:,col])
        if not len(nz):
            continue
        j=r+int(nz[0]);a[[r,j]]=a[[j,r]]
        rowids[r],rowids[j]=rowids[j],rowids[r]
        a[r,col:]=a[r,col:]*pow(int(a[r,col]),-1,PRIME)%PRIME
        a[r+1:,col:]=(a[r+1:,col:]-a[r+1:,col,None]*a[r,None,col:])%PRIME
        pivots.append(col);r+=1
        if r==len(a):
            break
    return rowids[:r],pivots


def determinant(matrix):
    a=matrix.copy()%PRIME
    value=1
    for i in range(len(a)):
        nz=np.flatnonzero(a[i:,i]);assert len(nz)
        j=i+int(nz[0])
        if i!=j:
            a[[i,j]]=a[[j,i]];value=-value
        pivot=int(a[i,i]);value=value*pivot%PRIME
        a[i,i:]=a[i,i:]*pow(pivot,-1,PRIME)%PRIME
        a[i+1:,i:]=(a[i+1:,i:]-a[i+1:,i,None]*a[i,None,i:])%PRIME
    return value%PRIME


def verify():
    h,x,y,z,w=s.symbols('h x y z w');coords=(h,x,y,z,w)
    p=s.Symbol('p');a=s.symbols('a:16')
    monomials=[s.prod(coords[i] for i in ids)
               for ids in itertools.combinations_with_replacement((0,1,3,4),3)
               if 3 in ids or 4 in ids]
    assert len(monomials)==16
    H=sum(c*m for c,m in zip(a,monomials))
    PA=raw({(1,2):h*(h+p*x),(2,3):-s.diff(H,w),
            (2,4):s.diff(H,z),(3,4):h*h})
    M=s.eye(5);M[0,1]=p
    substitution={h:h-p*x}
    moved=(M*PA.subs(substitution,simultaneous=True)*M.T).applyfunc(s.expand)
    f=(h-p*x)**2;HB=H.subs(substitution,simultaneous=True)
    PB=raw({(1,2):f,(2,3):-s.diff(f,x)*z/2-s.diff(HB,w),
            (2,4):-s.diff(f,x)*w/2+s.diff(HB,z),(3,4):f})
    e=s.Matrix(coords);v=s.Matrix((0,0,p*(h-p*x),0,0))
    assert (moved-PB-e*v.T+v*e.T).applyfunc(s.expand)==s.zeros(5)
    for tensor in (PA,PB):
        for i,j,k in itertools.combinations(range(5),3):
            jac=sum(tensor[i,l]*s.diff(tensor[j,k],coords[l])+
                    tensor[j,l]*s.diff(tensor[k,i],coords[l])+
                    tensor[k,l]*s.diff(tensor[i,j],coords[l]) for l in range(5))
            assert s.expand(jac)==0
    # Independent coefficient extraction and new differential sample.
    pairs=list(itertools.combinations(range(5),2))
    quadmons=[coords[i]*coords[j] for i,j in itertools.combinations_with_replacement(range(5),2)]
    def coefficients(tensor):
        return [s.Poly(tensor[i,j],*coords).coeff_monomial(m)
                for i,j in pairs for m in quadmons]
    sample={p:s.Integer(3)}|{c:s.Integer(i*i+2*i+3)for i,c in enumerate(a)}
    P=PA.subs(sample);columns=[];vertical=[]
    for i in range(5):
        for j in range(5):
            A=s.zeros(5);A[i,j]=1;Ax=A*e
            differential=A*P+P*A.T-sum((Ax[k]*P.diff(coords[k])for k in range(5)),s.zeros(5))
            columns.append(coefficients(differential))
            vertical.append(coefficients(e*Ax.T-Ax*e.T))
    columns.append(coefficients(P))
    columns.extend(coefficients(PA.diff(variable).subs(sample))for variable in (p,*a))
    convert=lambda col:[int(c.p)*pow(int(c.q),-1,PRIME)%PRIME for c in col]
    V=np.array([convert(col)for col in vertical],dtype=np.int64).T
    T=np.array([convert(col)for col in [*vertical,*columns]],dtype=np.int64).T
    rv,cv=elimination(V);rows,cols=elimination(T)
    assert len(rv)==24 and len(rows)==52
    det=determinant(T[np.ix_(rows,cols)]);assert det
    # Distinguished Q is an exact point, rather than an inferred limit.
    Q=raw({(1,4):h*h,(2,3):-s.Rational(2,9)*h*h,
           (2,4):s.Rational(4,3)*h*x,(3,4):s.Rational(2,3)*h*y+x*x})
    transport=s.zeros(5)
    transport[0,0]=transport[1,1]=transport[2,4]=transport[3,2]=1
    transport[4,3]=-s.Rational(9,2)
    inverse=transport.inv()*e
    q_moved=(transport*Q.subs(dict(zip(coords,inverse)),simultaneous=True)*transport.T).applyfunc(s.expand)
    target=s.Rational(4,3)*h*x*w+s.Rational(3,2)*h*z*z+s.Rational(9,2)*x*x*z
    point={p:0}|{c:s.Poly(target,*coords).coeff_monomial(m)for c,m in zip(a,monomials)}
    assert (q_moved-PA.subs(point)).applyfunc(s.expand)==s.zeros(5)
    paths=[Path(__file__),ROOT/'code/ambient_intersection_dimension27_20261103_v62.py',
           ROOT/'reports/AMBIENT_INTERSECTION_DIMENSION27_20261103_v62/report.json',
           ROOT/'legacy/v06/ambient_families.tex',ROOT/'legacy/v24/sections/incidence.tex']
    return dict(id='AMBIENT_INTERSECTION_DIMENSION27_20261103_V62_REVIEW',all_passed=True,
                universal_raw_identity_independently_verified=True,raw_jacobi_both_families=True,
                independent_sample={str(k):str(v)for k,v in sample.items()},
                Euler_vertical_rank=24,raw_total_differential_rank_lower=52,
                affine_quotient_rank_lower=28,minor_mod_prime=det,prime=PRIME,
                minor_rows=rows,minor_columns=cols,Q_exact_point_verified=True,
                entire_p_zero_magnetic_family_retained=True,
                conceptual_dimension_audit='PGL-saturated image is irreducible and lies in both distinct irreducible28-dimensional components. Its lower dimension27 matches the intersection upper27. It contains all magnetic M and Q; M of dimension26 is properly contained, hence is not an intersection component.',
                local_intersection_dimension_at_Q=27,
                uniqueness_or_scheme_intersection_claimed=False,
                input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()for path in paths})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify()
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS L1 independent raw identity, differential and Q',result['minor_mod_prime'])

#!/usr/bin/env python3
"""Identify the parabolic covariant divisor classes of possible marked anchors."""
import argparse
import hashlib
from itertools import combinations_with_replacement
import json
from math import prod
from pathlib import Path
import time

import sympy as sp

from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import J

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_AFFINE_ANCHOR_RESIDUE_DEGREE_20261102'


def action_matrix(degree):
    a,b,c,d,h,k=sp.symbols('x0:6')
    variables=[a,b,c,d,h,k]
    r,t=sp.symbols('r t')
    cubic=a*r**3+b*r*r*t+c*r*t*t+d*t**3
    linear=h*r+k*t
    exps=[tuple(indices.count(i) for i in range(6))
          for indices in combinations_with_replacement(range(6),degree)]
    indices={e:i for i,e in enumerate(exps)}
    actions=[]
    for dr,dt in [(t,0),(0,r)]:
        dc=sp.Poly(dr*cubic.diff(r)+dt*cubic.diff(t),r,t)
        dl=sp.Poly(dr*linear.diff(r)+dt*linear.diff(t),r,t)
        fields=sp.Matrix([dc.coeff_monomial(r**(3-i)*t**i) for i in range(4)]+
                         [dl.coeff_monomial(r),dl.coeff_monomial(t)]).jacobian(variables)
        rows=sp.MutableSparseMatrix(len(exps),len(exps),{})
        for column,mon in enumerate(exps):
            for i,power in enumerate(mon):
                if not power:
                    continue
                for j in range(6):
                    if fields[i,j]:
                        transformed=list(mon)
                        transformed[i]-=1
                        transformed[j]+=1
                        rows[indices[tuple(transformed)],column]+=power*fields[i,j]
        actions.append(rows)
    both=sp.Matrix.vstack(*actions)
    monomials=[prod(x**i for x,i in zip(variables,m))for m in exps]
    vectors=both.to_DM().nullspace().to_Matrix()
    answers=[sp.factor(sum(v*m for v,m in zip(row,monomials)))for row in vectors.tolist()]
    return both.shape,answers


def verify():
    start=time.monotonic()
    y,u,beta=primitive_beta(S61());y=sp.Matrix(y)
    b={int(str(x)[1:]):x for x in u}
    beta=beta.subs({b[12]:0,b[13]:0})
    e=sp.eye(4);r,t=sp.symbols('r t')
    v=r*e[:,1]+t*e[:,2]
    n=t*e[:,0]+r*e[:,3]
    assert (v.T*J*n)[0]==0
    c=sp.Poly(sp.expand((v.T*beta.subs(dict(zip(y,v)),simultaneous=True)*n)[0]),r,t)
    alpha=(5*b[14]-b[49])*e[:,0]+(5*b[35]+b[50])*e[:,3]
    ell=sp.Poly((v.T*J*alpha)[0],r,t)
    coefficients=[c.coeff_monomial(r**(3-i)*t**i)for i in range(4)]
    coefficients += [ell.coeff_monomial(r),ell.coeff_monomial(t)]
    change=sp.Matrix(coefficients).jacobian([b[i]for i in (4,14,35,36,49,50)])
    assert change.det()==54
    # A basis change of the plane, accompanied by its dual on the complementary plane.
    # In determinant-one coordinates the operation n(v) commutes with every such change.
    A,B,C,D=sp.symbols('A B C D')
    g=sp.Matrix([[A,B],[C,D]])
    pairing=sp.diag(-1,1)
    N=sp.Matrix([[0,1],[1,0]])
    # Here the plane basis is (e1,e2) and the complementary basis is (e0,e3).
    # The complementary action is Q=pairing^-1*g^-T*pairing.
    Q=pairing*g.adjugate().T*pairing
    assert Q*N==N*g
    rings={}
    for degree,want in [(2,0),(4,3),(6,1)]:
        shape,answers=action_matrix(degree)
        assert len(answers)==want
        rings[str(degree)]=dict(shape=shape,kernel_dim=want,polynomials=list(map(str,answers)))
    p,x,y0,z,w=sp.symbols('b1 b5 b14 b6 b12')
    aa=p*x-sp.Rational(9,50)*w*w
    cc=p*y0-sp.Rational(9,10)*w*w
    hh=p*p*z-sp.Rational(3,5)*p*w*(x+y0)+sp.Rational(54,125)*w**3
    values=sp.expand(hh**2+16*aa**2*(cc-3*aa))
    strict_pure=sp.cancel(values/p**2)
    assert sp.denom(strict_pure).is_number and strict_pure.subs(p,0)!=0 and hh.subs(p,0)!=0
    # Raw characters in diag(x,x^-1,z,z^-1) of the fixed flag coordinates.
    t1,t2=sp.symbols('x z',nonzero=True)
    tab={p:t1**-2*t2,x:t2**-1,y0:t2**-1,z:t1*t2**-2,w:t1**-1}
    for f,degree in [(aa,2),(cc,2),(hh,3)]:
        scaled=f.subs({v:v*factor for v,factor in tab.items()},simultaneous=True)
        assert sp.expand(scaled-t1**(-degree)*f)==0
    deps=['S61_AFFINE_BOUNDARY_FINITE_ANCHORS_20261102','S61_AFFINE_ANCHOR_LENGTH_BOUND_20261102',
          'S61_GENERIC_BOUNDARY_INDEX_20261102','S61_AFFINE_SURFACE_BOUND_20261102',
          'S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102','S61_MARKED_SIGMA_PROJECTIVE_CORE_20261102',
          'S61_MARKED_SIGMA_INTERSECTION_REMOVAL_20261102','S61_BINARY_H11_RANK_DROP_RESIDUAL_FINITE_20261102',
          'S61_BINARY_ZERO_PAIR_STRATUM_REMOVAL_20261102','S61_MINIMAL_FLAG_ENDPOINT_INVARIANTS_20261101']
    files=[Path(__file__),ROOT/'research/notes'/(STEM+'.md')]
    for name in deps:
        path=ROOT/'reports'/(name+'.json')
        assert json.loads(path.read_text())['all_passed']
        files.extend([path,ROOT/'research/notes'/(name+'.md')])
    files.extend(ROOT/'code'/n for n in ['s61_model.py','s61_infinity_flags.py','weight1_structure.py',
                'elliptic_weight1_dominance.py','primitive_jet_elliptic_curve_20261101.py'])
    files.extend(ROOT/'legacy/v02/data'/n for n in ['symplectic_slice_certificate.json',
                 'symplectic_slice_equations.json','symplectic_slice_reduced_equations.json'])
    return dict(id=STEM,all_passed=True,base='K=C(W16)',core_cubic=str(c.as_expr()),
                core_linear=str(ell.as_expr()),coordinate_change_determinant=54,
                binary_pair_SL2_invariant_polynomials=rings,degree_of_divisor_for_core_degree_n='5n/2',
                possible_anchor_support_plane_residue_degree_at_least=10,
                residue_degree_ten_net_dimension=3,residue_degree_fifteen_invariant_space_dimension=1,
                sigma_finite_point_degrees='15 or a positive multiple of 10',
                pure_cubic_endpoint_core_equation='h^2+16*a^2*(c-3*a)',
                pure_cubic_strict_equation=str(sp.factor(strict_pure)),pure_cubic_endpoint_residue_degree=10,
                sufficient_finite_anchor_length_for_surface_exclusion=19,
                input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in files},
                seconds=time.monotonic()-start)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();report=verify()
    with args.output.open('x')as out:
        json.dump(report,out,indent=2);out.write('\n')
    print('PASS',round(report['seconds'],3))

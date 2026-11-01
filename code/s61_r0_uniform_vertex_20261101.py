#!/usr/bin/env python3
"""Uniform QQ[k] radical for the R0 vertex model and symplectic translation.

The relative radical is certified before specializing k. The normal-form
translation applies when (3k-1)(3k+4) is nonzero. Other cases are not asserted.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sympy as s
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import COORDS,J

ROOT=Path(__file__).resolve().parents[1]


def maximal_faces(edges):
    faces=[]
    for n in range(1,10):
        for values in combinations(range(9),n):
            face=set(values)
            if not any(edge<=face for edge in edges):faces.append(face)
    return [face for face in faces if not any(face<other for other in faces)]


def run():
    source=ROOT/'reports/S61_R0_VERTEX_RADICAL_20261101.json'
    old=json.loads(source.read_text());k=s.Symbol('k');h=list(map(s.Symbol,old['variables']))
    F=list(map(s.sympify,old['nuisance_free_equations']))
    monomials=list(map(s.sympify,old['radical_squarefree_generators']))
    # This finite list is proved by the universal ideal memberships below.
    # It is not inferred from testing those specializations.
    special={s.Rational(-1,2):[(0,6),(1,7)],
        s.Rational(1,3):[(0,3),(0,7),(2,6),(3,7)],
        s.Rational(-4,3):[(1,3),(1,6),(3,6),(4,7)],
        s.Rational(2):[(0,4),(2,7),(3,8)],
        s.Rational(-3):[(1,2),(3,5),(4,6)],
        s.Rational(11,3):[(2,8)],s.Rational(-14,3):[(4,5)],
        s.Rational(7,6):[(0,8)],s.Rational(-13,6):[(1,5)]}
    pairs=[tuple(i for i,e in enumerate(s.Poly(m,*h).monoms()[0]) if e) for m in monomials]
    coefficients=[s.sympify(s.prod(k-value for value,allowed in special.items() if pair in allowed))
                  for pair in pairs]
    generators=[s.expand(c*m) for c,m in zip(coefficients,monomials,strict=True)]
    for f in F:
        for ex,co in s.Poly(f,*h).terms():
            if co:
                i=next(j for j,m in enumerate(monomials) if s.Poly(m,*h).monoms()[0]==ex)
                assert s.rem(co,coefficients[i],k)==0
    basis=s.groebner(F,*h,k,order='grevlex',domain=s.QQ)
    assert len(basis.polys)==125
    powers=[];quotients=[]
    for generator in generators:
        for power in range(1,6):
            q,remainder=basis.reduce(generator**power)
            if remainder==0:
                assert s.expand(sum(a*b.as_expr() for a,b in zip(q,basis.polys,strict=True))
                                -generator**power)==0
                powers.append(power);quotients.append(list(map(str,q)));break
        else:raise AssertionError(f'No membership for {generator}')
    assert max(powers)==5
    generic_faces=[set(v) for v in old['maximal_coordinate_faces']]
    fibers={value:maximal_faces([set(pair) for pair,c in zip(pairs,coefficients,strict=True)
                                 if c.subs(k,value)!=0]) for value in special}
    vertical={value:[face for face in faces if not any(face<=generic for generic in generic_faces)]
              for value,faces in fibers.items()}
    # Coefficient-by-coefficient proof of intersection with22 explicit primes.
    # It suffices to inspect squarefree supports, since all h-generators are squarefree.
    for n in range(10):
        for values in combinations(range(9),n):
            support=set(values)
            cs=[c for pair,c in zip(pairs,coefficients,strict=True) if set(pair)<=support]
            actual=s.S.Zero
            for c in cs:actual=s.gcd(actual,c)
            if actual!=0:actual=s.Poly(actual,k).monic().as_expr()
            expected=(s.S.Zero if any(support<=face for face in generic_faces) else
                s.prod(k-value for value,faces in vertical.items() if any(support<=face for face in faces)))
            assert s.expand(actual-expected)==0
    assert len(generic_faces)+sum(map(len,vertical.values()))==22

    # The actual symplectic coordinate change removes K,L uniformly away from resonance.
    zz,u,beta=primitive_beta(S61());allu=s.symbols('u1:62');K,L=s.symbols('K L')
    sub={allu[3]:0,allu[11]:0,allu[4]:0,allu[12]:0,allu[13]:k,
         allu[48]:2*k+1,allu[25]:0,allu[34]:0,allu[49]:K,allu[50]:-L}
    A=-K/(3*k-1);B=-L/(3*k+4)
    T=s.Matrix([[1,0,0,0],[0,1,A,-B],[B,0,1,0],[A,0,0,1]])
    assert T.T*J*T==J
    inverse=T.inv();bb=beta.subs(sub)
    pulled=(inverse.T*bb.subs(dict(zip(zz,inverse*s.Matrix(zz))),simultaneous=True)*inverse).applyfunc(s.expand)
    vec=s.Matrix([s.diff(beta[i,j],zz[t]) for (i,j),t in COORDS]);bbasis=vec.jacobian(u)
    rows=list(bbasis.T.rref()[1]);left=bbasis[rows,:].inv()
    newvec=s.Matrix([s.diff(pulled[i,j],zz[t]) for (i,j),t in COORDS])
    newu=(left*newvec[rows,:]).applyfunc(s.cancel)
    assert all(newu[u.index(allu[i])]==0 for i in [3,11,4,12])
    assert newu[u.index(allu[13])]==k and newu[u.index(allu[48])]==2*k+1
    assert s.cancel(5*newu[u.index(allu[34])]+newu[u.index(allu[49])])==0
    assert s.cancel(5*newu[u.index(allu[25])]-newu[u.index(allu[50])])==0

    return {'id':'S61_R0_UNIFORM_VERTEX_20261101','all_passed':True,'field':'QQ[k,h0,...,h8]',
        'variables':list(map(str,[*h,k])),'original_equations':list(map(str,F)),
        'radical_generators':list(map(str,generators)),
        'groebner_basis':list(map(lambda p:str(p.as_expr()),basis.polys)),
        'generator_powers_in_original_ideal':powers,
        'power_membership_quotients':quotients,
        'generic_coordinate_faces':[sorted(face) for face in generic_faces],
        'special_fiber_maximal_faces':{str(v):list(map(sorted,faces)) for v,faces in fibers.items()},
        'extra_vertical_prime_faces':{str(v):list(map(sorted,faces)) for v,faces in vertical.items()},
        'relative_radical_prime_count':22,
        'all_512_support_intersection_identities':True,
        'symplectic_translation':[[str(v) for v in row] for row in T.tolist()],
        'transformed_primitive_u':{str(a):str(b) for a,b in zip(u,newu,strict=True)},
        'uniform_translation_open':'(3*k-1)*(3*k+4)!=0',
        'scope':'Relative radical of the nuisance-free fixed-factor vertex model, uniformly in k. '
                'Symplectic translation covers arbitrary K,L off the two resonances. '
                'The22 relative coordinate primes are not global S61 components.',
        'proof_note':'research/notes/S61_R0_UNIFORM_VERTEX_20261101.md',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in [Path(__file__),source,ROOT/'code/primitive_jet_elliptic_curve_20261101.py']}}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();report=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'groebner_size':125,'max_power':5,
                      'relative_primes':22,'output':str(args.output)}))


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Generic R0 vertex fiber: an exact radical over QQ(k,K,L).

No specialization of k is used for the ideal calculation. The result is only
on the generic nondegenerate vertex stratum, not on its exceptional parameters.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def sm(matrix):
    return [[str(value) for value in row] for row in matrix.tolist()]


def run():
    source=ROOT/'reports/S61_R0_INFINITY_20261101.json'
    data=json.loads(source.read_text())
    g=list(map(s.Symbol,data['g_variables']))
    raw=list(map(s.sympify,data['leading_equations_70_through_133']))
    raw.append(s.sympify(data['leading_G6']))
    u=s.symbols('u1:62'); kk=s.symbols('k0:5')
    sub={u[3]:0,u[11]:0,u[4]:0,u[12]:0,u[13]:kk[0],u[48]:2*kk[0]+1,
         u[25]:kk[1],u[34]:kk[2],u[49]:kk[3],u[50]:kk[4]}
    F=[s.expand(f.subs(sub)) for f in raw]
    assert all(s.diff(f,kk[2])-5*s.diff(f,kk[3])==0 and
               s.diff(f,kk[1])+5*s.diff(f,kk[4])==0 for f in F)
    k,K,L=s.symbols('k K L')
    F=[s.expand(f.subs(dict(zip(kk,[k,0,0,K,-L],strict=True)))) for f in F]
    A=-K/(3*k-1);B=-L/(3*k+4)
    h=s.symbols('h0:9');x,z,w=s.symbols('x z w')
    cubics=[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3]
    G=sum(a*b.subs({z:z+B*x,w:w+A*x},simultaneous=True)
          for a,b in zip(h,cubics,strict=True))
    pG=s.Poly(G,x,z,w)
    substitution={v:pG.coeff_monomial(mon) for v,mon in zip(g,cubics,strict=True)}
    change=s.Matrix([substitution[v] for v in g]).jacobian(h)
    assert change.det()==1
    quadratics=[a*b for a,b in combinations_with_replacement(h,2)]
    translated=[s.Poly(s.expand(f.subs(substitution,simultaneous=True)),*h) for f in F]
    M=s.Matrix([[s.cancel(p.coeff_monomial(mon)) for mon in quadratics] for p in translated])
    M0=M.subs({K:0,L:0})
    D0,pivots=M0.to_DM().rref();R=D0.to_Matrix();rank=len(pivots)
    assert rank==15
    assert (M-M[:,list(pivots)]*R[:rank,:]).applyfunc(s.cancel)==s.zeros(65,45)
    rows=list(M0[:,list(pivots)].T.to_DM().rref()[1])
    rank_minor=M0.extract(rows,list(pivots))
    rank_det=s.factor(rank_minor.det())
    assert rank_det!=0
    # The preceding identity proves inclusion in the nuisance-free rowspace.
    # Its rank is at least15 because the displayed minor specializes to rank_det.
    zero_F=[s.expand(f.subs({K:0,L:0}).subs(dict(zip(g,h,strict=True)))) for f in F]
    ideal=s.groebner(zero_F,*h,order='grevlex',domain=s.QQ.frac_field(k))
    assert len(ideal.polys)==46
    faces=[{0,2,5},{1,4,8},{0,1},{3},{6},{7}]
    pairs=[(i,j) for i in range(9) for j in range(i+1,9)
           if not any({i,j}.issubset(face) for face in faces)]
    assert len(pairs)==29
    # I is contained in the squarefree monomial ideal J.
    for f in zero_F:
        for mon,coef in s.Poly(f,*h).terms():
            if coef:
                assert any(mon[i] and mon[j] for i,j in pairs)
    powers=[];quotients=[]
    for i,j in pairs:
        mon=h[i]*h[j]
        for power in range(1,4):
            q,remainder=ideal.reduce(mon**power)
            if remainder==0:
                assert s.cancel(sum(a*b.as_expr() for a,b in zip(q,ideal.polys,strict=True))
                                -mon**power)==0
                powers.append(power);quotients.append(list(map(str,q)));break
        else:
            raise AssertionError(f'No certified nilpotence power for {mon}')
    assert max(powers)==3
    return {'id':'S61_R0_VERTEX_RADICAL_20261101','all_passed':True,
        'field':'QQ(k,K,L)','generic_parameters_only':True,
        'parameters':{'k':'k0','K':'5*k2+k3','L':'5*k1-k4'},
        'coordinate_change':'Z=z-L*x/(3*k+4), W=w-K*x/(3*k-1)',
        'g_in_h':{str(a):str(b) for a,b in substitution.items()},
        'coordinate_change_determinant':'1','variables':list(map(str,h)),
        'quadratic_monomials':list(map(str,quadratics)),
        'translated_equation_matrix':sm(M),
        'nuisance_free_rref':sm(R[:rank,:]),'pivot_columns':list(pivots),
        'rank_minor_rows':rows,'rank_minor_determinant':str(rank_det),
        'all_65_rowspace_identities':True,'quadratic_rank':15,
        'nuisance_free_equations':list(map(str,zero_F)),
        'groebner_order':'grevlex','groebner_basis':list(map(lambda p:str(p.as_expr()),ideal.polys)),
        'radical_squarefree_generators':[str(h[i]*h[j]) for i,j in pairs],
        'radical_generator_powers_in_ideal':powers,
        'power_membership_quotients_in_groebner_basis':quotients,
        'maximal_coordinate_faces':[sorted(face) for face in faces],
        'projective_components':['P2: h0,h2,h5','P2: h1,h4,h8','P1: h0,h1',
                                'point h3','point h6','point h7'],
        'nonbinary_cubics_mod_x3':['x*Z*W','Z^2*W','Z*W^2'],
        'limitations':['Generic vertex stratum only; exceptional k and degenerate conics remain.',
            'Six fixed-v pieces are not a count of global irreducible components.',
            'Unsaturated infinity only; no affine-fiber exhaustion.'],
        'proof_note':'research/notes/S61_R0_VERTEX_RADICAL_20261101.md',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in [Path(__file__),source]}}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'quadratic_rank':15,'radical_generators':29,
                      'max_power':3,'output':str(args.output)}))


if __name__=='__main__':main()

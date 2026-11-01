#!/usr/bin/env python3
"""Exact primitive-jet section of U*(1) and its five-Pfaffian curve.

The generic smooth elliptic-quintic assertion is proved in the accompanying
note by global generation, Chern classes and Koszul cohomology. This program
certifies the linear construction and its S61 coordinate conventions only.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_infinity_flags import old_beta_basis
from weight1_structure import J,minor_witness
ROOT=Path(__file__).resolve().parents[1]


def primitive_beta(model):
    z,basis=old_beta_basis(model)
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    u=[s.Symbol('u'+str(i+1)) for i in one]
    old=sum((a*b for a,b in zip(u,basis,strict=True)),s.zeros(4))
    linear=-J*old*J;tau=linear[0,1]+linear[2,3]
    alpha=s.Matrix([s.diff(tau,v)/5 for v in z]);v=J*alpha;zz=s.Matrix(z)
    normalized=linear-2*alpha.dot(zz)*J+zz*v.T-v*zz.T
    beta=(-J*normalized*J).applyfunc(s.expand)
    assert beta[0,1]+beta[2,3]==0
    for i in range(4):
        for j in range(i+1,4):
            for k in range(j+1,4):
                assert s.expand(s.diff(beta[j,k],z[i])-s.diff(beta[i,k],z[j])+s.diff(beta[i,j],z[k]))==0
    return z,u,beta


def pfaffian4(matrix):
    assert matrix.shape==(4,4)
    return s.expand(matrix[0,1]*matrix[2,3]-matrix[0,2]*matrix[1,3]+matrix[0,3]*matrix[1,2])


def curve_model(beta,z):
    p=s.symbols('p12 p13 p14 p23 p24')
    a,b,c,d,e=p
    plucker=s.Matrix([[0,a,b,c],[-a,0,d,e],[-b,-d,0,-a],[-c,-e,a,0]])
    dual=s.Matrix([[0,-a,-e,d],[a,0,c,-b],[e,-c,0,a],[-d,b,-a,0]])
    quadric=pfaffian4(plucker)
    assert plucker*dual==-quadric*s.eye(4)
    covector=s.Matrix([s.expand(sum(s.diff(beta[i,j],z[k])*plucker[i,j]
                                   for i in range(4) for j in range(i+1,4))) for k in range(4)])
    matrix=dual.row_join(covector).col_join((-covector.T).row_join(s.zeros(1)))
    signed=s.Matrix([(-1)**i*pfaffian4(matrix.minor_submatrix(i,i)) for i in range(5)])
    expected=(plucker*covector).col_join(s.Matrix([quadric]))
    assert (signed-expected).applyfunc(s.expand)==s.zeros(5,1)
    assert (matrix*signed).applyfunc(s.expand)==s.zeros(5,1)
    return p,plucker,covector,matrix,signed


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    model=S61();z,u,beta=primitive_beta(model);p,plucker,a,matrix,pf=curve_model(beta,z)
    v=s.Matrix([0,-1,0,0]);w=s.Matrix([0,0,0,-1])
    form=s.expand((v.T*beta*w)[0]);restriction=s.Matrix([form.subs(dict(zip(z,v))),form.subs(dict(zip(z,w)))])
    evaluation=restriction.jacobian(u);assert evaluation.rank()==2
    # A section vanishes identically iff the four quadrics Pa vanish on Q3.
    # Reduce p12^2 modulo Q=-p12^2-p13*p24+p14*p23 and compute the
    # coefficient matrix of the resulting 4x14 coefficient linear forms.
    remainder=[s.rem(s.Poly(f,*p),s.Poly(pf[4],*p)).as_expr() for f in pf[:4,0]]
    monomials=sorted(set().union(*(set(s.Poly(f,*p).monoms()) for f in remainder)))
    coefficient=s.Matrix([s.Poly(f,*p).coeff_monomial(mon) for f in remainder for mon in monomials]).jacobian(u)
    witness=minor_witness(coefficient);assert witness['rank']==16
    report={'id':'PRIMITIVE_JET_ELLIPTIC_CURVE_20261101','field':'QQ','all_passed':True,
        'primitive_jet_variables':list(map(str,u)),'affine_variables':list(map(str,z)),
        'primitive_beta':[[str(v) for v in row] for row in beta.tolist()],
        'plucker_variables':list(map(str,p)),'plucker_relation':str(pf[4]),
        'plucker_matrix':[[str(v) for v in row] for row in plucker.tolist()],
        'contraction_covector':list(map(str,a)),
        'pfaffian_matrix':[[str(v) for v in row] for row in matrix.tolist()],
        'signed_pfaffians':list(map(str,pf)),
        'section_injectivity':{'matrix':[[str(v) for v in row] for row in coefficient.tolist()],
                               'coefficient_monomials':[list(m) for m in monomials],
                               'witness':witness},
        'representative_evaluation_rank':2,'representative_plane_restriction':list(map(str,restriction)),
        'representative_flag_covector':str(form),
        'identities':['P*dual(P)=-Pf(P)*I4','signed maximal Pfaffians(M)=(P*a,Pf(P))','M*Pfaffians=0'],
        'scope':'Constructs the intrinsic elliptic-quintic curve of a sufficiently general primitive jet; does not reconstruct a Poisson extension or prove exhaustion.',
        'proof_note':'research/notes/PRIMITIVE_JET_ELLIPTIC_CURVE_20261101.md',
        'source_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in [Path(__file__),ROOT/'code/s61_infinity_flags.py',ROOT/'code/s61_model.py']}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'section_injectivity_rank':16,'representative_evaluation_rank':2,'output':str(args.output)}))

if __name__=='__main__':main()

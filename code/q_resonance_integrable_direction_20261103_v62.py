#!/usr/bin/env python3
"""An exact integrable extra tangent at kappa=1/4, with its weight-block cone.

Constructs the actual canonical linearization directly in raw coordinates.
This is not the complete 35-variable formal local model.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
X=s.symbols('x0:5')
PAIRS=list(itertools.combinations(range(5),2))
MON2=[tuple(ids)for ids in itertools.combinations_with_replacement(range(5),2)]
MON3=[s.prod(X[i]for i in ids)for ids in itertools.combinations_with_replacement(range(5),3)]


def tensor(entries):
    P=s.zeros(5)
    for (i,j),value in entries.items():
        P[i,j]=value;P[j,i]=-value
    return P


def jacobi(P,Q=None):
    out=[]
    for i,j,k in itertools.combinations(range(5),3):
        value=0
        for u,v,w in ((i,j,k),(j,k,i),(k,i,j)):
            for l in range(5):
                value+=P[u,l]*s.diff((P if Q is None else Q)[v,w],X[l])
                if Q is not None:
                    value+=Q[u,l]*s.diff(P[v,w],X[l])
        poly=s.Poly(s.expand(value),*X)
        out.extend(poly.coeff_monomial(mon)for mon in MON3)
    return s.Matrix(out)


def divergence(P):
    return s.Matrix([s.expand(sum(s.diff(P[i,j],X[i])for i in range(5)))for j in range(5)])


def linear(P,R):
    div=divergence(R)
    first=[s.Poly(v,*X).coeff_monomial(x)for v in div for x in X]
    return s.Matrix(first+list(jacobi(P,R)))


def verify():
    started=time.monotonic();a=s.Symbol('a');eps=s.Symbol('eps')
    x0,x1,x2,x3,x4=X;astar=s.Rational(2,9)
    Q=tensor({(1,4):a*x0*x0,(2,3):-s.Rational(2,9)*x0*x0,
              (2,4):s.Rational(4,3)*x0*x1,(3,4):s.Rational(2,3)*x0*x2+x1*x1})
    basis=[];indices=[]
    for n,((i,j),mon)in enumerate(itertools.product(PAIRS,MON2)):
        if i+j-sum(mon)==-1:
            basis.append(tensor({(i,j):s.prod(X[k]for k in mon)}));indices.append(n)
    assert len(basis)==20
    full=s.Matrix.hstack(*(linear(Q,R)for R in basis))
    rows=[i for i in range(375)if any(full.row(i))];A=full[rows,:]
    generic=DomainMatrix.from_Matrix(A).convert_to(s.QQ.frac_field(a))
    grref,gpiv=generic.rref();assert len(gpiv)==17
    A0=A.subs(a,astar);K=A0.nullspace();assert len(K)==4 and A0.rank()==16
    Rs=[sum((v*B for v,B in zip(k,basis)),s.zeros(5))for k in K]
    R=tensor({(0,1):4*x0*x2-3*x1*x1,(0,2):4*x0*x3+12*x1*x2,
              (0,3):4*x0*x4+12*x1*x3+s.Rational(15,2)*x2*x2,
              (0,4):-18*x1*x4,(1,2):-s.Rational(20,3)*x0*x4+4*x1*x3+x2*x2,
              (1,3):-6*x1*x4+6*x2*x3,(1,4):-4*x2*x4,
              (2,3):4*x2*x4+6*x3*x3,(2,4):-4*x3*x4,(3,4):x4*x4})
    r=sum((v*k for v,k in zip((4,6,-4,1),K)),s.zeros(20,1))
    assert sum((v*B for v,B in zip(r,basis)),s.zeros(5))==R
    derivative=A.diff(a);tau=derivative*r
    witness=next(v for v in A0.T.nullspace()if (v.T*tau)[0]!=0)
    assert A0.row_join(tau).rank()==17
    late=[s.factor((witness.T*derivative*k)[0])for k in K]
    # Full second-order obstruction for a pure weight -1 first tangent.
    z=s.symbols('z:4');Rz=sum((v*B for v,B in zip(z,Rs)),s.zeros(5))
    T=tensor({(0,1):x4*x4})
    B=linear(Q.subs(a,astar),T)
    pivot=next(i for i,v in enumerate(B)if v)
    quadratic=s.Matrix([0]*25+list(jacobi(Rz)))
    obstruct=[s.expand(quadratic[i]-B[i]/B[pivot]*quadratic[pivot])for i in range(375)]
    obstruct=[v for v in obstruct if v]
    G=s.groebner(obstruct,*z,order='grevlex')
    w0,w1,w2,t=s.symbols('w0 w1 w2 t')
    changed=[s.Poly(v.subs(dict(zip(z,(w0+4*t,w1+6*t,w2-4*t,t)))),w0,w1,w2,t)for v in obstruct]
    mons=[w0*w0,w0*w1,w0*w2,w1*w1,w1*w2,w2*w2,w0*t,w1*t,w2*t]
    assert all(poly.coeff_monomial(t*t)==0 for poly in changed)
    assert s.Matrix([[poly.coeff_monomial(mon)for mon in mons]for poly in changed]).rank()==9
    # The full nonlinear curve, to every order, not only an obstruction class.
    P=Q.subs(a,astar)+eps*R-75*eps**2*T
    assert divergence(P)==s.zeros(5,1)and jacobi(P)==s.zeros(350,1)
    assert linear(Q.subs(a,astar),R)==s.zeros(375,1)
    F=0
    for i in range(5):
        j,k,l,m=[n for n in range(5)if n!=i]
        F+=(-1)**i*X[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])
    F=s.expand(F);assert s.Poly(F,*X).coeff_monomial(x0**5)==-s.Rational(4,81)
    paths=[Path(__file__),ROOT/'legacy/v24/sections/obstructions.tex',
           ROOT/'legacy/v24/proofs/obstructions/parameter_rank.json']
    return dict(id='Q_RESONANCE_INTEGRABLE_DIRECTION_20261103_V62',all_passed=True,
                a=str(astar),kappa='1/4',weight_convention='pair-index-sum minus coefficient monomial-index-sum',
                first_tangent_weight=-1,weight_block_columns=indices,active_rows=rows,
                generic_block_rank=17,resonant_block_rank=16,extra_not_generic_kernel_limit=True,
                derivative_cokernel_left_witness=list(map(str,witness)),
                derivative_cokernel_nonzero_value=str((witness.T*tau)[0]),
                late_extension_functional=list(map(str,late)),
                four_kernel_tensors=[{f'{i},{j}':str(B[i,j])for i,j in PAIRS if B[i,j]}for B in Rs],
                second_order_obstruction_groebner_basis=list(map(str,G)),
                obstruction_normal_form='(w0,w1,w2)*(w0,w1,w2,t), z=(w0+4t,w1+6t,w2-4t,t)',
                unique_projective_second_order_extendible_weight_block_direction=[4,6,-4,1],
                integrable_R={f'{i},{j}':str(R[i,j])for i,j in PAIRS},
                exact_curve='Q_(2/9)+eps*R-75*eps^2*x4^2*d0^d1',
                all_parameter_raw_Jacobi_zero=True,all_parameter_divergence_zero=True,
                all_parameter_generic_rank=4,degeneracy_polynomial=str(F),
                full_35_variable_formal_model_completed=False,new_component_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
                seconds=time.monotonic()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify()
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS resonance extra direction integrates to an exact quadratic curve',result['seconds'])

#!/usr/bin/env python3
"""Linear equations for the marked divisor and a finite local length bound for the fixed pair."""
import argparse
import hashlib
from itertools import combinations,combinations_with_replacement
import json
from pathlib import Path
import time

import sympy as sp

from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61, X

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_AFFINE_ANCHOR_LENGTH_BOUND_20261102'
PAIRS=list(combinations(range(4),2))
MONOMIALS=list(combinations_with_replacement(range(4),2))


def serialize_polynomials(matrix,variables):
    output=[]
    for row in matrix.tolist():
        result=[]
        for column,expr in enumerate(row):
            if not expr:
                continue
            terms=[]
            for exponents,coefficient in sp.Poly(expr,*variables,domain=sp.QQ).terms():
                monomial=[j for j,power in enumerate(exponents) for _ in range(power)]
                assert len(monomial)==2
                terms.append([*monomial,str(coefficient)])
            result.append([column,terms])
        output.append(result)
    return output


def raw_coordinates(model):
    q=sp.symbols('q0:45')
    iterator=iter(q)
    P=model.reconstruct([next(iterator) if w==2 else 0 for w in model.reduced['parameter_weights']])
    entries=[]
    for i,j in PAIRS:
        form=(P[1+i,1+j]-X[1+i]*P[0,1+j]+X[1+j]*P[0,1+i]).subs(X[0],1)
        polynomial=sp.Poly(form,*X[1:])
        entries.extend(polynomial.coeff_monomial(X[a+1]*X[b+1]) for a,b in MONOMIALS)
    E=sp.Matrix(entries).jacobian(q)
    assert E.shape==(60,45) and E.rank()==45
    return E


def linear_forms(beta,variables):
    # linear[i,k] is beta_{e_i}(p_k); star takes p to [e_i,e_j] = vol(p,e_i,e_j).
    star=sp.Matrix([[sp.LeviCivita(*(ab+ij)) for ab in PAIRS] for ij in PAIRS])
    assert star*star==sp.eye(6)
    linear=sp.Matrix([[beta[a,b].diff(y) for a,b in PAIRS] for y in variables])
    d=linear*star
    raw=sp.zeros(6,60)
    def index(i,j,h):
        return 10*h+MONOMIALS.index(tuple(sorted((i,j))))
    def polar_scale(i,j):
        return sp.Integer(1) if i==j else sp.Rational(1,2)
    for row,(a,b) in enumerate(PAIRS):
        for i,j,sign in [(a,b,1),(b,a,-1)]:
            for n,(k,l) in enumerate(PAIRS):
                coeff=25*sign*d[i,n]
                for h in range(6):
                    raw[row,index(j,l,h)]+=coeff*linear[k,h]*polar_scale(j,l)
                    raw[row,index(j,k,h)]-=coeff*linear[l,h]*polar_scale(j,k)
    return raw.applyfunc(sp.expand),linear,star


def check_standard_plane(raw,linear,star,u):
    sub={u[2]:0,u[3]:0}
    assert [str(x)for x in sub]==['u12','u13']
    alpha_x=5*u[4]-u[13]
    alpha_w=5*u[9]+u[14]
    p=sp.zeros(6,1)
    p[PAIRS.index((1,2))]=1
    alpha=sp.Matrix([alpha_x,0,0,alpha_w])
    brackets=sp.Matrix([sp.eye(4)[:,i].row_join(alpha).extract([0,3],[0,1]).det()
                        for i in range(4)])
    assert (brackets-5*linear.subs(sub)*p).applyfunc(sp.expand)==sp.zeros(4,1)
    f=sp.symbols('f0:3')
    q=sp.zeros(60,1)
    for idx,mon in enumerate([(0,0),(0,3),(3,3)]):
        q[10*PAIRS.index((1,2))+MONOMIALS.index(mon)]=f[idx]
    marked=f[0]*alpha_x**2+f[1]*alpha_x*alpha_w+f[2]*alpha_w**2
    assert (raw.subs(sub)*q-(star*p)*marked).applyfunc(sp.expand)==sp.zeros(6,1)
    return sp.sstr(marked)


def verify():
    start=time.monotonic()
    model=S61()
    spatial,u,beta=primitive_beta(model)
    E=raw_coordinates(model)
    raw,linear,star=linear_forms(beta,spatial)
    identity=check_standard_plane(raw,linear,star,u)
    qforms=(raw*E).applyfunc(sp.expand)
    deps=['S61_AFFINE_BOUNDARY_FINITE_ANCHORS_20261102','S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102',
          'S61_AFFINE_SURFACE_BOUND_20261102','S61_AFFINE_DEGREE_TEN_EXCLUSION_20261102',
          'S61_GENERIC_BOUNDARY_INDEX_20261102']
    files=[Path(__file__),ROOT/'research/notes'/(STEM+'.md')]
    for name in deps:
        f=ROOT/'reports'/(name+'.json')
        assert json.loads(f.read_text())['all_passed']
        files.extend([f,ROOT/'research/notes'/(name+'.md')])
    files.extend(ROOT/'code'/f for f in ['s61_model.py','primitive_jet_elliptic_curve_20261101.py',
                                        's61_infinity_flags.py','elliptic_weight1_dominance.py',
                                        'weight1_structure.py'])
    files.extend(ROOT/'legacy/v02/data'/f for f in ['symplectic_slice_certificate.json',
                'symplectic_slice_equations.json','symplectic_slice_reduced_equations.json'])
    return dict(id=STEM,all_passed=True,beta_variable_order=list(map(str,u)),
                weight2_u_indices_zero_based=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==2],
                bracket_order=PAIRS,spatial_monomial_order=MONOMIALS,volume_complement_matrix=star.tolist(),
                beta2_times_q_linear_coefficients=serialize_polynomials(qforms,u),
                auxiliary_chart_raw_tensor='g_{H,v}=sum_i H_i*(star*p)_i evaluated by pairing H with Q(v,v), i.e. H(Q(v,v))',
                marked_restriction_in_standard_plane=identity,
                all_six_standard_identities_over_14_variables=True,
                fixed_pair_generic_rank_equals='a+5e, where a is the marked intersection degree and e the degree over C of the boundary cycle',
                finite_anchor_inequality='0 <= a <= rank <= 2a',marked_intersection_degree_divisible_by=5,
                sufficient_finite_anchor_total_length_for_no_surface=9,
                no_new_numerical_upper_on_actual_rank_claimed=True,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files},
                seconds=time.monotonic()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=verify()
    with args.output.open('x')as stream:
        json.dump(report,stream,indent=2,default=int)
        stream.write('\n')
    print('PASS',round(report['seconds'],3))

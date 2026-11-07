#!/usr/bin/env python3
"""One finite, exact Rees-contraction ansatz for the entire quadratic_A source.

The SS expression is C+tau*A+tau²*B with H-weight-zero quadratic A,B.
All Laurent negative powers, endpoint and full Jacobi identities are imposed.
No jets are promoted to actual families, and tau degree is never enlarged.
"""
import argparse
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

import sympy as s
from classification_nil22_image_source_review import build_P,coords,p,F0,mag,x,y,q1,q2

ROOT=Path(__file__).resolve().parents[1]


def jacobi(P):
    return [s.expand(sum(P[i,l]*s.diff(P[j,k],coords[l])+P[j,l]*s.diff(P[k,i],coords[l])+P[k,l]*s.diff(P[i,j],coords[l])for l in range(4)))for i,j,k in it.combinations(range(4),3)]


def div(P):return s.Matrix([s.expand(sum(s.diff(P[i,j],coords[i])for i in range(4)))for j in range(4)])


def push(P,G):
    return (G*P*G.T).subs(dict(zip(coords,G.inv()*s.Matrix(coords))),simultaneous=True).applyfunc(s.expand)


def equations(expressions,variables):
    return [c for expression in expressions for e,c in s.Poly(s.expand(expression),*variables).terms()if c!=0]


def encode(P):return {f'{i},{j}':str(s.factor(P[i,j]))for i,j in it.combinations(range(4),2)if P[i,j]!=0}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    lam,tau=s.symbols('lam tau');a=lam-s.Rational(1,2);b=-lam-s.Rational(1,2)
    P,_,_=build_P(0,s.Matrix([[0,1],[1,0]]),s.zeros(2,1),s.diag(a,b))
    # Full magnetic affine normalizer, with its exact nonzero complement.
    I=mag[0]+mag[4]
    L1=mag[2]/b-mag[4]*x+mag[5]*y/(2*b);L2=-mag[1]/a-mag[3]*x/(2*a)
    changed_coords=s.Matrix([q1+L1,q2+L2,x,y]);G=changed_coords.jacobian(s.Matrix(coords))
    translated=(G*P*G.T).subs({q1:q1-L1,q2:q2-L2},simultaneous=True).applyfunc(s.cancel)
    canonical=P.subs(dict(zip(mag,[I,0,0,0,0,0])),simultaneous=True)
    assert (translated-canonical).applyfunc(s.cancel)==s.zeros(4)
    canonical_one=P.subs(dict(zip(mag,[1,0,0,0,0,0])))
    J=s.Symbol('I',nonzero=True);scaling=s.diag(1/J,1/J,1,1)
    normalized=push(P.subs(dict(zip(mag,[J,0,0,0,0,0]))),scaling)*J
    assert (normalized-canonical_one).applyfunc(s.cancel)==s.zeros(4)
    shift=s.eye(4);shift[0,2]=1
    target=push(canonical_one,shift)
    C=s.zeros(4);C[0,3]=C[1,2]=1;C[3,0]=C[2,1]=-1
    assert target.subs(dict.fromkeys(coords,0))==C
    assert not any(jacobi(target))
    N=s.zeros(4);N[0,2]=N[1,3]=1
    H=s.diag(-1,-1,1,1);g=s.eye(4)+N/(2*tau)
    assert g.det()==1 and g*C*g.T==C
    assert (g*(tau*H)*g.inv()-N-tau*H).applyfunc(s.cancel)==s.zeros(4)
    weights=[-1,-1,1,1];basis=[]
    for i,j in it.combinations(range(4),2):
        for k,l in it.combinations_with_replacement(range(4),2):
            if weights[k]+weights[l]!=weights[i]+weights[j]:continue
            B=s.zeros(4);B[i,j]=coords[k]*coords[l];B[j,i]=-B[i,j];basis.append(B)
    assert len(basis)==22
    unknown=s.symbols('r0:44');A=sum([co*B for co,B in zip(unknown[:22],basis)],s.zeros(4));B=sum([co*Q for co,Q in zip(unknown[22:],basis)],s.zeros(4))
    ss=C+tau*A+tau*tau*B
    family=push(ss,g)
    # Weightzeroquadratics have at most two powers of the unipotent shear.
    # Multiplication by tau clears every possible negative power here.
    assert all(s.denom(s.cancel(tau*value)).free_symbols<={lam} for value in family)
    linear_equations=equations(div(A)-H*s.Matrix(coords),coords)+equations(div(B),coords)
    linear_equations+=equations([v-z for v,z in zip(jacobi(C+A),jacobi(A))],coords)
    linear_equations+=equations([v-z for v,z in zip(jacobi(C+B),jacobi(B))],coords)
    endpoint=[]
    for value,target_value in zip(family,target):
        endpoint +=[s.expand(value).coeff(tau,-1),s.expand(value).coeff(tau,0)-target_value]
    linear_equations+=equations(endpoint,coords)
    matrix,rhs=s.linear_eq_to_matrix(linear_equations,unknown)
    solution_set=s.linsolve((matrix,rhs),unknown)
    assert solution_set is not s.EmptySet and bool(solution_set),'This fixed degree2 contraction has no linear lift.'
    solution=list(solution_set)[0]
    free=sorted(set().union(*(value.free_symbols for value in solution))-{lam},key=str)
    replacement=dict(zip(unknown,solution));A=A.subs(replacement,simultaneous=True).applyfunc(s.expand);B=B.subs(replacement,simultaneous=True).applyfunc(s.expand)
    ss=C+tau*A+tau*tau*B
    nonlinear=equations(jacobi(ss),(tau,*coords))
    field=s.QQ.frac_field(lam)
    if free:
        nonlinear=list({str(s.Poly(value,*free,domain=field).monic().as_expr()):s.Poly(value,*free,domain=field).monic().as_expr()for value in nonlinear}.values())
    progress=dict(status='full_linear_reduction_complete',tau_degree=2,linear_shape=list(matrix.shape),linear_rank=matrix.rank(),remaining_variables=list(map(str,free)),
        exact_A=encode(A),exact_B=encode(B),full_Jacobi_conditions=list(map(str,nonlinear)))
    Path('/tmp/nil22_quadratic_ss_limit_progress.json').write_text(json.dumps(progress,indent=2)+'\n')
    certificate=[];candidates=[]
    if not nonlinear:candidates=[{}]
    elif not free:certificate=['1']
    else:
        GB=s.groebner(nonlinear,*free,domain=field,order='grevlex');certificate=list(map(str,GB.polys))
        if not(GB.is_zero_dimensional is False and len(GB.polys)==1 and GB.polys[0].as_expr()==1):
            if not(len(GB.polys)==1 and GB.polys[0].as_expr()==1):
                candidates=s.solve([v.as_expr()for v in GB.polys],free,dict=True,check=False,simplify=False)
    accepted=None
    for candidate in candidates:
        values={key:s.cancel(value)for key,value in candidate.items()}
        if any(not value.is_rational_function(lam,*free)for value in values.values()):continue
        if any(s.cancel(value.subs(values,simultaneous=True))!=0 for value in nonlinear):continue
        accepted=values;break
    found=accepted is not None;result_family=None
    if found:
        A=A.subs(accepted,simultaneous=True).applyfunc(s.cancel);B=B.subs(accepted,simultaneous=True).applyfunc(s.cancel)
        ss=C+tau*A+tau*tau*B;family=push(ss,g).applyfunc(s.cancel)
        assert not any(jacobi(ss)) and not any(jacobi(family))
        assert (div(ss)-tau*H*s.Matrix(coords)).applyfunc(s.cancel)==s.zeros(4,1)
        assert (div(family)-(N+tau*H)*s.Matrix(coords)).applyfunc(s.cancel)==s.zeros(4,1)
        assert (family.subs(tau,0)-target).applyfunc(s.cancel)==s.zeros(4)
        assert all(s.denom(value).has(tau)is False for value in family)
        assert all(s.Poly(value,tau).degree()<=2 for value in family if value!=0)
        assert family.subs(dict.fromkeys(coords,0))==C
        result_family=dict(source_SS_A=encode(A),source_SS_B=encode(B),whole_regular_contracted_family=encode(family),
            exact_endpoint=encode(target),all_tau_strictJacobi_zero=True,SS_modular='tau*diag(-1,-1,1,1)',
            contracted_modular='N+tau*diag(-1,-1,1,1)',constant_symplectic_Pfaffian='1',
            every_nonzero_tau_in_accepted_SS20_rank2_source=True,
            accepted_coefficient_substitution={str(k):str(v)for k,v in accepted.items()})
    inputs=[Path(__file__),ROOT/'code/classification_nil22_image_source_review.py',ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_SOURCE_REVIEW_V2_20261107.json',ROOT/'research/classification/HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_SUPPORT.md']
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,exact_fixed_ansatz_completed=True,tau_degree_not_enlarged=2,
        full_seven_original_source_parameters=[*map(str,mag),'lam'],magnetic_invariant='I=m0+m4',
        exact_all_magnetic_affine_normalizer=[str(L1),str(L2)],normalizer_open='(lam-1/2)*(lam+1/2)*I !=0',
        complement='I0 explicitF0=tau Poissoncurve; lambda±1/2 and every other genericparametercomplement lie in closure of the originalirreducible A7source.',
        linear_reduction=progress,Groebner_basis_over_QQ_lambda=certificate,actual_curve_found=found,
        entire_quadratic_A_source_MPS_covered=found,actual_family=result_family,
        scope='Positive means literal all-tau Poissonfamily and wholeoriginalsource closedness; negative excludes only this specified degree2 SS20contraction, not otherfamilies or components.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_source_bytes=Path(__file__).read_text(),
        seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:result[k]for k in ['all_passed','actual_curve_found','entire_quadratic_A_source_MPS_covered','seconds']}))


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Finite whole-source support reduction for the simple-H nil31 short translation."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import time
from itertools import combinations
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json'

def ranks_with_exceptional_parameters(matrix, parameter):
    generic_rank=matrix.rank()
    minors=[]
    for rows in combinations(range(matrix.rows),generic_rank):
        for cols in combinations(range(matrix.cols),generic_rank):
            minor=s.factor(matrix.extract(rows,cols).det())
            if minor:
                minors.append(s.fraction(minor)[0])
    g=s.factor(s.polys.polytools.terms_gcd(s.polys.polytools.gcd_list(minors,parameter)))
    exceptional=[]
    for value in s.solve(g,parameter):
        if value.is_rational:
            exceptional.append(dict(value=str(value),rank=matrix.subs(parameter,value).rank()))
        else:
            raise ValueError('Handle algebraic exceptional parameters explicitly: '+str(value))
    return dict(generic_rank=generic_rank,maximum_minor_gcd=str(g),exceptional=exceptional,
                matrix=[[str(q)for q in row]for row in matrix.tolist()])

def whole_core_branches(n,core):
    c=n[17];b=n[16];t=n[15];k=n[9];a=n[11];j=n[10]
    branches=[]
    for value in [s.Rational(1,3),s.Rational(1,6),-s.Rational(1,12),s.Rational(1,15)]:
        for solution in s.solve([q.subs(b,value) for q in core],[j,a,k,t],dict=True):
            assert set(solution)=={j,a,k,t},solution
            sub={**solution,b:value}
            branches.append(dict(name='n17_nonzero_b_'+str(value)+'_t_'+str(sub[t]),
                parameter=c,nonzero=True,sub=sub))
    branches += [
        dict(name='n17_zero_n15_nonzero',parameter=t,nonzero=True,sub={c:0,j:-2*t,a:-s.Rational(1,2),b:s.Rational(1,3),k:-3*t*t/2}),
        dict(name='n17_zero_n15_zero_n9_nonzero',parameter=k,nonzero=True,sub={c:0,j:0,a:s.Rational(1,3),b:-s.Rational(1,3),t:0}),
        dict(name='n17_zero_n15_n9_zero_b_third',parameter=a,nonzero=False,sub={c:0,j:0,b:s.Rational(1,3),t:0,k:0}),
        dict(name='n17_zero_n15_n9_zero_a_minus_b',parameter=b,nonzero=False,sub={c:0,j:0,a:-b,t:0,k:0})]
    for branch in branches:
        assert all(s.cancel(q.subs(branch['sub'],simultaneous=True))==0 for q in core),branch
    return branches

def source_dimension_bounds(n,equations,core):
    lower=[n[i]for i in [3,4,8,12,13,14]]
    magnetic=[n[i]for i in [0,1,2,5,6,7]]
    out=[]
    for branch in whole_core_branches(n,core):
        parameter=branch['parameter'];sub=branch['sub']
        eq=[s.cancel(q.subs(sub,simultaneous=True))for q in equations]
        eq=[q for q in eq if q]
        leq=[q for q in eq if not(q.free_symbols&set(magnetic))]
        meq=[q for q in eq if q.free_symbols&set(magnetic)]
        lm=s.linear_eq_to_matrix(leq,lower)[0]
        mm=s.linear_eq_to_matrix(meq,magnetic)[0]
        assert not (set().union(*(q.free_symbols for q in mm))&set(lower))
        lrank=ranks_with_exceptional_parameters(lm,parameter)
        mrank=ranks_with_exceptional_parameters(mm,parameter)
        dimension=13-lrank['generic_rank']-mrank['generic_rank']
        values=set(q['value']for q in lrank['exceptional']+mrank['exceptional'])
        special=[]
        for value in sorted(values):
            v=s.Rational(value)
            if branch['nonzero']and not v:continue
            lr=lm.subs(parameter,v).rank();mr=mm.subs(parameter,v).rank()
            special.append(dict(parameter=value,lower_rank=lr,magnetic_rank=mr,fixed_source_dimension_upper_bound=12-lr-mr))
        upper=max([dimension]+[x['fixed_source_dimension_upper_bound']for x in special])
        out.append(dict(name=branch['name'],parameter=str(parameter),parameter_nonzero=branch['nonzero'],
            substitution={str(k):str(v)for k,v in sub.items()},all_core_relations_zero=True,
            lower_variables=list(map(str,lower)),magnetic_variables=list(map(str,magnetic)),
            lower_rank=lrank,magnetic_rank=mrank,generic_fixed_source_dimension_upper_bound=dimension,
            specializations=special,fixed_source_dimension_upper_bound=upper))
    return out

def main(output):
    resource.setrlimit(resource.RLIMIT_CPU, (110,115))
    resource.setrlimit(resource.RLIMIT_AS, (1536*1024**2,1536*1024**2))
    begun=time.monotonic()
    original=json.loads(SOURCE.read_text())['full_nilpotent_sources'][3]
    n=s.symbols('n0:18')
    equations=[s.sympify(q,locals={str(v):v for v in n}) for q in original['distinct_nonzero_jacobi_coefficients']]
    equations=list(dict.fromkeys(s.Poly(q,*n).monic().as_expr() for q in equations))
    corevars=[n[10],n[11],n[9],n[15],n[16],n[17]]
    core=[q for q in equations if q.free_symbols<=set(corevars)]
    # The localization/complement is explicit; no full-support conclusion is
    # inferred solely from the fraction-field branch.
    nonzero=s.groebner(core,*corevars[:-1],domain=s.QQ.frac_field(n[17]))
    zero=s.groebner([q.subs(n[17],0) for q in core],*corevars[:-1],domain=s.QQ)
    bounds=source_dimension_bounds(n,equations,core)
    result=dict(task_id='G',scope='The whole fixed-cofactor quadratic Poisson source has dimension at most7; its projective PGL saturation has dimension at most23. No component membership or global exhaustion is claimed.',
        field=original['field_linear_matrix'],translation=original['field_constant'],
        full_source_parameter_count=18,all_reduced_jacobi_equations=list(map(str,equations)),
        core_variables=list(map(str,corevars)),core_equations=list(map(str,core)),
        n17_nonzero_core_groebner=list(map(str,nonzero.polys)),
        n17_zero_core_groebner=list(map(str,zero.polys)),
        n17_nonzero_core_expressions=[str(p.as_expr()) for p in nonzero.polys],
        n17_zero_core_expressions=[str(p.as_expr()) for p in zero.polys],
        n17_nonzero_last_core_polynomial_factorization=str(s.factor(nonzero.polys[-1].as_expr())),
        original_pfaffian=original['actual_affine_pfaffian'],
        whole_source_dimension_bounds=bounds,
        fixed_source_dimension_upper_bound=max(q['fixed_source_dimension_upper_bound']for q in bounds),
        affine_cofactor_orbit_dimension=13,hyperplane_parameter_dimension=4,projective_scale_dimension=1,
        saturated_projective_dimension_upper_bound=23,
        arithmetic_passed=True,independent_review='pending',
        input_sha256={str(SOURCE.relative_to(ROOT)):hashlib.sha256(SOURCE.read_bytes()).hexdigest()},
        source_bytes=Path(__file__).read_text(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seconds=time.monotonic()-begun)
    with output.open('x') as f:
        json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k] for k in ['scope','fixed_source_dimension_upper_bound','seconds']}))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    main(p.parse_args().output)

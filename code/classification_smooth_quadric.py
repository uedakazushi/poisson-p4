#!/usr/bin/env python3
"""Full invariant quadratic lifts for a smooth quadric normal Hamiltonian.

This finite entrance calculation does not assert component membership.
All ten homogeneous brackets and all 350 Jacobi coefficient slots are kept.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s


def run_nil(kind, fixed_spec=None):
    """Complete fixed-cofactor sources; optionally use a rank3/4 table row."""
    started=time.monotonic()
    x=s.symbols('x0:5')
    N=s.zeros(5)
    if fixed_spec is not None:
        local=dict(zip(map(str,x),x))
        q=s.sympify(fixed_spec['quadric'],locals=local)
        components=s.Matrix([s.sympify(f,locals=local) for f in fixed_spec['cofactor_components']])
        N=components.jacobian(x)
        assert components==N*s.Matrix(x)
        assert s.hessian(q,x).rank()==fixed_spec['quadric_rank']
    elif kind == 'nil5':
        q=x[0]*x[4]-4*x[1]*x[3]+3*x[2]**2
        for i in range(4):N[i+1,i]=i+1
    elif kind == 'nil311':
        q=x[0]*x[2]-x[1]**2+x[3]**2+x[4]**2
        N[1,0]=1;N[2,1]=2
    elif kind == 'nil221':
        q=x[0]*x[3]-x[1]*x[2]+x[4]**2
        N[1,0]=1;N[3,2]=1
    else:
        raise ValueError(kind)
    E=s.Matrix(x);V=N*E
    assert s.expand(sum(V[i]*s.diff(q,x[i]) for i in range(5)))==0
    assert N.trace()==0 and N!=s.zeros(5)
    if fixed_spec is None:
        assert s.hessian(q,x).det()!=0
    pairs=list(itertools.combinations(range(5),2))
    mons=[x[i]*x[j] for i in range(5) for j in range(i,5)]
    unknowns=s.symbols('p0:150')
    P=s.zeros(5)
    for at,(i,j) in enumerate(pairs):
        P[i,j]=sum(unknowns[15*at+k]*mon for k,mon in enumerate(mons))
        P[j,i]=-P[i,j]
    constraints=[sum(s.diff(P[i,j],x[i]) for i in range(5)) for j in range(5)]
    constraints += [sum(s.diff(q,x[i])*P[i,j] for i in range(5))-q*V[j] for j in range(5)]
    lie=P.applyfunc(lambda f:sum(V[i]*s.diff(f,x[i]) for i in range(5)))-N*P-P*N.T
    constraints += [lie[i,j] for i,j in pairs]
    equations=[c for f in constraints for c in s.Poly(s.expand(f),*x).coeffs() if c!=0]
    M,b=s.linear_eq_to_matrix(equations,unknowns)
    solutions=s.linsolve((M,b),unknowns)
    if solutions is s.EmptySet:
        source=Path(__file__).read_text()
        return {'scope':'empty full150 linear fixed-cofactor source','case':kind,
                'fixed_specification':fixed_spec,'linear_matrix_shape':list(M.shape),
                'linear_rank':M.rank(),'augmented_rank':M.row_join(b).rank(),
                'elapsed_seconds':time.monotonic()-started,
                'producer_sha256':hashlib.sha256(source.encode()).hexdigest(),'producer_source':source}
    sol=next(iter(solutions))
    oldparams=sorted(set().union(*(f.free_symbols for f in sol)),key=str)
    parameters=s.symbols(f'n0:{len(oldparams)}')
    rename=dict(zip(oldparams,parameters))
    substitution={p:f.subs(rename,simultaneous=True) for p,f in zip(unknowns,sol)}
    P=P.subs(substitution,simultaneous=True).applyfunc(s.expand)
    assert all(s.expand(f.subs(substitution,simultaneous=True))==0 for f in constraints)
    jacobi=[]
    for i,j,k in itertools.combinations(range(5),3):
        f=s.expand(sum(P[i,r]*s.diff(P[j,k],x[r])+P[j,r]*s.diff(P[k,i],x[r])
                        +P[k,r]*s.diff(P[i,j],x[r]) for r in range(5)))
        jacobi += [s.factor(c) for c in s.Poly(f,*x).coeffs() if c!=0]
    jacobi=sorted(set(jacobi),key=str)
    pf=[]
    for missing in range(5):
        i,j,k,l=[r for r in range(5) if r != missing]
        pf.append(s.expand((-1)**missing*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])))
    potential=s.factor(sum(x[i]*pf[i] for i in range(5))/5)
    source=Path(__file__).read_text()
    return {'scope':'complete linear fixed-cofactor source; all residual Jacobi retained',
        'case':kind,'fixed_specification':fixed_spec,'q':str(q),'normal_hamiltonian_matrix':str(N),
        'linear_matrix_shape':list(M.shape),'linear_rank':150-len(parameters),
        'source_dimension':len(parameters),'parameters':list(map(str,parameters)),
        'brackets':{f'{i}{j}':str(P[i,j]) for i,j in pairs},
        'residual_jacobi_coefficients':list(map(str,jacobi)),
        'euler_contracted_pfaffian_not_yet_actual_divisor':str(potential),
        'elapsed_seconds':time.monotonic()-started,
        'producer_sha256':hashlib.sha256(source.encode()).hexdigest(),'producer_source':source}


def run_rank2(mixed=False):
    """All lifts for the semisimple normal Hamiltonian with a zero3 block."""
    started = time.monotonic()
    x = s.symbols("x0:5")
    weights = (0,0,0,1,-1)
    q = sum(x[i]**2 for i in range(3))+x[3]*x[4]
    pairs = list(itertools.combinations(range(5),2))
    monomials = [(x[i]*x[j],weights[i]+weights[j]) for i in range(5) for j in range(i,5)]
    unknowns=[]
    universal=s.zeros(5)
    for i,j in pairs:
        for mon,weight in monomials:
            if weight == weights[i]+weights[j]:
                p=s.Symbol(f"p{len(unknowns)}")
                unknowns.append(p)
                universal[i,j]+=p*mon
        universal[j,i]=-universal[i,j]
    linear_polys=[sum(s.diff(universal[i,j],x[i]) for i in range(5)) for j in range(5)]
    linear_polys += [sum(s.diff(q,x[i])*universal[i,j] for i in range(5))
                     -q*weights[j]*x[j] for j in range(5)]
    equations=[c for f in linear_polys for c in s.Poly(s.expand(f),*x).coeffs()]
    M,rhs=s.linear_eq_to_matrix(equations,unknowns)
    linear_rank=M.rank()
    assert M.row_join(rhs).rank()==linear_rank
    ell=s.symbols('l0:3')
    mm=s.symbols('m0:8')
    A=s.Matrix([[mm[0],mm[1],mm[2]],[mm[3],mm[4],mm[5]],[mm[6],mm[7],-1-mm[0]-mm[4]]])
    K=s.Matrix([[0,ell[2],-ell[1]],[-ell[2],0,ell[0]],[ell[1],-ell[0],0]])
    U=s.Matrix(x[:3])
    linear_ell=sum(ell[i]*x[i] for i in range(3))
    J=s.Matrix([[0,x[2],-x[1]],[-x[2],0,x[0]],[x[1],-x[0],0]])
    P=s.zeros(5)
    P[:3,:3]=linear_ell*J-x[3]*x[4]*K/2
    for i in range(3):
        P[i,3]=x[3]*(A*U)[i]
        P[3,i]=-P[i,3]
        P[i,4]=x[4]*((K-A)*U)[i]
        P[4,i]=-P[i,4]
    P[3,4]=(2*(U.T*A*U)[0]-sum(x[i]**2 for i in range(3)))-x[3]*x[4]
    P[4,3]=-P[3,4]
    assert all(s.expand(sum(s.diff(P[i,j],x[i]) for i in range(5)))==0 for j in range(5))
    assert all(s.expand(sum(s.diff(q,x[i])*P[i,j] for i in range(5))-q*weights[j]*x[j])==0 for j in range(5))
    parameters=ell+mm
    coefficient_columns=[s.Poly(P[i,j],*x).as_dict() for i,j in pairs]
    source_vector=s.Matrix([d.get(mon,0) for d in coefficient_columns
        for mon in sorted(set().union(*(e.keys() for e in coefficient_columns)))])
    assert source_vector.jacobian(parameters).rank()==11
    assert len(unknowns)-linear_rank==11
    mixed_linear_data=None
    if mixed:
        N=s.zeros(5)
        N[:3,:3]=s.Matrix([[0,0,-s.I],[0,0,1],[s.I,-1,0]])
        E=s.Matrix(x)
        V=N*E
        P=P-5*q*N/6-(E*V.T-V*E.T)/3
        assert all(s.expand(sum(s.diff(P[i,j],x[i]) for i in range(5)))==0 for j in range(5))
        assert all(s.expand(sum(s.diff(q,x[i])*P[i,j] for i in range(5))
                            -q*(weights[j]*x[j]+V[j]))==0 for j in range(5))
        lie=P.applyfunc(lambda f:sum(V[i]*s.diff(f,x[i]) for i in range(5)))-N*P-P*N.T
        lie_eq=[c for i,j in pairs for c in s.Poly(s.expand(lie[i,j]),*x).coeffs() if c!=0]
        LM,lb=s.linear_eq_to_matrix(lie_eq,parameters)
        sol=next(iter(s.linsolve((LM,lb),parameters)))
        sub=dict(zip(parameters,sol))
        P=P.subs(sub,simultaneous=True).applyfunc(s.expand)
        mixed_linear_data={'nilpotent_matrix':str(N),'invariance_rank':LM.rank(),
            'substitution':{str(k):str(v) for k,v in sub.items()},
            'remaining_parameters':sorted(map(str,set().union(*(f.free_symbols for f in P))-set(x)))}
    jacobi=[]
    for i,j,k in itertools.combinations(range(5),3):
        poly=s.expand(sum(P[i,r]*s.diff(P[j,k],x[r])+P[j,r]*s.diff(P[k,i],x[r])
                          +P[k,r]*s.diff(P[i,j],x[r]) for r in range(5)))
        jacobi += [s.factor(c) for c in s.Poly(poly,*x).coeffs() if c!=0]
    jacobi=sorted(set(jacobi),key=str)
    if mixed:
        reduced_sub={ell[1]:-5*s.I/3,mm[4]:-s.Rational(1,3)}
        assert all(s.expand(f.subs(reduced_sub))==0 for f in jacobi)
        final_P=P.subs(reduced_sub).applyfunc(s.expand)
        pfs=[]
        for missing in range(5):
            i,j,k,l=[r for r in range(5) if r != missing]
            pfs.append(s.expand((-1)**missing*(final_P[i,j]*final_P[k,l]-final_P[i,k]*final_P[j,l]+final_P[i,l]*final_P[j,k])))
        actual=s.factor(sum(x[i]*pfs[i] for i in range(5))/5,extension=s.I)
        assert all(s.expand(s.diff(actual,x[i])-pfs[i])==0 for i in range(5))
        assert s.expand(actual-s.Rational(5,6)*(x[0]+s.I*x[1])*sum(t*t for t in x[:3])*q)==0
        source=Path(__file__).read_text()
        return {'scope':'complete mixed normal-Hamiltonian source and all residual Jacobi',
            'case':'semisimple zero3 block plus nilpotent Jordan3 on that block',
            'q':str(q),'weights':list(weights),'linear_rank_before_invariance':linear_rank,
            'mixed_invariance':mixed_linear_data,
            'brackets':{f'{i}{j}':str(P[i,j]) for i,j in pairs},
            'residual_jacobi_coefficients':list(map(str,jacobi)),
            'full_poisson_source_substitution':{str(k):str(v) for k,v in reduced_sub.items()},
            'full_poisson_source_free_parameter':str(mm[7]),
            'actual_divisor':str(actual),'irreducible_cubic_factor_impossible':True,
            'elapsed_seconds':time.monotonic()-started,
            'producer_sha256':hashlib.sha256(source.encode()).hexdigest(),'producer_source':source}
    support=[]
    for name,ell_value in [('zero',(0,0,0)),('nonisotropic',(1,0,0)),('isotropic',(1,s.I,0))]:
        # Nonzero isotropic vectors form one SO3 orbit.  In the nonisotropic
        # case a nonzero norm scalar factors out of every Jacobi equation.
        ell_sub=dict(zip(ell,ell_value))
        specialized=[s.expand(f.subs(ell_sub)) for f in jacobi]
        CM,cb=s.linear_eq_to_matrix(specialized,mm)
        solution=s.linsolve((CM,cb),mm)
        if solution is s.EmptySet:
            support.append({'case':name,'empty':True,'rank':CM.rank(),'augmented_rank':CM.row_join(cb).rank()})
            continue
        values=next(iter(solution))
        msub=dict(zip(mm,values))
        pp=P.subs(ell_sub).subs(msub,simultaneous=True).applyfunc(s.expand)
        assert all(s.expand(f.subs(msub,simultaneous=True))==0 for f in specialized)
        pfs=[]
        for missing in range(5):
            i,j,k,l=[r for r in range(5) if r != missing]
            pfs.append(s.expand((-1)**missing*(pp[i,j]*pp[k,l]-pp[i,k]*pp[j,l]+pp[i,l]*pp[j,k])))
        actual=s.factor(sum(x[i]*pfs[i] for i in range(5))/5,extension=s.I)
        assert all(s.expand(s.diff(actual,x[i])-pfs[i])==0 for i in range(5))
        cubic=s.cancel(actual/q)
        assert s.denom(cubic).free_symbols.isdisjoint(x)
        support.append({'case':name,'empty':False,'substitution':{str(k):str(v) for k,v in msub.items()},
                        'matrix_A':str(A.subs(msub,simultaneous=True)),
                        'actual_divisor':str(actual),'cubic_factor':str(s.factor(cubic,extension=s.I)),
                        'brackets':{f'{i}{j}':str(pp[i,j]) for i,j in pairs}})
    source=Path(__file__).read_text()
    return {'scope':'complete linear source; residual Jacobi retained, no full support conclusion',
        'case':'smooth quadric, semisimple normal Hamiltonian (0,0,0,1,-1)',
        'q':str(q),'weights':list(weights),'invariant_unknowns':len(unknowns),
        'linear_shape':list(M.shape),'linear_rank':linear_rank,'source_dimension':11,
        'parameters':list(map(str,parameters)),'matrix_A':str(A),'matrix_K':str(K),
        'brackets':{f'{i}{j}':str(s.expand(P[i,j])) for i,j in pairs},
        'residual_jacobi_coefficients':list(map(str,jacobi)),
        'complete_ell_orbit_cases':support,
        'elapsed_seconds':time.monotonic()-started,
        'producer_sha256':hashlib.sha256(source.encode()).hexdigest(),'producer_source':source}


def run(ambient=False):
    started = time.monotonic()
    x = s.symbols("x0:5")
    weights = (0, 1, -1, 2, -2)
    q = x[0] ** 2 + x[1] * x[2] + x[3] * x[4]
    pairs = list(itertools.combinations(range(5), 2))
    monomials = [(x[i] * x[j], weights[i] + weights[j])
                 for i in range(5) for j in range(i, 5)]
    unknowns = []
    P = s.zeros(5)
    for i, j in pairs:
        for mon, weight in monomials:
            if weight == weights[i] + weights[j]:
                c = s.Symbol(f"p{len(unknowns)}")
                unknowns.append(c)
                P[i, j] += c * mon
        P[j, i] = -P[i, j]

    linear_polynomials = [sum(s.diff(P[i, j], x[i]) for i in range(5))
                          for j in range(5)]
    # Convention: X(f) = {q,f}/q, hence dq_i P^{ij}=q X^j.
    linear_polynomials += [sum(s.diff(q, x[i]) * P[i, j] for i in range(5))
                          - q * weights[j] * x[j] for j in range(5)]
    equations = []
    for poly in linear_polynomials:
        equations += s.Poly(s.expand(poly), *x).coeffs()
    M, b = s.linear_eq_to_matrix(equations, unknowns)
    solution = next(iter(s.linsolve((M, b), unknowns)))
    substitution = dict(zip(unknowns, solution))
    P = P.applyfunc(lambda f: s.expand(f.subs(substitution, simultaneous=True)))
    parameters = sorted(set().union(*(f.free_symbols for f in P)) - set(x), key=str)
    jacobi = []
    for i, j, k in itertools.combinations(range(5), 3):
        poly = s.expand(sum(P[i, r] * s.diff(P[j, k], x[r])
                            + P[j, r] * s.diff(P[k, i], x[r])
                            + P[k, r] * s.diff(P[i, j], x[r])
                            for r in range(5)))
        for coefficient in s.Poly(poly, *x).coeffs():
            if coefficient != 0:
                jacobi.append(s.factor(coefficient))
    jacobi = sorted(set(jacobi), key=str)
    pfaffians = []
    for missing in range(5):
        a, b, c, d = [i for i in range(5) if i != missing]
        pfaffians.append(s.factor((-1) ** missing * (
            P[a, b] * P[c, d] - P[a, c] * P[b, d] + P[a, d] * P[b, c])))
    a, b, u, v, t = (s.Symbol(name) for name in ("p12", "p14", "p15", "p16", "p18"))
    conditions = [a*b, a*(u+1), b*(v+t+3), t*(u+v-1),
                  8*(u+v)*t + 30*(u+v) + 7*t + 20]
    coefficient_monomials = sorted(set().union(*(
        set(s.Poly(f, *parameters).monoms()) for f in jacobi + conditions)))
    def coefficient_matrix(polys):
        rows = [s.Poly(f, *parameters).as_dict() for f in polys]
        return s.Matrix([[row.get(mon, 0) for mon in coefficient_monomials] for row in rows])
    assert coefficient_matrix(jacobi).rank() == 5
    assert coefficient_matrix(conditions).rank() == 5
    assert coefficient_matrix(jacobi + conditions).rank() == 5
    z = s.Symbol("z")
    branches = [
        ("I_central", {a:0,b:0,u:z,v:-s.Rational(2,3)-z,t:0}),
        ("I_a", {a:z,b:0,u:-1,v:s.Rational(1,3),t:0}),
        ("I_b", {a:0,b:z,u:s.Rational(7,3),v:-3,t:0}),
        ("II_central", {a:0,b:0,u:z,v:1-z,t:-s.Rational(10,3)}),
        ("II_a", {a:z,b:0,u:-1,v:2,t:-s.Rational(10,3)}),
        ("II_b", {a:0,b:z,u:s.Rational(2,3),v:s.Rational(1,3),t:-s.Rational(10,3)}),
    ]
    branch_reports = []
    for name, sub in branches:
        assert all(s.expand(f.subs(sub, simultaneous=True)) == 0 for f in jacobi)
        pf = [s.expand(f.subs(sub, simultaneous=True)) for f in pfaffians]
        D = s.factor(sum(x[i]*pf[i] for i in range(5))/5)
        assert all(s.expand(s.diff(D,x[i])-pf[i]) == 0 for i in range(5))
        cubic = s.cancel(D/q)
        assert s.denom(cubic).free_symbols.isdisjoint(x)
        degeneration_weights = ((0,0,3,1,1) if name == "I_b" else (0,3,0,1,1)
                                if name.startswith("I_") else (0,0,2,2,1)
                                if name == "II_b" else (0,2,0,1,2))
        limit = s.zeros(5)
        exponents = []
        branch_P = P.subs(sub, simultaneous=True)
        for i,j in pairs:
            for mon,coefficient in s.Poly(branch_P[i,j],*x).terms():
                power = degeneration_weights[i]+degeneration_weights[j]-sum(
                    mon[k]*degeneration_weights[k] for k in range(5))
                assert power >= 0
                exponents.append(power)
                if power == 0:
                    limit[i,j] += coefficient*s.prod(x[k]**mon[k] for k in range(5))
            limit[j,i] = -limit[i,j]
            assert s.cancel(limit[i,j]/(x[i]*x[j])).free_symbols.isdisjoint(x)
        limit_pf = []
        for missing in range(5):
            i,j,k,ell = [r for r in range(5) if r != missing]
            limit_pf.append((-1)**missing*(limit[i,j]*limit[k,ell]-limit[i,k]*limit[j,ell]+limit[i,ell]*limit[j,k]))
        limit_D = s.factor(sum(x[i]*limit_pf[i] for i in range(5))/5)
        assert s.expand(limit_D-s.Rational(10,3)*s.prod(x)) == 0
        branch_reports.append({"name":name, "substitution":{str(k):str(v) for k,v in sub.items()},
                               "actual_divisor":str(D), "cubic_factor":str(s.factor(cubic)),
                               "full_jacobi_zero":True, "all_pfaffian_minors_equal_gradient":True,
                               "diagonal_degeneration_weights":list(degeneration_weights),
                               "nonnegative_parameter_powers":sorted(set(exponents)),
                               "limit_is_log_canonical":True,"limit_actual_divisor":str(limit_D)})
    ambient_reports = []
    if ambient:
        from sympy.polys.matrices import DomainMatrix
        from classification_rank3_mixed_component_review import (
            COLUMNS, jacobi as full_jacobi, sparse, vector, minor_certificate)
        for name, sub in branches:
            if name.endswith("central"):
                continue
            point_sub = {k:s.sympify(v).subs(z,1) for k,v in sub.items()}
            point = P.subs(point_sub, simultaneous=True)
            scale = s.ilcm(*[s.denom(c) for f in point for c in s.Poly(f,*x).coeffs()])
            point = point.applyfunc(lambda f:s.expand(scale*f))
            assert all(s.denom(c) == 1 for f in point for c in s.Poly(f,*x).coeffs())
            ps = sparse(point, x)
            assert full_jacobi(ps) == s.zeros(350,1)
            columns = []
            for ij, exponents in COLUMNS:
                variation = {key:dict(value) for key,value in ps.items()}
                variation.setdefault(ij,{})[exponents] = variation.get(ij,{}).get(exponents,0)+1
                columns.append(full_jacobi(variation))
            D = s.Matrix.hstack(*columns)
            Xcoord = s.Matrix(x)
            vertical, orbit = [], []
            for i,j in itertools.product(range(5),repeat=2):
                field = s.zeros(5,1)
                field[i] = x[j]
                vertical.append(vector(sparse(Xcoord*field.T-field*Xcoord.T,x)))
                E = s.zeros(5)
                E[i,j] = 1
                orbit.append(vector(sparse(E*point+point*E.T-x[j]*point.diff(x[i]),x)))
            F = s.Matrix.hstack(*vertical,*orbit,vector(ps))
            assert D*F == s.zeros(350,F.cols)
            ranks = {label:len(DomainMatrix.from_Matrix(mat).convert_to(s.QQ).rref()[1])
                     for label,mat in (("D",D),("orbit_vertical_scale",F))}
            ambient_reports.append({"source":name,"point_parameter":1,"representative_scale":str(scale),"ranks_QQ":ranks,
                "orbit_dimension":ranks["orbit_vertical_scale"]-25,
                "projective_tangent_dimension":125-ranks["D"],
                "minors":{label:minor_certificate(mat,ranks[label])
                          for label,mat in (("D",D),("orbit_vertical_scale",F))}})
    source = Path(__file__).read_text()
    return {
        "scope": "complete affine linear invariant-lift source; all residual Jacobi retained",
        "case": "smooth q, normal Hamiltonian diag(0,1,-1,2,-2)",
        "coefficient_field": "QQ",
        "q": str(q), "weights": list(weights),
        "invariant_quadratic_unknowns": len(unknowns),
        "linear_matrix_shape": list(M.shape), "linear_rank": M.rank(),
        "parameters": list(map(str, parameters)),
        "brackets": {f"{i}{j}": str(P[i, j]) for i, j in pairs},
        "residual_jacobi_coefficients": list(map(str, jacobi)),
        "equivalent_five_conditions": list(map(str, conditions)),
        "condition_span_rank":5,
        "six_complete_sources":branch_reports,
        "four_nonzero_orbit_ambient_tests":ambient_reports,
        "signed_pfaffian_minors": list(map(str, pfaffians)),
        "elapsed_seconds": time.monotonic() - started,
        "producer_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "producer_source": source,
    }


def run_degenerate_regular_support():
    """All-source SNC limits and boundary pencils for three regular types."""
    started=time.monotonic()
    root=Path(__file__).resolve().parents[1]
    x=s.symbols('x0:5');n=s.symbols('n0:6');t,z=s.symbols('t z')
    pairs=list(itertools.combinations(range(5),2))
    outputs=[]
    for name in ['Q4_REGULAR_RATIO2','Q3_REGULAR_A1_C2','Q3_REGULAR_A2_C1']:
        path=root/f'reports/CLASSIFICATION_DEGENERATE_QUADRIC_{name}_20261107.json'
        raw=json.loads(path.read_text());params=n[:len(raw['parameters'])]
        local={str(v):v for v in x+params}
        P=s.zeros(5)
        for i,j in pairs:
            P[i,j]=s.sympify(raw['brackets'][f'{i}{j}'],locals=local);P[j,i]=-P[i,j]
        jac=[s.sympify(f,locals=local) for f in raw['residual_jacobi_coefficients']]
        if name.startswith('Q4'):
            cond=[(n[2]+4)*(n[4]+2),n[0]*(n[2]-n[3]+1),n[1]*(n[3]+n[4])]
            arms=[({n[2]:-4},(1,1,0,0,3),2*(n[4]+2)),
                  ({n[4]:-2},(0,0,1,1,3),n[2]+4)]
            boundaries=[{n[0]:0,n[1]:0,n[2]:-4,n[3]:z,n[4]:-2+t},
                        {n[1]:0,n[2]:-4,n[3]:-3,n[4]:-2+t},
                        {n[0]:0,n[2]:-4+t,n[3]:2,n[4]:-2}]
        elif name.endswith('A1_C2'):
            total=n[4]+n[5]
            cond=[n[0]*(n[5]+3),n[1]*(n[4]+1),-2*n[0]*n[1]+n[2]*(total-1),5*n[2]-(3*total+2)*n[3]]
            arms=[({},(1,1,0,3,3),-2*total-3)]
            boundaries=[{n[1]:0,n[2]:0,n[3]:0,n[4]:s.Rational(3,2)+t,n[5]:-3},
                        {n[0]:0,n[2]:0,n[3]:0,n[4]:-1,n[5]:-s.Rational(1,2)+t},
                        {n[0]:0,n[1]:0,n[2]:0,n[3]:0,n[4]:z,n[5]:-s.Rational(3,2)-z+t}]
        else:
            total=n[4]+n[5]
            cond=[n[0]*n[1],n[0]*(n[5]+2),n[1]*(n[4]+1),n[2]*(total-2),(3*total+4)*n[3]-10*n[2]]
            arms=[({n[1]:0},(1,1,0,4,2),-total-3),
                  ({n[0]:0},(1,1,0,2,4),-total-3)]
            boundaries=[{n[1]:0,n[2]:0,n[3]:0,n[4]:-1+t,n[5]:-2},
                        {n[0]:0,n[2]:0,n[3]:0,n[4]:-1,n[5]:-2+t},
                        {n[0]:0,n[1]:0,n[2]:0,n[3]:0,n[4]:z,n[5]:-3-z+t}]
        polynomial_rows=[s.Poly(f,*params).as_dict() for f in jac+cond]
        monomials=sorted(set().union(*(r.keys() for r in polynomial_rows)))
        M=s.Matrix([[r.get(m,0) for m in monomials] for r in polynomial_rows])
        ranks=[M[:len(jac),:].rank(),M[len(jac):,:].rank(),M.rank()]
        assert ranks==[len(cond)]*3
        arm_reports=[]
        for sub,weights,expected in arms:
            PP=P.subs(sub,simultaneous=True).applyfunc(s.expand);limit=s.zeros(5)
            for i,j in pairs:
                for powers,coeff in s.Poly(PP[i,j],*x).terms():
                    exponent=weights[i]+weights[j]-sum(powers[k]*weights[k] for k in range(5))
                    assert exponent>=0,(name,(i,j),powers,exponent)
                    if exponent==0:limit[i,j]+=coeff*s.prod(x[k]**powers[k] for k in range(5))
                limit[j,i]=-limit[i,j]
                assert s.cancel(limit[i,j]/(x[i]*x[j])).free_symbols.isdisjoint(x)
            pf=[]
            for missing in range(5):
                i,j,k,l=[r for r in range(5) if r!=missing]
                pf.append(s.expand((-1)**missing*(limit[i,j]*limit[k,l]-limit[i,k]*limit[j,l]+limit[i,l]*limit[j,k])))
            actual=s.factor(sum(x[i]*pf[i] for i in range(5))/5)
            assert all(s.expand(s.diff(actual,x[i])-pf[i])==0 for i in range(5))
            assert s.expand(actual-expected*s.prod(x))==0,(name,actual,expected)
            arm_reports.append({'substitution':{str(k):str(v) for k,v in sub.items()},'weights':weights,'all_weights_nonnegative':True,'actual_limit':str(actual),'SNC_condition':str(expected)+' != 0'})
        for sub in boundaries:
            assert all(s.expand(f.subs(sub,simultaneous=True))==0 for f in jac)
        outputs.append({'name':name,'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'full_jacobi_equivalent_conditions':list(map(str,cond)),'coefficient_span_ranks':ranks,'SNC_arms':arm_reports,'boundary_actual_pencils':[{str(k):str(v) for k,v in sub.items()} for sub in boundaries],'all_boundary_jacobi_zero':True})
    source=Path(__file__).read_text()
    return {'scope':'three regular degenerate-quadric types: full support MPS via all-parameter SNC limits and complete boundary pencils','cases':outputs,'elapsed_seconds':time.monotonic()-started,'producer_sha256':hashlib.sha256(source.encode()).hexdigest(),'producer_source':source}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--ambient", action="store_true")
    parser.add_argument("--case", choices=("regular","rank2","rank2mixed","nil5","nil311","nil221","fixed","degenerate-regular-support"), default="regular")
    parser.add_argument("--cofactor-input")
    parser.add_argument("--cofactor-name")
    args = parser.parse_args()
    if args.case == 'degenerate-regular-support':
        result=run_degenerate_regular_support()
    elif args.case == 'fixed':
        if not args.cofactor_input or not args.cofactor_name:
            parser.error('fixed case requires --cofactor-input and --cofactor-name')
        table=json.loads(Path(args.cofactor_input).read_text())
        rows=[r for r in table['cases'] if r['name']==args.cofactor_name]
        assert len(rows)==1
        result=run_nil(args.cofactor_name,rows[0])
        result['cofactor_input_sha256']=hashlib.sha256(Path(args.cofactor_input).read_bytes()).hexdigest()
    elif args.case.startswith('nil'):
        result=run_nil(args.case)
    else:
        result = run(ambient=args.ambient) if args.case == "regular" else run_rank2(mixed=args.case == "rank2mixed")
    Path(args.output).write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in {"producer_source", "brackets", "signed_pfaffian_minors",
                                   "residual_jacobi_coefficients", "four_nonzero_orbit_ambient_tests"}},
                     ensure_ascii=False, indent=2))
    if args.ambient:
        print(json.dumps([{k:v for k,v in report.items() if k != "minors"}
                          for report in result["four_nonzero_orbit_ambient_tests"]], indent=2))

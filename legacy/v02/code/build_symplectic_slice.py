#!/usr/bin/env python3
"""Exact universal transverse equations for generic-rank-four P4 Poisson.

The constant affine symplectic point is x0^2*(d1^d2+d3^d4).
Fix its value at [e0], remove four projective perspective directions, and
eliminate 55 linear Jacobi conditions using the positive affine-degree grading.
All arithmetic and polynomial elimination is over Q.  This constructs an
equivalent 61-variable scheme, not a prime decomposition.
"""
from collections import defaultdict, Counter
from fractions import Fraction as F
from pathlib import Path
import itertools
import json
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HERE.parent/'data'
import numpy as np
import poisson_p4_components_check as p
import verify_rank2_components as v
from build_universal_slice import rref, trans, clean, add_to, substitute, encode, polynomial_row_basis

W = (4,-1,-1,-1,-1)
X0S = (2,0,0,0,0)
P0 = {(1,2):{X0S:1}, (3,4):{X0S:1}}

def weight(raw):
    pair,alpha = p.BIVECTOR_BASIS[raw]
    n=10+sum(W[i] for i in pair)-sum(a*w for a,w in zip(alpha,W))
    assert n%5==0
    return n//5

def equations_from_raw(raw):
    """Expand affine Schouten and E wedge, keeping parameter polynomials."""
    entries=defaultdict(list)
    for pair,poly in P0.items():
        for alpha,c in poly.items():
            entries[pair].append((alpha,(),F(c)))
            entries[pair[::-1]].append((alpha,(),F(-c)))
    for k,r in enumerate(raw):
        pair,alpha=p.BIVECTOR_BASIS[r]
        entries[pair].append((alpha,(k,),F(1)))
        entries[pair[::-1]].append((alpha,(k,),F(-1)))
    out=[defaultdict(F) for _ in p.FOUR_VECTOR_BASIS]
    for triple in p.TRIPLES:
        for a,b,c in [triple,triple[1:]+triple[:1],triple[2:]+triple[:2]]:
            for h in range(5):
                for alpha,ma,ca in entries[a,h]:
                    for beta,mb,cb in entries[b,c]:
                        if beta[h]==0:continue
                        gamma=[x+y for x,y in zip(alpha,beta)]
                        gamma[h]-=1
                        monomial=tuple(sorted(ma+mb))
                        for e in range(5):
                            sign=p.permutation_sign((e,)+triple)
                            if not sign:continue
                            delta=gamma.copy();delta[e]+=1
                            row=p.FOUR_VECTOR_INDEX[(tuple(sorted((e,)+triple)),tuple(delta))]
                            out[row][monomial]+=sign*ca*cb*beta[h]
    return [clean(q) for q in out]

def main():
    ROOT.mkdir(parents=True,exist_ok=True)
    t0=time.time()
    vertical=p.vertical_bivector_matrix()
    vr,vp,_=rref(vertical.T.tolist())
    active=[k for k in range(150) if k not in vp]
    def quotient(vec):
        return [F(int(vec[k]))-sum(vr[a][k]*int(vec[h]) for a,h in enumerate(vp)) for k in active]
    constants=[k for k in active if weight(k)==0]
    assert len(constants)==6
    perspectives=[]
    for j in range(1,5):
        a=np.zeros((5,5),dtype=np.int64);a[0,j]=1
        perspectives.append(quotient(p.bivector_vector(v.linear_action_derivative(P0,a))))
    _,op,_=rref(perspectives)
    orbit_raw=[active[k] for k in op]
    assert len(orbit_raw)==4 and all(weight(k)==1 for k in orbit_raw)
    raw=[k for k in active if k not in constants and k not in orbit_raw]
    assert len(raw)==116
    wt=[weight(k) for k in raw]
    Jfull=p.linearized_jacobi_matrix(P0)
    J=Jfull[:,raw].tolist()
    _,indrows,_=rref(trans(J))
    jr,jp,tr=rref([J[i] for i in indrows],transform=True)
    free=[k for k in range(len(raw)) if k not in jp]
    assert len(jp)==55 and len(free)==61
    equations=equations_from_raw(raw)
    for i,q in enumerate(equations):
        assert not q.get((),0)
        assert all(q.get((j,),0)==J[i][j] for j in range(len(raw)))
        assert len({sum(wt[k] for k in m) for m in q})<=1
    pivot_eq=[]
    for row in tr:
        q=defaultdict(F)
        for a,c in enumerate(row):add_to(q,equations[indrows[a]],c)
        pivot_eq.append(clean(q))
    fi={k:i for i,k in enumerate(free)}
    expressions={k:{(fi[k],):F(1)} for k in free}
    for a in sorted(range(len(jp)),key=lambda a:wt[jp[a]]):
        k=jp[a];q=pivot_eq[a].copy()
        assert q.pop((k,))==1
        assert all(j in expressions for m in q for j in m),(k,wt[k])
        expressions[k]={m:-c for m,c in substitute(q,expressions).items()}
    reduced=[substitute(q,expressions) for q in equations]
    assert all(not reduced[i] for i in indrows)
    assert all(len(m)>=2 for q in reduced for m in q)
    selected,rowbasis=polynomial_row_basis(reduced)
    # Preserve sparse original rows; the selected original rows span all rows.
    generators=[reduced[i] for i in selected]
    quadratic=[{m:c for m,c in q.items() if len(m)==2} for q in generators]
    qsel,qbasis=polynomial_row_basis(quadratic)
    weights=[wt[k] for k in free]
    report={
        'base_bivector':'x0^2*(d1 wedge d2+d3 wedge d4)',
        'one_parameter_weights':W,
        'vertical_pivots':vp,
        'constant_raw_indices':constants,
        'perspective_orbit_pivots_raw':orbit_raw,
        'linear_jacobi_rank_exact':len(jp),
        'affine_tangent_dimension':71,
        'projective_tangent_dimension':70,
        'projective_orbit_dimension':9,
        'transverse_tangent_dimension':61,
        'slice_ambient_count':len(raw),
        'slice_ambient_weights':dict(sorted(Counter(wt).items())),
        'eliminated_weight_counts':dict(sorted(Counter(wt[k] for k in jp).items())),
        'free_weight_counts':dict(sorted(Counter(weights).items())),
        'raw_slice_indices':raw,
        'raw_free_indices':[raw[k] for k in free],
        'free_variable_basis':[{'pair':p.BIVECTOR_BASIS[raw[k]][0],
                              'exponents':p.BIVECTOR_BASIS[raw[k]][1],
                              'weight':wt[k]} for k in free],
        'elimination_rows':indrows,
        'remaining_equation_rows':selected,
        'generator_count':len(generators),
        'generator_max_ordinary_degree':max(len(m) for q in generators for m in q),
        'generator_weight_counts':dict(sorted(Counter(next(iter({sum(weights[k] for k in m) for m in q})) for q in generators).items())),
        'generator_terms':sum(len(q) for q in generators),
        'independent_quadratic_obstructions':len(qbasis),
        'elimination_expression_max_degree':max(len(m) for k in jp for m in expressions[k]),
        'runtime_seconds':round(time.time()-t0,3),
        'claim_scope':'Exact transverse equations and component-surjective parameterization; no prime decomposition.'
    }
    data={
        'parameter_count':61,
        'parameter_weights':weights,
        'eliminated_raw_variables':[{'raw_index':raw[k],'polynomial':encode(expressions[k])} for k in jp],
        'equations':[encode(q) for q in generators],
        'quadratic_obstruction_basis':[encode(q) for q in qbasis],
    }
    (ROOT/'symplectic_slice_certificate.json').write_text(json.dumps(report,indent=2)+'\n')
    (ROOT/'symplectic_slice_equations.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    print(json.dumps({k:val for k,val in report.items() if not isinstance(val,list)},indent=2))

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Exact coefficient-field S61 normalization, avoiding expression-tree expansion.

This implements the same finite perspective action as s61_model, with arithmetic
in QQ(parameters). It does not modify the existing normalizer or certificates.
"""
import itertools
from functools import lru_cache
import sympy as s
from sympy.polys.rings import ring

from s61_model import S61,X,J,PAIRS,coefficients,pushforward,vertical_rref

MON_IDS=list(itertools.combinations_with_replacement(range(5),2))
EXPONENTS=[tuple(ids.count(i) for i in range(5)) for ids in MON_IDS]


def normalized_field(P,parameters):
    model=S61()
    K=s.QQ.frac_field(*parameters)
    R,*x=ring(','.join(map(str,X)),K)
    def from_expression(matrix):
        return {pair:R.from_expr(matrix[pair[0],pair[1]]) for pair in PAIRS}
    def entries(polynomials):
        return [polynomials[pair].get(exp,K.zero) for pair in PAIRS for exp in EXPONENTS]
    reduced,pivots=vertical_rref()
    def quotient(values):
        return [value-sum((K.from_sympy(reduced[i,k])*values[p] for i,p in enumerate(pivots)
                          if reduced[i,k]),K.zero) for k,value in enumerate(values)]
    # First move the all-ones projective point to e0 and normalize its value.
    move=s.eye(5)
    for i in range(1,5):move[i,0]=-1
    shifted=pushforward(P,move)
    point={X[0]:1,**dict.fromkeys(X[1:],0)}
    C=shifted[1:5,1:5].subs(point)
    pivot=next((i,j) for i,j in itertools.combinations(range(4),2) if C[i,j]!=0)
    perm=[*pivot,*[i for i in range(4) if i not in pivot]]
    reorder=s.eye(4)[perm,:];changed=reorder*C*reorder.T
    a,b,c,d,e,f=[changed[i,j] for i,j in itertools.combinations(range(4),2)]
    delta=s.factor(a*f-b*e+c*d)
    F=s.Matrix([[1,0,0,0],[0,1/a,0,0],[d/a,-b/a,1,0],
                [e/delta,-c/delta,0,a/delta]])*reorder
    F=F.applyfunc(s.cancel)
    assert (F*C*F.T-J[1:5,1:5]).applyfunc(s.cancel)==s.zeros(4)
    homogeneous=s.eye(5);homogeneous[1:5,1:5]=F
    based=pushforward(shifted,homogeneous)
    polynomials=from_expression(based)
    initial=quotient(entries(polynomials))
    selected=[initial[i] for i in model.certificate['perspective_orbit_pivots_raw']]
    inverse=model.perspective_matrix().inv()
    alpha=[-sum((K.from_sympy(inverse[i,j])*selected[j] for j in range(4)),K.zero)
           for i in range(4)]
    L=[x[0]-sum((v*x[j+1] for j,v in enumerate(alpha)),R.zero),*x[1:]]
    monomials={exp:L[i]*L[j] for exp,(i,j) in zip(EXPONENTS,MON_IDS)}
    pulled={pair:sum((c*monomials[mon] for mon,c in polynomial.items()),R.zero)
            for pair,polynomial in polynomials.items()}
    def entry(i,j):
        if i<j:return pulled[(i,j)]
        if i>j:return -pulled[(j,i)]
        return R.zero
    transformed=dict(pulled)
    for j in range(1,5):
        transformed[(0,j)]=entry(0,j)+sum((alpha[h-1]*entry(h,j) for h in range(1,5)),R.zero)
    normal=quotient(entries(transformed))
    assert all(not normal[i] for i in model.certificate['perspective_orbit_pivots_raw'])
    base=coefficients(X[0]**2*J)
    assert all(normal[i]==K.from_sympy(base[i]) for i in model.certificate['constant_raw_indices'])
    values=[normal[i] for i in model.certificate['raw_free_indices']]
    metadata=dict(point=[1]*5,pivot_indices=list(pivot),value_pivot=str(s.factor(a)),
                  value_pfaffian=str(delta),frame=[[str(v) for v in row] for row in F.tolist()],
                  perspective_parameters=[str(v) for v in alpha])
    return model,K,values,normal,metadata


def verify_field_map(model,K,parameters,normal):
    # Respect the actual parameter weights: using a q denominator for every
    # beta variable would unnecessarily double the polynomial degrees.
    weights=model.reduced['parameter_weights']
    D=K.one.denom
    for value,weight in zip(parameters,weights):
        if weight==1:D=D.lcm(value.denom)
    numerators=[]
    for value,weight in zip(parameters,weights):
        quotient,remainder=divmod(D**weight,value.denom)
        assert not remainder,'The selected common denominator is not weight compatible'
        numerators.append(value.numer*quotient)
    @lru_cache(None)
    def monomial(ids):
        if not ids:return D.ring.one
        return numerators[ids[0]]*monomial(ids[1:])
    # The 134 equations are weighted homogeneous, while the reconstruction
    # entries have weights <=3. No expansion of D^4 or D^5 is needed.
    powers=[D**i for i in range(4)]
    def evaluate(polynomial):
        degree=max((sum(weights[i] for i in mon) for mon,_,_ in polynomial),default=0)
        value=D.ring.zero
        for mon,numerator,denominator in polynomial:
            monomial_weight=sum(weights[i] for i in mon)
            assert monomial_weight==degree,'Unexpected nonhomogeneous S61 input'
            value+=s.QQ(numerator,denominator)*monomial(tuple(mon))*powers[degree-monomial_weight]
        return value,degree
    for i,equation in enumerate(model.reduced['equations']):
        value,_=evaluate(equation)
        assert not value,('reduced equation',i)
    assigned=set(model.certificate['raw_free_indices'])
    for item in model.data['eliminated_raw_variables']:
        value,degree=evaluate(item['polynomial'])
        index=item['raw_index'];expected=normal[index]
        assert value*expected.denom==expected.numer*powers[degree],('reconstruction',index)
        assigned.add(index)
    base=coefficients(X[0]**2*J)
    assert all(normal[i]==K.from_sympy(base[i]) for i in range(150) if i not in assigned)
    return dict(reduced_equations=134,reconstruction_coefficients=150,
                arithmetic='exact common-denominator polynomial arithmetic over QQ(parameters)',
                common_denominator_total_degree=max(map(sum,D.keys()),default=0))

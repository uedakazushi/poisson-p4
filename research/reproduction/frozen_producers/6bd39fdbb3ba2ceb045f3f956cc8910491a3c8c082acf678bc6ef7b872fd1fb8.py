#!/usr/bin/env python3
"""Complete affine normal forms for nil211 with long-top translation.

This is a finite reduction of full source type10, not a component theorem.
No producer or CAS imports: exact sparse arithmetic and Fraction matrices.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

from classification_minimal_quartic_review import LaurentRing
from classification_quadratic_cofactor_review import MONS, Ring

ROOT = Path(__file__).resolve().parents[1]
PAIRS = tuple(itertools.combinations(range(4),2))
PARAMETERS = ('alpha','A00','A10','A01','A20','A11','A02','d0','dz','dw','e0','ez',
              'delta','beta','gamma','eta','g10','g01','g20','g11','g02','g30','g21','g12','g03')


def matrix(ring, entries):
    P = [[{} for _ in range(4)] for _ in range(4)]
    for (i,j), value in zip(PAIRS,entries):
        P[i][j],P[j][i] = value,ring.scale(value,-1)
    return P


def jacobi(ring,P):
    return [ring.add(*(ring.mul(P[a][k],ring.derivative(P[b][d],k))
                        for a,b,d in [(i,j,l),(j,l,i),(l,i,j)] for k in range(4)))
            for i,j,l in itertools.combinations(range(4),3)]


def divergence(ring,P):
    return [ring.add(*(ring.derivative(P[i][j],i) for i in range(4))) for j in range(4)]


def pfaffian(ring,entries):
    return ring.add(ring.mul(entries[0],entries[5]),ring.scale(ring.mul(entries[1],entries[4]),-1),
                    ring.mul(entries[2],entries[3]))


def generic_data(ring,r,s,original_coordinates):
    A = ring.parse('A00+A10*z+A01*w+A20*z**2+A11*z*w+A02*w**2')
    d = ring.parse('d0+dz*z+dw*w')
    e = ring.parse('e0+ez*z+(alpha-1-dz)*w')
    L = ring.parse('eta+beta*w-gamma*z')
    G = ring.parse('g10*z+g01*w+g20*z**2+g11*z*w+g02*w**2+g30*z**3+g21*z**2*w+g12*z*w**2+g03*w**3')
    a = ring.add(ring.mul(ring.variable('alpha'),r),A)
    phi = ring.add(ring.scale(ring.mul(ring.variable('delta'),r,r),Q(1,2)),ring.mul(r,L),G)
    p = ring.add(ring.mul(ring.variable('beta'),r),ring.derivative(G,3))
    q = ring.add(ring.mul(ring.variable('gamma'),r),ring.scale(ring.derivative(G,2),-1))
    entries = [a,ring.add(p,ring.mul(s,d)) if original_coordinates else p,
               ring.add(q,ring.mul(s,e)) if original_coordinates else q,d,e,
               ring.add(ring.mul(ring.variable('delta'),r),L)]
    return dict(A=A,d=d,e=e,L=L,G=G,a=a,phi=phi,entries=entries)


def inverse(matrix_):
    n = len(matrix_)
    A = [[Q(value) for value in row]+[Q(i==j) for j in range(n)] for i,row in enumerate(matrix_)]
    determinant = Q(1)
    for j in range(n):
        pivot = next(i for i in range(j,n) if A[i][j])
        if pivot != j:
            A[pivot],A[j] = A[j],A[pivot]
            determinant *= -1
        value = A[j][j]
        determinant *= value
        A[j] = [t/value for t in A[j]]
        for i in range(n):
            if i != j and A[i][j]:
                value = A[i][j]
                A[i] = [a-value*b for a,b in zip(A[i],A[j])]
    assert [row[:n] for row in A] == [[Q(i==j) for j in range(n)] for i in range(n)]
    return [row[n:] for row in A],determinant


def affine_vectors(ring,entries):
    result = []
    for value in entries:
        for monomial in MONS:
            vector = [Q(0)]*(ring.n+1)
            for powers,coefficient in value.items():
                if powers[-4:] == monomial:
                    parameter_powers = powers[:-4]
                    assert min(parameter_powers)>=0 and sum(parameter_powers)<=1
                    index = 0 if not any(parameter_powers) else parameter_powers.index(1)+1
                    vector[index] += coefficient
            result.append(vector)
    return result


def substitute_constant(ring,polynomial,name,value):
    position = ring.names.index(name)
    result = {}
    for powers,coefficient in polynomial.items():
        coefficient *= value**powers[position]
        target = list(powers)
        target[position] = 0
        target = tuple(target)
        result[target] = result.get(target,Q(0))+coefficient
    return {e:c for e,c in result.items() if c}


def group_r_power(ring,polynomial,degree):
    result = {}
    for powers,coefficient in polynomial.items():
        if powers[ring.n] == degree:
            target = list(powers)
            target[ring.n] = 0
            result[tuple(target)] = coefficient
    return result


def verify_complete_source(source):
    ring = LaurentRing(PARAMETERS,('x','y','z','w'))
    r,s = ring.parse('x-y**2/2'),ring.variable('y')
    data = generic_data(ring,r,s,True)
    entries = data['entries']
    assert all(sum(e[-4:])<=2 for value in entries for e in value)
    P = matrix(ring,entries)
    assert divergence(ring,P) == [ring.variable('y'),ring.constant(1),{},{}]
    lie = []
    for i,j in PAIRS:
        value = ring.add(ring.mul(ring.variable('y'),ring.derivative(P[i][j],0)),ring.derivative(P[i][j],1))
        if i == 0:
            value = ring.add(value,ring.scale(P[1][j],-1))
        if j == 0:
            value = ring.add(value,ring.scale(P[i][1],-1))
        lie.append(value)
    assert lie == [{}]*6
    vectors = affine_vectors(ring,entries)
    free = source['free_coefficient_indices']
    assert len(free) == source['free_parameter_count'] == 25
    selected = [vectors[i][1:] for i in free]
    offsets = [vectors[i][0] for i in free]
    backwards,determinant = inverse(selected)
    old_ring = Ring(25)
    # Transport the entire ninety-coefficient affine solution via its free
    # coordinates. The inverse map is explicit, not a dimension-only claim.
    for actual,reported in zip(vectors,source['all_ninety_coefficients']):
        coefficients = [sum(actual[k+1]*backwards[k][j] for k in range(25)) for j in range(25)]
        constant = actual[0]-sum(a*b for a,b in zip(coefficients,offsets))
        polynomial = old_ring.constant(constant)
        for j,value in enumerate(coefficients):
            polynomial = old_ring.add(polynomial,old_ring.scale(old_ring.variable(f'n{j}'),value))
        assert polynomial == old_ring.parse(reported)
    Vphi = ring.add(ring.scale(ring.mul(data['a'],ring.derivative(data['phi'],0)),-1),
                    ring.mul(data['d'],ring.derivative(data['phi'],2)),
                    ring.mul(data['e'],ring.derivative(data['phi'],3)))
    assert pfaffian(ring,entries) == ring.scale(Vphi,-1)
    delta_zero = substitute_constant(ring,pfaffian(ring,entries),'delta',Q(0))
    assert not {e:c for e,c in delta_zero.items() if sum(e[-4:])==4}
    return dict(all_original_ninety_coefficients_matched=True,source_parameter_count=25,
                all_original_degree_two_constraints_preserved=True,
                full_divergence=['y','1','0','0'],full_cofactor_invariance_checked=True,
                free_indices=free,affine_parameter_matrix_determinant=str(determinant),
                inverse_parameter_matrix=[[str(a) for a in row] for row in backwards],
                original_free_coordinate_offsets=list(map(str,offsets)),
                raw_pfaffian_is_minus_V_phi=True,delta_zero_implies_highest_pfaffian_zero=True)


def verify_jacobi_and_square_completion():
    ring = LaurentRing(PARAMETERS,('r','s','z','w'),invert=('delta',))
    r,s = ring.variable('r'),ring.variable('s')
    data = generic_data(ring,r,s,False)
    P = matrix(ring,data['entries'])
    Vphi = ring.add(ring.scale(ring.mul(data['a'],ring.derivative(data['phi'],0)),-1),
                    ring.mul(data['d'],ring.derivative(data['phi'],2)),
                    ring.mul(data['e'],ring.derivative(data['phi'],3)))
    psi = ring.add(Vphi,data['phi'])
    assert jacobi(ring,P) == [ring.scale(ring.derivative(psi,3),-1),ring.derivative(psi,2),{},ring.derivative(psi,0)]
    assert group_r_power(ring,psi,2) == ring.parse('delta*(1/2-alpha)')
    half = lambda polynomial: substitute_constant(ring,polynomial,'alpha',Q(1,2))
    d,e = data['d'],half(data['e'])
    base_L = ring.add(ring.mul(d,ring.derivative(data['L'],2)),ring.mul(e,ring.derivative(data['L'],3)))
    assert group_r_power(ring,half(psi),1) == ring.add(ring.scale(data['L'],Q(1,2)),base_L,
                                                                      ring.scale(ring.mul(ring.variable('delta'),data['A']),-1))
    A = ring.mul(ring.parse('1/delta'),ring.add(ring.scale(data['L'],Q(1,2)),base_L))
    assert all(sum(exponents[-4:])<=1 for exponents in A)
    R = ring.add(r,ring.mul(data['L'],ring.parse('1/delta')))
    H = ring.add(data['G'],ring.scale(ring.mul(data['L'],data['L'],ring.parse('1/delta')),Q(-1,2)))
    Vr = ring.add(ring.scale(r,Q(-1,2)),ring.scale(A,-1))
    assert ring.add(Vr,ring.mul(base_L,ring.parse('1/delta'))) == ring.scale(R,Q(-1,2))
    assert ring.add(data['phi'],ring.scale(ring.mul(ring.variable('delta'),R,R),Q(-1,2))) == H
    normalized_psi = ring.add(ring.mul(Vr,ring.derivative(data['phi'],0)),
                               ring.mul(d,ring.derivative(data['phi'],2)),
                               ring.mul(e,ring.derivative(data['phi'],3)),data['phi'])
    assert normalized_psi == ring.add(ring.mul(d,ring.derivative(H,2)),ring.mul(e,ring.derivative(H,3)),H)
    return dict(full_Jacobi_equivalent_to_psi_constant=True,
                psi_definition='V(phi)+phi',r_squared_coefficient='delta*(1/2-alpha)',
                delta_nonzero_forces_alpha='1/2',
                r_coefficient='L/2+Vbase(L)-delta*A2',A2_is_affine=True,
                R='r+L/delta',H='G-L^2/(2delta)',V_R='-R/2',
                complete_square_identity_checked=True,
                remaining_equation='Vbase(H)+H=constant, trace(linear Vbase)=-1/2')


def weighted_table():
    monomials = [(i,j) for i in range(4) for j in range(4-i) if i+j]
    spectra = set()
    for i,j in monomials:
        if i != j:
            lam = Q(j-2,2*(i-j))
            spectra.add(tuple(sorted((lam,Q(-1,2)-lam))))
        else:
            assert -Q(i,2) != -1
    assert spectra == {(Q(-1),Q(1,2)),(Q(-1,2),Q(0)),(Q(-1,3),Q(-1,6))}
    table = []
    for lam,mu in sorted(spectra):
        allowed = [(i,j) for i,j in monomials if i*lam+j*mu == -1]
        table.append(dict(eigenvalues=[str(lam),str(mu)],all_degree_at_most_three_eigenmonomials=[list(t) for t in allowed]))
    assert not [(i,j) for i,j in monomials if Q(-1,4)*(i+j) == -1]
    zero_ring = LaurentRing(('a','b','kappa'),('z','w','p','q'))
    H = zero_ring.parse('z**2*(a+b*w)')
    assert zero_ring.mul(zero_ring.variable('kappa'),zero_ring.derivative(H,1)) == zero_ring.parse('kappa*b*z**2')
    return dict(all_trace_minus_half_spectra=table,
                repeated_eigenvalue_minus_quarter_has_no_eigenpolynomial=True,
                zero_eigen_translation_equation='kappa*b=0 for H=z^2(a+bw)',
                irreducible_F_forces_b_nonzero=True,
                resulting_normalized_H=['z','z^3','z^2w'])


def actual_models():
    ring = LaurentRing((),('x','y','z','w'))
    models = [dict(name='linear',H='z',lam='-1',mu='1/2',
                   brackets=['x/2-y**2/4','-y*z','-1+y*w/2','-z','w/2','x-y**2/2']),
              dict(name='cubic',H='z**3',lam='-1/3',mu='-1/6',
                   brackets=['x/2-y**2/4','-y*z/3','-3*z**2-y*w/6','-z/3','-w/6','x-y**2/2']),
              dict(name='mixed',H='z**2*w',lam='-1/2',mu='0',
                   brackets=['x/2-y**2/4','z**2-y*z/2','-2*z*w','-z/2','0','x-y**2/2'])]
    for model in models:
        entries = [ring.parse(t) for t in model['brackets']]
        P = matrix(ring,entries)
        assert jacobi(ring,P) == [{}]*4
        assert divergence(ring,P) == [ring.variable('y'),ring.constant(1),{},{}]
        F = ring.parse('(x-y**2/2)**2/2+'+model['H'])
        assert pfaffian(ring,entries) == F
        top = {e:c for e,c in F.items() if sum(e)==4}
        assert top == ring.parse('y**4/8')
        model.update(all_four_Jacobi_zero=True,full_divergence=['y','1','0','0'],
                     affine_pfaffian='(x-y^2/2)^2/2+'+model['H'],highest_pfaffian='y^4/8',
                     no_component_membership_or_componenthood_claim=True)
    return models


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(170,175))
    resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    start = time.monotonic()
    input_path = ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json'
    source = json.loads(input_path.read_text())['full_nilpotent_sources'][9]
    assert source['variant'] == dict(jordan_blocks=[2,1,1],translation_top_length=2)
    bridge = verify_complete_source(source)
    reduction = verify_jacobi_and_square_completion()
    table = weighted_table()
    models = actual_models()
    inputs = [Path(__file__),Path(__file__).with_name('classification_minimal_quartic_review.py'),
              Path(__file__).with_name('classification_quadratic_cofactor_review.py'),input_path]
    result = dict(task_id='E_reduced_nonnormal',all_passed=True,no_producer_import=True,no_CAS_import=True,
                  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                  complete_original_source_bridge=bridge,full_Jacobi_reduction=reduction,
                  complete_weighted_table=table,three_actual_Poisson_models=models,
                  normalizing_coordinate_changes='r=x-y^2/2 is auxiliary only; r->r+L/delta is the actual affine x->x+L/delta. Base diagonalization/translations are affine. Final z->tz,w->uw gives delta_new=delta*t*u, b_new=b*t^(1-i)*u^(1-j).',
                  scale_normalization={'H=bz':['t=b/delta','u=1/b'],
                                       'H=bz^3':['choose t^3=b/delta','u=1/(delta*t)'],
                                       'H=bz^2w':['t=b','u=1/(delta*b)']},
                  inverted_elements=['delta!=0 forced by simpleH','nonzero H coefficientb after excluding reducibleF','nonzero determinant base affine coordinate transforms'],
                  complements=['delta=0 forces highestPfaffian0, hence H is not simple','repeated eigenvalues -1/4 have no degree<=3 eigenpolynomial','zero-eigen translation with kappa!=0 forces b=0 and a reducibleF','H coefficientb=0 makesF reducible'],
                  scope='All type10 simpleH/irreducible-quartic points reduce to three actual affine/PGL orbits. Ambient componenthood, known-component membership, pairwise/global distinctness remain open.',
                  seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  source_bytes=Path(__file__).read_text())
    with output.open('x') as handle:
        json.dump(result,handle,indent=2)
        handle.write('\n')
    print(json.dumps({key:result[key] for key in ('all_passed','seconds','maximum_RSS_KiB')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)

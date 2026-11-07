#!/usr/bin/env python3
"""Independent homogeneous-coordinate review of the full cubic collision.

No producer imports. Start with an arbitrary ten-coefficient cubic, compute
the inverse five-coordinate linear map, and compare its regular polynomial
limit with the immutable full six-magnetic source.
"""
import argparse
import hashlib
import itertools
import json
import resource
import time
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cubic_tensor(f, lam, variables):
    H, A, B, C, W = variables
    result = sp.zeros(5)
    entries = {(1, 2): sp.diff(f, C), (1, 3): -sp.diff(f, B),
               (2, 3): sp.diff(f, A)}
    entries.update({(i, 4): -lam*W*variables[i] for i in range(1, 4)})
    for (i, j), value in entries.items():
        result[i, j], result[j, i] = value, -value
    return result


def strict_jacobi(tensor, variables):
    values = []
    for a, b, c in itertools.combinations(range(5), 3):
        value = sum(tensor[a, k]*sp.diff(tensor[b, c], variables[k])
                    + tensor[b, k]*sp.diff(tensor[c, a], variables[k])
                    + tensor[c, k]*sp.diff(tensor[a, b], variables[k])
                    for k in range(5))
        values.append(sp.expand(value))
    return values


def pfaffian(tensor):
    return sp.expand(tensor[1, 2]*tensor[3, 4]
                     - tensor[1, 3]*tensor[2, 4]
                     + tensor[1, 4]*tensor[2, 3])


def bracket_strings(tensor):
    return {f'{i},{j}': str(tensor[i, j])
            for i, j in itertools.combinations(range(1, 5), 2)}


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU, (100, 105))
    resource.setrlimit(resource.RLIMIT_AS, (2048*1024**2, 2048*1024**2))
    started = time.monotonic()
    H, A, B, C, W = old = sp.symbols('H A B C W')
    h, x, y, z, w = new = sp.symbols('h x y z w')
    e, lam = sp.symbols('epsilon lambda')
    r = sp.symbols('r0:6')
    g = sp.symbols('g0:4')
    free = sp.symbols('c0:10')
    monomials = [A**i*B**j*C**(3-i-j) for i in range(4) for j in range(4-i)]
    arbitrary_f = sum(c*m for c, m in zip(free, monomials))
    universal = cubic_tensor(arbitrary_f, lam, old)
    assert strict_jacobi(universal, old) == [0]*10
    assert sp.expand(pfaffian(universal)+3*lam*W*arbitrary_f) == 0

    # All six magnetic terms are recovered by integrating in B, with zero
    # integration constant. No hand-copied primitive from the producer.
    magnetic = (r[0]*C*C+r[1]*A*C+r[2]*B*C
                + r[3]*A*A+r[4]*A*B+r[5]*B*B)
    primitive = sp.integrate(3*magnetic, B)
    assert sp.expand(sp.diff(primitive, B)-3*magnetic) == 0
    f0 = -lam*A*B*C+sum(g[i]*A**(3-i)*C**i for i in range(4))
    f = f0+e*primitive
    tensor = cubic_tensor(f, lam, old)

    # Rows give OLD homogeneous coordinates in terms of NEW coordinates.
    transform = sp.Matrix([[1,0,0,0,0], [0,1,0,0,0], [0,0,1,0,0],
                           [1,0,0,e,0], [0,1,0,0,e]])
    assert transform.det() == e**2
    inverse = transform.inv()
    substitution = dict(zip(old, transform*sp.Matrix(new)))
    pulled = e*inverse*tensor.subs(substitution, simultaneous=True)*inverse.T
    family = pulled.applyfunc(lambda v: sp.expand(sp.cancel(v)))
    assert all(sp.denom(v) == 1 or not sp.denom(v).free_symbols for v in family)
    assert all(v == 0 or sp.Poly(v, *new).total_degree() == 2 for v in family)
    assert family.row(0) == sp.zeros(1, 5)
    assert strict_jacobi(family, new) == [0]*10
    expected_pf = -3*lam*(x+e*w)*f.subs(substitution, simultaneous=True)
    assert sp.expand(pfaffian(family)-expected_pf) == 0
    limit = family.subs(e, 0)
    assert sp.expand(limit[3, 4]+lam*h*w
                     +3*magnetic.subs({A:x, B:y, C:h}, simultaneous=True)) == 0

    producer_path = ROOT/'reports/CLASSIFICATION_PLANE_CUBIC_TORIC_COLLISION_20261107.json'
    producer = json.loads(producer_path.read_text())
    env = {str(q): q for q in (*old, *new, e, lam, *r, *g)}
    env['lam'] = lam
    parse_saved = lambda value: sp.sympify(value.replace('lambda', 'lam'), locals=env)
    for key, value in producer['full_polynomial_collision_bivector'].items():
        i, j = map(int, key.split(','))
        assert sp.expand(family[i, j]-parse_saved(value)) == 0
    assert sp.expand(primitive-parse_saved(producer['magnetic_primitive_G'])) == 0

    node_parameters = {lam:-1, **dict(zip(g, [-sp.Rational(1,8),sp.Rational(3,8),
                                           -sp.Rational(3,8),sp.Rational(1,8)]))}
    Z, V = sp.symbols('Z V')
    final_coordinates = (h, x, y, Z, V)
    fiber = sp.eye(5)
    fiber[3,3] = -sp.Rational(8,3)
    fiber[4,3] = sp.Rational(1,8)
    fiber[4,4] = sp.Rational(1,8)
    assert fiber.det() == -sp.Rational(1,3)
    back = dict(zip(new, fiber.inv()*sp.Matrix(final_coordinates)))
    final_family = (24*fiber*family.subs(node_parameters)*fiber.T).subs(back, simultaneous=True).applyfunc(sp.expand)
    final_limit = final_family.subs(e, 0)
    target_path = ROOT/'reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_20261107.json'
    target_report = json.loads(target_path.read_text())
    case = next(v for v in target_report['cases'] if v['case'] == 'node_two_node_tangents')
    target_env = {'h':h,'x':x,'y':y,'z':Z,'w':V, **{str(q):q for q in r}}
    target = sp.zeros(5)
    for key, value in case['full_six_magnetic_source'].items():
        i, j = map(int, key.strip('()').split(','))
        target[i,j] = sp.sympify(value, locals=target_env)
        target[j,i] = -target[i,j]
    assert final_limit == target
    assert strict_jacobi(final_family, final_coordinates) == [0]*10
    node_f = x*h*y-(x-h)**3/8
    assert sp.expand(pfaffian(final_limit)+576*x*node_f) == 0
    # The arbitrary magnetic parameters survive with independent coefficients:
    # each r_i appears as 24 times a different degree-two monomial in P_ZV.
    magnetic_columns = [sp.diff(final_limit[3,4], q) for q in r]
    expected_columns = [24*h*h,24*h*x,24*h*y,24*x*x,24*x*y,24*y*y]
    assert [sp.expand(v) for v in magnetic_columns] == expected_columns

    # Reconstruct the SNC witness from f=ABC, rather than copying its bracket
    # list. The H row stays Casimir and all five projective factors are the
    # independent original coordinates.
    toric = cubic_tensor(A*B*C, -1, old)
    assert strict_jacobi(toric, old) == [0]*10
    assert pfaffian(toric) == 3*A*B*C*W

    inputs = [Path(__file__), producer_path, target_path,
              ROOT/'research/classification/KNOWN_COMPONENT_AUDIT.md']
    result = dict(task_id='E_reduced_nonnormal', all_passed=True,
                  method='Independent arbitrary ten-coefficient cubic, exact inverse five-dimensional homogeneous coordinate map, full strict Jacobi and all-six-magnetic comparison; producer not imported.',
                  coefficient_ring='QQ[ten cubic coefficients,lambda,g0..g3,r0..r5,epsilon,h,x,y,z,w]',
                  input_sha256={str(p.relative_to(ROOT)):digest(p) for p in inputs},
                  arbitrary_cubic_whole_Jacobi_zero=True, arbitrary_cubic_pfaffian_identity=True,
                  homogeneous_coordinate_determinant='epsilon**2',
                  inverse_PGL_pushforward_is_polynomial=True,
                  complete_family_homogeneous_degree_two=True,
                  complete_family_strict_Jacobi_zero=True,
                  full_polynomial_family=bracket_strings(family),
                  complete_family_pfaffian_identity=True,
                  magnetic_primitive_derivative='G_B=3R',
                  limit_B_coefficient=str(limit[3,4]),
                  exact_existing_six_magnetic_source_match=True,
                  six_independent_magnetic_columns=list(map(str, magnetic_columns)),
                  node_limit_pfaffian=str(-576*x*node_f),
                  toric_brackets=bracket_strings(toric), toric_pfaffian='3*A*B*C*W',
                  toric_projective_divisor='3*H*A*B*C*W',
                  inverted_elements=['epsilon for the invertible PGL interpretation only'],
                  complements=['epsilon=0 is the actual regular limit, with all six r coefficients independent and unrestricted'],
                  geometric_review='The universal arbitrary-cubic source (lambda!=0) is irreducible and contains the rank-four SNC toric point f=ABC,lambda=-1. An irreducible component containing its image closure meets the SNC locus, hence is one of the accepted MPS40. The finite closed MPS union is PGL-invariant and closed under nonzero tensor scaling; regular specialization therefore includes every six-magnetic epsilon=0 target. Uniqueness of the component at the toric point is not required.',
                  external_dependency='MPS Theorem1.4 identifying all components meeting the SNC log-symplectic locus, as audited in KNOWN_COMPONENT_AUDIT.md section4. Not reproved here.',
                  limitation='Closes this entire node_two_node_tangents source only; does not cover the other three plane configurations or all1+1+3 supports.',
                  reviewer_source_bytes=Path(__file__).read_text(),
                  seconds=time.monotonic()-started,
                  maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k:result[k] for k in ('all_passed','seconds','maximum_RSS_KiB')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

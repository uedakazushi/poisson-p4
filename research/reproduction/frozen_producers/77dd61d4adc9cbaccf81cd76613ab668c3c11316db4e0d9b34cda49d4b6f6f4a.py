#!/usr/bin/env python3
"""Full ambient ranks at one actual eight-parameter residual-source point."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import (
    COLUMNS, jacobi, sparse, vector, minor_certificate)

ROOT = Path(__file__).resolve().parents[1]


def projective_source():
    h, x, y, z, w = xx = s.symbols('h x y z w')
    a0, a1, *rs = params = s.symbols('a0 a1 r0:6')
    f = x*(x-h)
    q = (a0*h+a1*x)*y+f
    magnetic = sum(c*m for c, m in zip(rs, (h*h, h*x, h*y, x*x, x*y, y*y)))
    P = s.zeros(5)
    for ij, value in {(1, 3): f, (2, 3): y*(a1*y+2*x-h),
                      (2, 4): q, (3, 4): a1*y*w-magnetic}.items():
        P[ij] = s.expand(value)
        P[ij[::-1]] = -P[ij]
    return P, xx, tuple(params), f, q


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU, (55, 56))
    resource.setrlimit(resource.RLIMIT_AS, (1536*1024**2, 1536*1024**2))
    begun = time.monotonic()
    source, xx, params, f, q = projective_source()
    universal_jacobi = [s.expand(sum(
        source[i, ell]*s.diff(source[j, k], xx[ell])
        + source[j, ell]*s.diff(source[k, i], xx[ell])
        + source[k, ell]*s.diff(source[i, j], xx[ell])
        for ell in range(5))) for i, j, k in itertools.combinations(range(5), 3)]
    assert universal_jacobi == [0]*10
    pf = s.expand(source[1, 2]*source[3, 4]-source[1, 3]*source[2, 4]
                  + source[1, 4]*source[2, 3])
    assert s.expand(pf+f*q) == 0
    values = dict(zip(params, (1, 1, 1, 2, 3, 5, 7, 11)))
    P = source.subs(values)
    ps = sparse(P, xx)
    assert all(not c for c in jacobi(ps, raw=True).values())
    assert jacobi(ps) == s.zeros(350, 1)
    assert (s.hessian(q.subs(values), xx)/2).rank() == 3
    directions = [source.diff(parameter).subs(values) for parameter in params]
    dcols = []
    for ij, exponents in COLUMNS:
        variation = {key: dict(value) for key, value in ps.items()}
        variation.setdefault(ij, {})[exponents] = variation.get(ij, {}).get(exponents, 0)+1
        dcols.append(jacobi(variation))
    D = s.Matrix.hstack(*dcols)
    X = s.Matrix(xx)
    vertical, orbit = [], []
    for a, b in itertools.product(range(5), repeat=2):
        field = s.zeros(5, 1)
        field[a] = xx[b]
        vertical.append(vector(sparse(X*field.T-field*X.T, xx)))
        E = s.zeros(5)
        E[a, b] = 1
        orbit.append(vector(sparse(E*P+P*E.T-xx[b]*P.diff(xx[a]), xx)))
    V = s.Matrix.hstack(*vertical)
    M = V.row_join(s.Matrix.hstack(*orbit)).row_join(vector(ps))
    F = M.row_join(s.Matrix.hstack(*(vector(sparse(R, xx)) for R in directions)))
    assert D*F == s.zeros(350, F.cols)
    mats = {'full350_D': D, 'vertical': V, 'orbit_vertical_scale': M,
            'actual_full_source': F}
    ranks = {name: len(DomainMatrix.from_Matrix(A).convert_to(s.QQ).rref()[1])
             for name, A in mats.items()}
    certificates = {name: minor_certificate(A, ranks[name]) for name, A in mats.items()}
    lower = ranks['actual_full_source']-25
    upper = 125-ranks['full350_D']
    paths = ['code/classification_three_hyperplanes_quadric_dependent.py',
             'code/classification_rank3_mixed_component_review.py',
             'research/classification/THREE_HYPERPLANES_QUADRIC_DEPENDENT.md']
    report = dict(
        task_id='E_reduced_nonnormal', all_passed=True,
        source_parameters=list(map(str, params)),
        source_open='a0*(a0+a1) != 0; omitted factor choices retained by permutation of the three hyperplanes',
        universal_all_ten_raw_Jacobi_zero=True,
        universal_affine_pfaffian=str(-f*q),
        universal_actual_quintic=str(-xx[0]*f*q),
        all_eight_actual_source_directions_retained=True,
        point_substitution={str(k): str(v) for k, v in values.items()},
        point_brackets={str(ij): str(P[ij]) for ij in itertools.combinations(range(5), 2) if P[ij]},
        all8_actual_parameter_directions=[{str(ij): str(R[ij])
            for ij in itertools.combinations(range(5), 2) if R[ij]} for R in directions],
        full350_projective_Jacobi_zero=True,
        full350_times_all_actual_directions_zero=True,
        ranks_QQ=ranks, certificates=certificates,
        projective_family_dimension_lower=lower,
        projective_tangent_dimension_upper=upper,
        true_transverse_tangent_dimension=upper-lower,
        smooth_full_component_proved=(lower == upper),
        known_component_membership_not_checked=True,
        global_exhaustion_not_claimed=True,
        input_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        matrix_sha256={name: hashlib.sha256(json.dumps(A.tolist(), default=str).encode()).hexdigest()
                       for name, A in mats.items()},
        seconds=time.monotonic()-begun)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x') as handle:
        json.dump(report, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(ranks_QQ=ranks, lower=lower, upper=upper,
                         transverse=upper-lower, seconds=report['seconds'])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

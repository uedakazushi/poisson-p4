#!/usr/bin/env python3
"""Complete quadratic W sources for the eleven nilpotent cofactor normal forms."""
import argparse
import hashlib
import itertools
import json
import resource
import time
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU, (55, 56))
    resource.setrlimit(resource.RLIMIT_AS, (1536 * 1024**2,) * 2)
    start = time.monotonic()
    h, x, u, v, w = xx = s.symbols('h x u v w')
    mons = [xx[i] * xx[j] for i, j in itertools.combinations_with_replacement(range(5), 2)]
    unknowns = s.symbols('a0:45')
    forms = [('j3_zero', (v, w, 0)), ('j3_h', (v, w, h)),
             ('j3_hx', (v, w, h + x))]
    for name, L, M in [('zero', 0, 0), ('0_h', 0, h), ('0_hx', 0, h+x),
                       ('h_0', h, 0), ('hx_0', h+x, 0),
                       ('h_x', h, x), ('hx_h', h+x, h)]:
        forms.append(('j21_' + name, (v, L, M)))
    forms.append(('j111_rank2', (x, h, 0)))
    records = []
    for name, V in forms:
        P = s.zeros(5)
        for i, value in enumerate(V, 2):
            P[1, i] = x * value
        for k, (i, j) in enumerate([(2, 3), (2, 4), (3, 4)]):
            P[i, j] = sum(a*m for a, m in zip(unknowns[15*k:15*(k+1)], mons))
        for i, j in itertools.combinations(range(5), 2):
            P[j, i] = -P[i, j]

        def jac(i, j, k, matrix):
            return s.expand(sum(matrix[i, l]*s.diff(matrix[j, k], xx[l]) +
                                matrix[j, l]*s.diff(matrix[k, i], xx[l]) +
                                matrix[k, l]*s.diff(matrix[i, j], xx[l]) for l in range(5)))

        linear = []
        for i, j in itertools.combinations(range(2, 5), 2):
            linear += s.Poly(jac(1, i, j, P), *xx).coeffs()
        matrix, rhs = s.linear_eq_to_matrix(linear, unknowns)
        aug, pivots = matrix.row_join(rhs).rref()
        if 45 in pivots:
            records.append(dict(case=name, V=[str(f) for f in V], inconsistent=True))
            print(json.dumps(dict(case=name, inconsistent=True)), flush=True)
            continue
        free = [j for j in range(45) if j not in pivots]
        params = s.symbols('n0:' + str(len(free)))
        solution = {unknowns[j]: z for j, z in zip(free, params)}
        for i, pivot in enumerate(pivots):
            solution[unknowns[pivot]] = aug[i, 45] - sum(aug[i, j]*solution[unknowns[j]] for j in free)
        Q = P.subs(solution, simultaneous=True).applyfunc(s.expand)
        assert all(jac(1, i, j, Q) == 0 for i, j in itertools.combinations(range(2, 5), 2))
        residual = jac(2, 3, 4, Q)
        equations = list(dict.fromkeys(str(s.factor(f)) for f in s.Poly(residual, *xx).coeffs() if f))
        F = s.expand(V[0]*Q[3, 4]-V[1]*Q[2, 4]+V[2]*Q[2, 3])
        assert s.expand(sum(V[i]*s.diff(F, xx[i+2]) for i in range(3))) == 0
        records.append(dict(case=name, V=[str(f) for f in V], inconsistent=False,
            linear_matrix_shape=list(matrix.shape), linear_rank=len(pivots),
            free_parameters=list(map(str, params)), dimension=len(free),
            original_unknown_monomials=[str(m) for m in mons], free_columns=free,
            linear_rref=[[str(f) for f in aug.row(i)] for i in range(len(pivots))],
            brackets={','.join(map(str, ij)): str(Q[ij]) for ij in itertools.combinations(range(1,5),2)},
            actual_cubic=str(s.factor(F)), remaining_full_Jacobi=str(residual),
            remaining_coefficient_equations=equations))
        print(json.dumps(dict(case=name, dimension=len(free), equations=len(equations), F=str(s.factor(F)))), flush=True)
    source = Path(__file__).read_text()
    result = dict(task_id='E_reduced_nonnormal', cases=records,
        all_45_vertical_coefficients_retained=True,
        all_four_affine_Jacobi_equations_retained=True,
        normal_form_exhaustion_requires_the_written_centralizer_argument=True,
        irreducibility_and_MPS_membership_not_inferred_from_this_linear_calculation=True,
        seconds=time.monotonic()-start, producer_source=source,
        input_sha256={'code/classification_two_hyperplanes_cubic_nilpotent.py': hashlib.sha256(source.encode()).hexdigest()})
    with output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

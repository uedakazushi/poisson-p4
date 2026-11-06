#!/usr/bin/env python3
"""One full ambient component test for each remaining nilJordan4 source.

The sources already exhaust this modular stratum. No sampled parameter farm
or truncated Jacobi system is used: all 150 ambient variables and 350 rows.
"""
import hashlib
import itertools
import json
from pathlib import Path
import resource
import sys
import time

import sympy as sp
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'legacy/v03/code'))
import poisson_p4_components_check as ambient


def dm(mat):
    return DomainMatrix.from_Matrix(mat).convert_to(sp.QQ)


def nonzero_minor(mat):
    columns = list(dm(mat).rref()[1])
    rows = list(dm(mat[:, columns].T).rref()[1])
    det = dm(mat[rows, columns]).det()
    assert det != 0
    return {'rank': len(columns), 'rows': rows, 'columns': columns,
            'determinant': str(det)}


def main():
    resource.setrlimit(resource.RLIMIT_CPU, (90, 90))
    resource.setrlimit(resource.RLIMIT_AS, (1536*1024**2, 1536*1024**2))
    started = time.monotonic()
    source_path = ROOT/'reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_NILPOTENT_SUPPORT_20261107.json'
    source = json.loads(source_path.read_text())
    h = sp.Symbol('h')
    xx = sp.symbols('x1:5')
    yy = (h,)+xx
    nn = sp.symbols('n0:9')
    env = {str(v): v for v in xx+nn}
    base = [sp.sympify(f, locals=env) for f in source['charts'][0]['all_six_brackets']]
    results = {}
    for name, branch in source['modular_zero_chart']['five_full_sources'].items():
        if name == 's_1':
            continue  # Already identified with the old smooth triple-hyperplane component.
        substitutions = {env[k]: sp.sympify(v, locals=env)
                         for k, v in branch['all_point_substitution'].items()}
        parameters = [env[v] for v in branch['free_parameters']]
        nonzero = env[branch['rank4_iff_nonzero']]
        point = {v: int(v == nonzero) for v in parameters}
        family = sp.zeros(5)
        for (i, j), entry in zip(itertools.combinations(range(1, 5), 2), base):
            entry = sp.expand(entry.subs(substitutions, simultaneous=True))
            family[i, j] = sum(c*h**(2-sum(mon))*sp.prod(x**e for x, e in zip(xx, mon))
                               for mon, c in sp.Poly(entry, *xx).terms())
            family[j, i] = -family[i, j]
        witness = family.subs(point).applyfunc(sp.expand)
        denominator = sp.ilcm(*[c.q for i, j in ambient.PAIRS
                               for mon, c in sp.Poly(witness[i, j], *yy).terms()])
        witness *= denominator
        family *= denominator

        def vector(mat):
            return sp.Matrix([sp.Poly(mat[i, j], *yy).coeff_monomial(mon)
                              for (i, j), mon in ambient.BIVECTOR_BASIS])

        sparse = {(i, j): {mon: int(c) for mon, c in sp.Poly(witness[i, j], *yy).terms() if c}
                  for i, j in ambient.PAIRS if witness[i, j] != 0}
        D = sp.Matrix(ambient.linearized_jacobi_matrix(sparse).tolist())
        kernel = dm(D).nullspace().to_Matrix().T
        assert D*kernel == sp.zeros(350, kernel.cols)
        orbit = []
        for i, j in itertools.product(range(5), repeat=2):
            A = sp.zeros(5)
            A[i, j] = 1
            X = A*sp.Matrix(yy)
            orbit.append(vector(sum([X[k]*witness.diff(yy[k]) for k in range(5)], sp.zeros(5))
                                - A*witness-witness*A.T))
        M = sp.Matrix(ambient.vertical_bivector_matrix().tolist()).row_join(
            sp.Matrix.hstack(*orbit, vector(witness)))
        F = sp.Matrix.hstack(*[vector(family.diff(v).subs(point)) for v in parameters])
        assert D*M == sp.zeros(350, M.cols)
        assert D*F == sp.zeros(350, F.cols)
        dcert, mcert, fcert = nonzero_minor(D), nonzero_minor(M), nonzero_minor(M.row_join(F))
        assert dcert['rank']+kernel.cols == 150
        results[name] = {
            'full_witness_brackets': {f'{i},{j}': str(witness[i, j])
                                     for i, j in ambient.PAIRS if witness[i, j] != 0},
            'source_point': {str(k): v for k, v in point.items()},
            'integer_scale': str(denominator),
            'full_derivative': dcert,
            'orbit_vertical_scale': mcert,
            'source_orbit_vertical_scale': fcert,
            'full_kernel_dimension': kernel.cols,
            'kernel_sparse_columns': [{str(i): str(kernel[i, j]) for i in range(150) if kernel[i, j]}
                                      for j in range(kernel.cols)],
            'all_kernel_and_actual_directions_annihilated': True,
            'projective_tangent_dimension': 125-dcert['rank'],
            'projective_orbit_dimension': mcert['rank']-25,
            'projective_actual_source_rank': fcert['rank']-25,
            'smooth_component_test_passes': 150-dcert['rank'] == fcert['rank'],
        }
    report = {
        'status': 'EXACT_FINITE_COMPONENT_TESTS',
        'producer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'producer_bytes': Path(__file__).read_text(),
        'input_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
        'ambient_helper_sha256': hashlib.sha256(Path(ambient.__file__).read_bytes()).hexdigest(),
        'all_ambient_variables': 150, 'all_jacobi_rows': 350,
        'results': results,
        'seconds': round(time.monotonic()-started, 4),
    }
    output = ROOT/'reports/CLASSIFICATION_NILPOTENT_SOURCE_COMPONENTS_20261107.json'
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'report': str(output), 'seconds': report['seconds'],
                      'ranks': {k: [v['full_derivative']['rank'],
                                    v['orbit_vertical_scale']['rank'],
                                    v['source_orbit_vertical_scale']['rank']]
                                for k, v in results.items()}}))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact primitive-beta differential of an actual component source.

Use constant evaluation, all first-jet coefficients, and the four finite
perspective-gauge tangent directions. No moving Darboux formula is omitted:
the constant-evaluation kernel is exactly the fixed-value tangent space.
"""
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as sp
from sympy.polys.matrices import DomainMatrix


def dm(M):
    return DomainMatrix.from_Matrix(M).convert_to(sp.QQ)


def certificate(M):
    columns = list(dm(M).rref()[1])
    rows = list(dm(M[:, columns].T).rref()[1])
    determinant = dm(M[rows, columns]).det()
    kernel = dm(M).nullspace().to_Matrix().T
    assert M*kernel == sp.zeros(M.rows, kernel.cols)
    assert len(columns)+kernel.cols == M.cols
    return {'rank': len(columns), 'minor_rows': rows, 'minor_columns': columns,
            'minor_determinant': str(determinant),
            'kernel_sparse_columns': [{str(i): str(kernel[i, j]) for i in range(M.cols) if kernel[i, j]}
                                      for j in range(kernel.cols)]}


def main():
    started = time.monotonic()
    source_path = Path('reports/CLASSIFICATION_RANK2_NIL_31_COMPONENT_20261107.json')
    data = json.loads(source_path.read_text())
    h = sp.Symbol('h')
    xx = sp.symbols('x1:5')
    x = sp.Matrix(xx)
    nn = sp.symbols('n0:15')
    local = {str(v): v for v in (h,)+xx+nn}
    parameters = [local[v] for v in data['full_source_parameters']]
    at = {local[k]: sp.Rational(v) for k, v in data['witness'].items()}
    p = dict(zip(xx, [0, 0, 1, 1]))
    pairs = tuple(itertools.combinations(range(4), 2))
    source = sp.zeros(4)
    for key, entry in data['universal_brackets'].items():
        i, j = [int(v)-1 for v in key.split(',')]
        source[i, j] = sp.sympify(entry, locals=local).subs(h, 1)
        source[j, i] = -source[i, j]
    P = source.subs(at).applyfunc(sp.expand)
    C = P.subs(p)
    pf = sp.factor(C[0, 1]*C[2, 3]-C[0, 2]*C[1, 3]+C[0, 3]*C[1, 2])
    assert pf != 0

    def jet(M):
        return sp.Matrix([sp.expand(M[i, j]).subs(p) for i, j in pairs]
                         + [sp.diff(M[i, j], xx[k]).subs(p) for i, j in pairs for k in range(4)])

    directions = [source.diff(q).subs(at).applyfunc(sp.expand) for q in parameters]
    labels = [str(q) for q in parameters]
    for i, j in itertools.product(range(4), repeat=2):
        A = sp.zeros(4)
        A[i, j] = 1
        X = A*x
        directions.append(A*P+P*A.T-sum([X[k]*P.diff(xx[k]) for k in range(4)], sp.zeros(4)))
        labels.append(f'GL_{i}_{j}')
    for k in range(4):
        directions.append(P.diff(xx[k]))
        labels.append(f'translation_{k}')
    directions.append(P)
    labels.append('scale')
    actual = sp.Matrix.hstack(*[jet(v) for v in directions])
    eval_matrix = actual[:6, :]
    constant_kernel = dm(eval_matrix).nullspace().to_Matrix().T
    fixed_linear = actual[6:, :]*constant_kernel

    centered = x-sp.Matrix([p[v] for v in xx])
    gauge = []
    for k in range(4):
        a = sp.eye(4)[:, k]
        ca = C*a
        G = -2*centered[k]*C+centered*ca.T-ca*centered.T
        gj = jet(G)
        assert gj[:6, :] == sp.zeros(6, 1)
        gauge.append(gj)
    gauge = sp.Matrix.hstack(*gauge)
    all_jets = actual.row_join(gauge)
    ec, gc, jc, kc = (certificate(eval_matrix), certificate(gauge),
                      certificate(all_jets), certificate(gauge[6:, :].row_join(fixed_linear)))
    assert ec['rank'] == 6 and gc['rank'] == 4
    beta_rank = jc['rank']-ec['rank']-gc['rank']
    assert beta_rank == kc['rank']-4
    assert beta_rank <= 16
    report = {
        'status': 'PRIMITIVE_BETA_DOMINANCE_PROVED' if beta_rank == 16 else 'POINT_DIFFERENTIAL_ONLY',
        'producer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'producer_bytes': Path(__file__).read_text(),
        'source_report_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
        'full_source_parameters': data['full_source_parameters'],
        'source_witness': data['witness'],
        'affine_basepoint': {str(k): str(v) for k, v in p.items()},
        'constant_bivector': [[str(v) for v in row] for row in C.tolist()],
        'constant_pfaffian': str(pf),
        'actual_direction_labels': labels,
        'jet_matrix': [[str(v) for v in row] for row in all_jets.tolist()],
        'evaluation_certificate': ec,
        'perspective_certificate': gc,
        'all_jet_certificate': jc,
        'fixed_constant_linear_plus_perspective_certificate': kc,
        'primitive_beta_image_rank': beta_rank,
        'logical_scope': 'A rank16 differential proves dominance of the full16 primitive beta target. A smaller rank at this point is not a generic upper bound.',
        'seconds': round(time.monotonic()-started, 4),
    }
    output = Path('reports/CLASSIFICATION_NIL31_PRIMITIVE_BETA_DIFFERENTIAL_20261107.json')
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('status', 'primitive_beta_image_rank', 'constant_pfaffian', 'seconds')}))


if __name__ == '__main__':
    main()

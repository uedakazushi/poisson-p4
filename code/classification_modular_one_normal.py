#!/usr/bin/env python3
"""A single true-normal support test at two already certified source points."""
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
import classification_rank3_mixed_component_review as engine


def dm(M):
    return DomainMatrix.from_Matrix(M).convert_to(s.QQ)


def rank(M):
    return len(dm(M).rref()[1])


def run_case(path, name):
    data = json.loads(path.read_text())
    data = next(v for v in data['branches'] if v['name'] == name)
    x = s.symbols('h x1 x2 x3 x4')
    nn = s.symbols('n0:15')
    env = {str(v): v for v in x+nn}
    family = s.zeros(5)
    for key, expression in data['universal_brackets'].items():
        i, j = map(int, key.split(','))
        family[i, j] = s.sympify(expression, locals=env)
        family[j, i] = -family[i, j]
    parameters = [env[v] for v in data['full_source_specification']['free_parameters']]
    at = {env[k]: v for k, v in data['full_source_specification']['witness'].items()}
    P = family.subs(at)
    ps = engine.sparse(P, x)
    assert engine.jacobi(ps) == s.zeros(350, 1)
    columns = []
    for ij, mon in engine.COLUMNS:
        v = {key: dict(value) for key, value in ps.items()}
        v.setdefault(ij, {})[mon] = v.get(ij, {}).get(mon, 0)+1
        columns.append(engine.jacobi(v))
    D = s.Matrix.hstack(*columns)
    K = dm(D).nullspace().to_Matrix().T
    assert K.cols == 48
    V, O = [], []
    for a, b in itertools.product(range(5), repeat=2):
        A = s.zeros(5)
        A[a, b] = 1
        field = A*s.Matrix(x)
        E = s.Matrix(x)
        V.append(engine.vector(engine.sparse(E*field.T-field*E.T, x)))
        O.append(engine.vector(engine.sparse(A*P+P*A.T-x[b]*P.diff(x[a]), x)))
    G = s.Matrix.hstack(*V, *O, engine.vector(ps))
    M = G.row_join(s.Matrix.hstack(*[engine.vector(engine.sparse(family.diff(v).subs(at), x)) for v in parameters]))
    assert D*M == s.zeros(350, M.cols)
    assert rank(M) == 47
    pivots = dm(M.row_join(K)).rref()[1]
    new = [v-M.cols for v in pivots if v >= M.cols]
    assert len(new) == 1
    rv = K[:, new[0]]
    rv *= s.ilcm(*[v.q for v in rv])
    rv /= s.igcd(*[int(v) for v in rv])

    def tensor(vector):
        R = s.zeros(5)
        for ((i, j), mon), coefficient in zip(engine.COLUMNS, vector):
            R[i, j] += coefficient*s.prod(v**k for v, k in zip(x, mon))
        for i, j in engine.PAIRS:
            R[j, i] = -R[i, j]
        return R

    R = tensor(rv)
    assert D*rv == s.zeros(350, 1)
    Qjac = engine.jacobi(engine.sparse(R, x))
    cols = list(dm(D).rref()[1])
    J = D[:, cols]
    rows = list(dm(J.T).rref()[1])
    inverse = dm(J[rows, :]).inv()
    correction = inverse.matmul(dm(Qjac[rows, :])).to_Matrix()
    residual = Qjac-J*correction
    result = {
        'source_report': str(path),
        'source_report_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'branch': name,
        'full_D_rank': len(cols), 'full_kernel_dimension': K.cols,
        'actual_group_rank': rank(G), 'actual_source_rank': rank(M),
        'all_actual_directions_annihilated': True,
        'normal_brackets': {f'{i},{j}': str(R[i, j]) for i, j in engine.PAIRS if R[i, j] != 0},
        'selected_D_columns': cols, 'selected_D_rows': rows,
        'normal_outside_actual_tangent': True,
    }
    nz = next((i for i, v in enumerate(residual) if v != 0), None)
    if nz is not None:
        functional = s.zeros(1, 350)
        functional[0, nz] = 1
        selected = dm(J[nz, :]).matmul(inverse).to_Matrix()
        for k, row in enumerate(rows):
            functional[0, row] -= selected[0, k]
        assert functional*D == s.zeros(1, 150)
        value = (functional*Qjac)[0]
        assert value != 0
        result.update(status='NONZERO_QUADRATIC_OBSTRUCTION',
                      left_cokernel_functional={str(i): str(functional[0, i]) for i in range(350) if functional[0, i]},
                      quadratic_obstruction=str(value),
                      true_one_normal_support_dimension=0,
                      projective_component_dimension=22)
        return result
    qv = s.zeros(150, 1)
    for i, col in enumerate(cols):
        qv[col] = -correction[i]
    Q = tensor(qv)
    assert D*qv+Qjac == s.zeros(350, 1)
    tau = s.Symbol('tau')
    T = P+tau*R+tau*tau*Q
    cubic = {}
    for i, j, k in engine.TRIPLES:
        cubic[i, j, k] = s.expand(sum(T[i, ell]*s.diff(T[j, k], x[ell])
                                     + T[j, ell]*s.diff(T[k, i], x[ell])
                                     + T[k, ell]*s.diff(T[i, j], x[ell]) for ell in range(5)))
    remaining = []
    for inds in engine.QUADS:
        value = s.expand(sum((-1)**at*x[ell]*cubic[inds[:at]+inds[at+1:]] for at, ell in enumerate(inds)))
        remaining.extend(v for v in s.Poly(value, *x).coeffs() if v != 0)
    result.update(status='EXACT_POISSON_CURVE' if not remaining else 'QUADRATIC_RANGE_ONLY_UNRESOLVED',
                  second_order_correction={f'{i},{j}': str(Q[i, j]) for i, j in engine.PAIRS if Q[i, j] != 0},
                  corrected_curve_all_350_jacobi_zero=not remaining,
                  remaining_coefficients=sorted(set(map(str, remaining))))
    if not remaining and result['actual_group_rank'] == 47:
        result['projective_component_dimension'] = 23
    return result


if __name__ == '__main__':
    resource.setrlimit(resource.RLIMIT_CPU, (90, 90))
    resource.setrlimit(resource.RLIMIT_AS, (1536*1024**2, 1536*1024**2))
    started = time.monotonic()
    cases = [run_case(Path('reports/CLASSIFICATION_NIL31_SHORT_COMPONENTS_20261107.json'), 'r_one'),
             run_case(Path('reports/CLASSIFICATION_NIL22_TOP_COMPONENTS_20261107.json'), 'A_minus_one')]
    report = {'producer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'producer_bytes': Path(__file__).read_text(),
              'engine_sha256': hashlib.sha256(Path(engine.__file__).read_bytes()).hexdigest(),
              'cases': cases, 'seconds': round(time.monotonic()-started, 4)}
    out = Path('reports/CLASSIFICATION_MODULAR_ONE_NORMAL_20261107.json')
    out.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'seconds': report['seconds'], 'cases': [{'branch': v['branch'], 'status': v['status'],
                                                            'dimension': v.get('projective_component_dimension')}
                                                           for v in cases]}))

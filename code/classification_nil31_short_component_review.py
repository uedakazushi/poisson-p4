#!/usr/bin/env python3
"""Independent full-ambient check of the nil31 short-top new22 source.

Reuses the independent sparse arithmetic engine, not the producer helper.
"""
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s
import classification_rank3_mixed_component_review as audit


def main():
    started = time.monotonic()
    path = Path('reports/CLASSIFICATION_NIL31_SHORT_COMPONENTS_20261107.json')
    data = json.loads(path.read_text())
    source = next(v for v in data['branches'] if v['name'] == 'r_minus_one_third')
    x = s.symbols('h x1 x2 x3 x4')
    nn = s.symbols('n0:14')
    env = {str(v): v for v in x+nn}
    family = s.zeros(5)
    for key, f in source['universal_brackets'].items():
        i, j = map(int, key.split(','))
        family[i, j] = s.sympify(f, locals=env)
        family[j, i] = -family[i, j]
    at = {env[k]: v for k, v in source['full_source_specification']['witness'].items()}
    parameters = [env[v] for v in source['full_source_specification']['free_parameters']]
    P = family.subs(at)
    for i, j, k in audit.TRIPLES:
        assert s.expand(sum(family[i, ell]*s.diff(family[j, k], x[ell])
                            + family[j, ell]*s.diff(family[k, i], x[ell])
                            + family[k, ell]*s.diff(family[i, j], x[ell]) for ell in range(5))) == 0
    ps = audit.sparse(P, x)
    assert audit.jacobi(ps) == s.zeros(350, 1)
    columns = []
    for ij, mon in audit.COLUMNS:
        q = {key: dict(value) for key, value in ps.items()}
        q.setdefault(ij, {})[mon] = q.get(ij, {}).get(mon, 0)+1
        columns.append(audit.jacobi(q))
    D = s.Matrix.hstack(*columns)
    V, O = [], []
    for a, b in itertools.product(range(5), repeat=2):
        v = s.zeros(5, 1)
        v[a] = x[b]
        E = s.Matrix(x)
        V.append(audit.vector(audit.sparse(E*v.T-v*E.T, x)))
        A = s.zeros(5)
        A[a, b] = 1
        O.append(audit.vector(audit.sparse(A*P+P*A.T-x[b]*P.diff(x[a]), x)))
    vertical = s.Matrix.hstack(*V)
    group = vertical.row_join(s.Matrix.hstack(*O)).row_join(audit.vector(ps))
    full = group.row_join(s.Matrix.hstack(*[audit.vector(audit.sparse(family.diff(v).subs(at), x)) for v in parameters]))
    assert D*full == s.zeros(350, full.cols)
    certs = {name: audit.minor_certificate(M, rank) for name, M, rank in
             [('full_D', D, 103), ('vertical', vertical, 24), ('group_scale', group, 47), ('full_source', full, 47)]}
    Z = s.Matrix([s.expand(sum(s.diff(P[i, j], x[i]) for i in range(5))) for j in range(5)])
    expected = s.Matrix([0, 24*x[2], 24*x[3], 0, 24*x[0]])
    assert Z == expected
    A = Z.jacobian(x)
    assert [A.rank(), (A*A).rank(), (A**3).rank()] == [3, 1, 0]
    E = s.Matrix(x)
    canonical = P-(E*Z.T-Z*E.T)/5
    assert all(s.expand(canonical[0, j]+x[0]*Z[j]/5) == 0 for j in range(5))
    pf = s.expand(P[1, 2]*P[3, 4]-P[1, 3]*P[2, 4]+P[1, 4]*P[2, 3])
    cubic = s.cancel(pf/x[0])
    assert s.Poly(cubic, *x).total_degree() == 3
    leading = s.Poly(cubic, x[1]).coeff_monomial(x[1])
    remainder = cubic.subs(x[1], 0)
    assert s.Poly(cubic, x[1]).degree() == 1
    assert s.Poly(s.gcd(leading, remainder), *x).total_degree() == 0
    result = {
        'status': 'PASS',
        'producer_imported': False,
        'source_report_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'checker_bytes': Path(__file__).read_text(),
        'independent_engine_sha256': hashlib.sha256(Path(audit.__file__).read_bytes()).hexdigest(),
        'universal_raw_jacobi_all_zero': True,
        'all_actual_directions_annihilated': True,
        'certificates': certs,
        'projective_component_dimension': 22,
        'projective_orbit_dimension': 22,
        'canonical_h_cofactor': [str(-v/5) for v in Z],
        'cofactor_nilpotent_power_ranks': [3, 1, 0],
        'cofactor_Jordan_type': [3, 2],
        'actual_divisor': str(s.factor(x[0]*pf)),
        'cubic_factor_irreducible_by_primitive_linear_form': True,
        'seconds': round(time.monotonic()-started, 4),
    }
    output = Path('reports/CLASSIFICATION_NIL31_SHORT_COMPONENT_REVIEW_20261107.json')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'projective_component_dimension', 'cofactor_Jordan_type', 'seconds')}))


if __name__ == '__main__':
    main()

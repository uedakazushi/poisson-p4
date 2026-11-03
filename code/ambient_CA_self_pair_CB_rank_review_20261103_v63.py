#!/usr/bin/env python3
"""Independent differentiation of the actual nonlinear CB source on the CA pair."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_CA_SELF_PAIR_CB_SOURCE_RANK_20261103_v63'


def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def same(A, B):
    assert A.shape == B.shape and (A-B).applyfunc(s.cancel) == s.zeros(*A.shape)


def main(output, replay):
    started = time.monotonic()
    folder = ROOT/'reports'/STEM
    report = read(folder/'report.json')
    assert report['all_passed']
    for path, digest in report['input_sha256'].items():
        assert sha(ROOT/path) == digest, path
    for path, digest in report['evidence_sha256'].items():
        assert sha(folder/path) == digest, path
    subprocess.run([str(ROOT/'.venv/bin/python'),
                    str(ROOT/'code/ambient_CA_self_pair_CB_source_rank_20261103_v63.py'),
                    '--output', str(replay)], check=True, timeout=30)
    assert sha(folder/'certificate.json.gz') == sha(replay/'certificate.json.gz')
    data = read(folder/'certificate.json.gz')
    h, x, y, z, w, a, mu = s.symbols('h x y z w a mu')
    X = s.Matrix([h, x, y, z, w])
    variables = list(s.symbols('u0:19'))
    symbols = {str(q): q for q in [*X, a, mu, *variables]}
    parse = lambda q: s.sympify(q, locals=symbols)
    matrix = lambda rows: s.Matrix([[parse(q) for q in row] for row in rows])
    pairs = list(itertools.combinations(range(5), 2))
    mons = list(itertools.combinations_with_replacement(range(5), 2))

    def vec(P):
        return s.Matrix([s.Poly(P[i, j], *X).coeff_monomial(X[u]*X[v])
                         for i, j in pairs for u, v in mons])

    def canonical(P):
        d = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5))
                      for j in range(5)])
        return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)

    def raw(f, g, H):
        P = s.zeros(5)
        P[1, 2], P[3, 4] = f, g
        P[2, 3] = -s.diff(f, x)*z/2-s.diff(H, w)
        P[2, 4] = -s.diff(f, x)*w/2+s.diff(H, z)
        return P-P.T

    cubics = [s.prod((h, x, z, w)[i] for i in mon)
              for mon in itertools.combinations_with_replacement(range(4), 3)
              if any(i >= 2 for i in mon)]
    f = sum(c*q for c, q in zip(variables[:3], [h*h, h*x, x*x]))
    H = sum(c*q for c, q in zip(variables[3:], cubics))
    U = canonical(raw(f, mu*f, H))
    source = vec(U)
    same(source, s.Matrix(list(map(parse, data['universal_actual_CB150']))))
    HB = -a*x*w*z+x*w*(h+x)/(3*a)+(h+x)**2*z/(9*a*a)+z**3/3
    values = [0, 0, 3*a, *[s.Poly(HB, h, x, z, w).coeff_monomial(q)
                            for q in cubics], 0]
    point = dict(zip([*variables, mu], values))
    P = U.subs(point).applyfunc(s.cancel)
    J = source.jacobian([*variables, mu]).subs(point).applyfunc(s.cancel)
    same(J, matrix(data['all150_by20_actual_source_derivative']))
    mixed = s.Matrix.hstack(*(source.diff(mu).diff(q) for q in variables[:3]))
    same(mixed, matrix(data['all_three_mu_times_f_mixed_derivatives150_by3']))
    same(J[:, 19], vec(canonical(raw(0, 3*a*x*x, 0))))
    globalpath = ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json'
    parabolic = set(map(tuple, read(globalpath)['parabolic_GL_positions']))
    flags = [q for q in itertools.product(range(5), repeat=2) if q not in parabolic]
    assert flags == list(map(tuple, data['full_flag9_positions'])) and len(flags) == 9
    lie_columns = []
    for i, j in flags:
        E = s.zeros(5)
        E[i, j] = 1
        V = E*X
        L = E*P+P*E.T-sum((V[k]*P.diff(X[k]) for k in range(5)), s.zeros(5))
        lie_columns.append(vec(L.applyfunc(s.expand)))
    full = J.row_join(s.Matrix.hstack(*lie_columns))
    same(full, matrix(data['full150_by29_actual_CB_source_derivative']))
    F = J[:, :19]
    rows = data['fixed19_selected_rows']
    same(F[rows, :].inv(), matrix(data['fixed19_inverse']))
    quot = (full[:, 19:]-F*F[rows, :].inv()*full[rows, 19:]).applyfunc(s.cancel)
    same(quot, matrix(data['all150_by10_quotient_after_fixed19']))
    field = s.QQ.frac_field(a)
    assert len(DomainMatrix.from_Matrix(full).convert_to(field).rref()[1]) == 29
    assert len(DomainMatrix.from_Matrix(quot).convert_to(field).rref()[1]) == 10
    selected = quot[data['quotient_minor_rows'], data['quotient_minor_columns']]
    dm = DomainMatrix.from_Matrix(selected).convert_to(field)
    determinant = s.factor(field.to_sympy(dm.det()))
    assert s.cancel(determinant-parse(data['quotient_minor_determinant'])) == 0
    assert determinant == s.Rational(53747712, 15625)*a**10
    incidence = ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63_REVIEW.json'
    assert read(incidence)['all_passed']
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}_REVIEW.md',
             ROOT/'code/ambient_CA_self_pair_CB_source_rank_20261103_v63.py',
             ROOT/f'research/notes/{STEM}.md', folder/'report.json',
             folder/'certificate.json.gz', globalpath, incidence]
    result = dict(id=STEM+'_REVIEW', all_passed=True,
        fresh_author_certificate_byte_identical=True,
        independently_built_full_actual_nonlinear_CB_source_mu_times_f=True,
        all150_by20_source_derivatives_and_all_three_mixed_columns_checked=True,
        actual_ratio_column_is_canonical_raw_g_3a_xsquared=True,
        independent_physical_full9flag_derivatives_and_full150_by29_matrix_checked=True,
        exact_QQa_source_rank=29, fixed_fiber_rank=19, quotient_rank=10,
        quotient_minor=str(determinant),
        every_nonzero_a_specified_CB_lift_is_immersion=True,
        whole_CB_fiber_count_smoothness_or_completed_type_claimed=False,
        input_sha256={str(path.relative_to(ROOT)): sha(path) for path in paths},
        seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['all_passed', 'exact_QQa_source_rank',
                                           'quotient_minor', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--replay', type=Path, required=True)
    args = parser.parse_args()
    main(args.output.resolve(), args.replay.resolve())

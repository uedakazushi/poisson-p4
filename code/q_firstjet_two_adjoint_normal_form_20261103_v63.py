#!/usr/bin/env python3
"""Symbolic two-adjoint seven-coordinate hypersurface chart, no new beta points."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import sys
import time

import sympy as s

from q_firstjet_spin4_new25_point_20261103_v63 import K, mat, mm, eye
from weight1_structure import J

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_FIRSTJET_TWO_ADJOINT_NORMAL_FORM_20261103_v63'
sys.set_int_max_str_digits(0)


def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(M):
    return [[str(x) for x in row] for row in M.tolist()]


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    out.mkdir(parents=True, exist_ok=False)
    pointpath = ROOT / 'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/certificate.json.gz'
    point = read(pointpath)
    assert point['same_point_index'] == 0
    x, y, a, b, c, d, e, f, g, h, i, j = s.symbols('x y a b c d e f g h i j')
    C = s.Matrix([[a, b, g, h], [c, -a, i, j], [-j, h, d, e], [i, -g, f, -d]])
    D = s.diag(x, -x, y, -y)
    assert C.T * J + J * C == s.zeros(4)
    assert D.T * J + J * D == s.zeros(4)
    q1, q2, q3 = s.trace(C * C), s.trace(C * D), s.trace(D * D)
    C4, D4 = s.trace(C ** 4), s.trace(D ** 4)
    f1 = s.expand(q1 - 2 * q3)
    f2 = s.expand(q2 + q3)
    f3 = s.expand(C4 + 4 * D4 - 3 * q3 ** 2)
    f4 = s.expand(-C4 + 6 * s.trace(C * C * D * D) + 12 * s.trace(C * D ** 3) + 8 * D4)
    delta = x * x - y * y
    detdiff = s.expand(C.det(method='domain-ge') - delta ** 2)
    assert s.expand(f3 + 4 * detdiff - f1 * (q1 + 2 * q3) / 2) == 0
    H = s.expand((f3 + f4) / 12)
    d0 = -(a * x + x * x + y * y) / y
    M, L = g * j - h * i, (a + x) ** 2 - y ** 2
    c0, f0 = M - L, M - x ** 2 * L / y ** 2
    gauge = {b: s.Integer(1), e: s.Integer(1)}
    graph = {**gauge, d: d0, c: c0, f: f0}
    Z1, Z2, Z3 = d - d0, c - c0, f - f0
    F1 = f1.subs(gauge) / 2
    HH = H.subs(gauge)
    assert s.cancel(f2 - 2 * y * Z1) == 0
    assert s.cancel(F1 - Z2 - Z3 - (2 * d0 + Z1) * Z1) == 0
    assert s.cancel(HH - x ** 2 * Z2 - y ** 2 * Z3
                    - (2 * y ** 2 * d0 + 2 * y ** 3 + y ** 2 * Z1) * Z1) == 0
    XX = F1 - (2 * d0 + Z1) * Z1
    YY = HH - (2 * y ** 2 * d0 + 2 * y ** 3 + y ** 2 * Z1) * Z1
    assert s.cancel((YY - y ** 2 * XX) / delta - Z2) == 0
    assert s.cancel((x ** 2 * XX - YY) / delta - Z3) == 0
    graph_detdiff = s.cancel(detdiff.subs(graph))
    numerator, denominator = s.fraction(graph_detdiff)
    assert denominator == y ** 2
    variables7 = [x, y, a, g, h, i, j]
    polynomial = s.Poly(numerator, *variables7)
    assert len(polynomial.terms()) == 77 and polynomial.total_degree() == 8
    for equation, multiple in [(f1, 0), (f2, 0), (f3, -4), (f4, 4)]:
        assert s.cancel(equation.subs(graph) - multiple * numerator / y ** 2) == 0
    # Explicit telescoping determinant cofactors; no Groebner elimination.
    P = detdiff.subs(gauge)
    telescopes = []
    residual = P
    for variable, value in [(d, d0), (c, c0), (f, f0)]:
        next_residual = residual.subs(variable, value)
        quotient = s.cancel((residual - next_residual) / (variable - value))
        assert s.cancel(residual - next_residual - (variable - value) * quotient) == 0
        telescopes.append(quotient)
        residual = next_residual
    assert s.cancel(residual - numerator / y ** 2) == 0
    # Check the two torus gauges and its effective residual relative sign.
    t, u = s.symbols('t u', nonzero=True)
    T = s.diag(t, 1 / t, u, 1 / u)
    TC = T * C * T.inv()
    assert TC[0, 1] == t ** 2 * b and TC[2, 3] == u ** 2 * e
    relative = s.diag(-1, -1, 1, 1)
    assert relative * C * relative.inv() == C.subs({g: -g, h: -h, i: -i, j: -j})
    # Nonemptiness certificate on the SAME old pair, using D^2 projectors only.
    actual_C, actual_D = s.Matrix(point['actual_C4']), s.Matrix(point['actual_D4'])
    S = s.trace(actual_D * actual_D) / 2
    P0 = actual_D.det()
    discriminant = S ** 2 - 4 * P0
    assert P0 and discriminant
    K.q = Fraction(str(discriminant))
    root = K(0, 1)
    X, Y = (K(S) + root) / K(2), (K(S) - root) / K(2)
    Ik, Ck, Dk = eye(4), mat(actual_C), mat(actual_D)
    D2 = mm(Dk, Dk)
    projector = [[(D2[r][c] - Y * Ik[r][c]) / root for c in range(4)] for r in range(4)]
    assert mm(projector, projector) == projector
    conjugate_projector = [[z.conjugate() for z in row] for row in projector]
    assert [[projector[r][c] + conjugate_projector[r][c] for c in range(4)] for r in range(4)] == Ik
    def trace(A):
        return sum((A[k][k] for k in range(4)), K(0))
    trDC = trace(mm(mm(Dk, Ck), projector))
    bc = trace(mm(mm(mm(projector, Ck), projector), Ck)) / K(2) - trDC * trDC / (K(4) * X)
    ef = bc.conjugate()
    norm = bc * ef
    assert norm and norm.b == 0
    assert norm.a.denominator % 32003
    cert = dict(generic_C=rows(C), generic_D=rows(D),
                four_original_constraint_polynomials=list(map(str, [f1, f2, f3, f4])),
                determinant_replacement_identity='f3+4(detC-(x^2-y^2)^2)=f1*(q1+2q3)/2',
                linear_constraint_H=str(H), gauge_b_e=[1, 1],
                graph=dict(d=str(d0), c=str(c0), f=str(f0), M=str(M), L=str(L)),
                eliminated_graph_difference_variables=list(map(str, [Z1, Z2, Z3])),
                determinant_telescoping_cofactors=list(map(str, telescopes)),
                final_coordinates=list(map(str, variables7)), final_polynomial=str(numerator),
                final_polynomial_sparse=[[list(m), str(c)] for m, c in polynomial.terms()],
                final_graph_detdiff_denominator=str(denominator),
                required_known_unit_conditions='x*y*(x^2-y^2)*b*e !=0 before gauge; x*y*(x^2-y^2)!=0 after gauge',
                normal_form_ideal_equals_graph_plus_one_equation_after_localization=True,
                same_point_index=0, actual_D_squared_trace_half=str(S), actual_detD=str(P0),
                actual_D_squared_discriminant=str(discriminant),
                quadratic_etale_algebra_encoding='[constant,coefficient of s], s^2=D2_discriminant',
                D_squared_X_Y=[X.save(), Y.save()],
                D_squared_projector=[[[str(z.a), str(z.b)] for z in row] for row in projector],
                actual_block_bc=bc.save(), actual_block_ef=ef.save(),
                actual_block_bc_times_ef=norm.save(),
                actual_block_product_mod32003=norm.a.numerator * pow(norm.a.denominator, -1, 32003) % 32003)
    with gzip.open(out / 'certificate.json.gz', 'wt') as f:
        json.dump(cert, f, indent=2)
        f.write('\n')
    result = dict(id=STEM, all_passed=True, exact_field='QQ',
                  symbolic_four_constraint_graph_and_reverse_identities=True,
                  no_Groebner_basis_used=True, final_coordinates=7, final_equations=1,
                  final_polynomial_terms=77, final_polynomial_total_degree=8,
                  actual_image_component_meets_regular_D_and_gauge_open=True,
                  same_previously_used_point_only=True,
                  effective_torus_gauge_residual='relative sign on all four cross-block coordinates',
                  finite_normalized_cross_section_over_generic_orbits=True,
                  quotient_rationality_or_global_image_exhaustion_claimed=False,
                  NEW25_equation_claimed=False, no_new_beta_points_or_degree22plus_probe=True,
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__),
                      ROOT / 'code/weight1_structure.py', ROOT / 'code/q_firstjet_spin4_new25_point_20261103_v63.py', pointpath]},
                  seconds=time.monotonic() - started)
    (out / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'final_polynomial_terms',
          'actual_image_component_meets_regular_D_and_gauge_open', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output.resolve())

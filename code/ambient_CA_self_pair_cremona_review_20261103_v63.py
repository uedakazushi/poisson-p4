#!/usr/bin/env python3
"""Independent all-nonzero-parameter CA Cremona, center, leaf and log-source audit."""
import argparse
import ast
import gzip
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_CA_SELF_PAIR_CREMONA_20261103_v63_REVIEW'
PAIRS = list(it.combinations(range(5), 2))
MONS = list(it.combinations_with_replacement(range(5), 2))


def load(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(path):
    report = load(path)
    assert report['all_passed']
    for name, digest in report.get('input_sha256', {}).items():
        assert sha(ROOT / name) == digest, name
    for name, digest in report.get('evidence_sha256', {}).items():
        assert sha(path.parent / name) == digest, name
    return report


def parity(indices):
    return (-1) ** sum(indices[i] > indices[j] for i in range(len(indices)) for j in range(i + 1, len(indices)))


def main(output):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(120)
    output.mkdir(parents=True, exist_ok=False)
    producer = ROOT / 'reports/AMBIENT_CA_SELF_PAIR_CREMONA_20261103_v63'
    check(producer / 'report.json')
    evidence = load(producer / 'certificate.json')
    parent = ROOT / 'reports/AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63'
    check(parent / 'report.json')
    actual = load(parent / 'certificate.json.gz')
    physicalpath = ROOT / 'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    physical = load(physicalpath)
    h, x, y, z, w, a = s.symbols('h x y z w a')
    coords = s.Matrix([h, x, y, z, w])
    variables = {str(v): v for v in (*coords, a)}

    def parse(value, more=None):
        return s.sympify(value, locals=variables | (more or {}))

    def canonical(P, X):
        div = s.Matrix([sum(s.diff(P[i, j], X[i]) for i in range(5)) for j in range(5)])
        return (P - (X * div.T - div * X.T) / 5).applyfunc(s.expand)

    def raw_CA(f, g, Hamiltonian):
        P = s.zeros(5)
        P[0, 2], P[3, 4] = f, g * x * x
        P[2, 3], P[2, 4] = -s.diff(Hamiltonian, w), s.diff(Hamiltonian, z)
        return P - P.T

    def flatten(P, X=coords):
        return s.Matrix([s.Poly(s.expand(P[i, j]), *X).coeff_monomial(X[u] * X[v]) for i, j in PAIRS for u, v in MONS])

    # Bind the polynomial paired source to an independently reconstructed raw CA.
    Hamiltonian = -x * h * w + h * h * z + z ** 3 / 3 - a * x * z * w
    P = canonical(raw_CA(-x * x - 3 * a * h * x, 0, Hamiltonian), coords)
    stored = s.Matrix([sum(s.Rational(c) * a ** int(d) for d, c in q) for q in actual['first_source_physical150']])
    assert flatten(P) == stored
    cubic = [s.prod((x, h, z, w)[i] for i in inds) for inds in it.combinations_with_replacement(range(4), 3) if any(i >= 2 for i in inds)]
    C = s.Matrix.hstack(*[flatten(canonical(raw_CA((x * x, h * x, h * h)[i] if i < 3 else 0, int(i == 3), cubic[i - 4] if i >= 4 else s.Integer(0)), coords)) for i in range(20)])
    assert C == s.Matrix([[s.Rational(c) for c in row] for row in physical['fixed_fiber_linear_basis150_by20']])
    fixed = s.Matrix([sum(s.Rational(c) * a ** int(d) for d, c in q) for q in actual['first_source_fixed20']])
    assert C * fixed == stored and canonical(P, coords) == P
    ey = s.Matrix([0, 0, 1, 0, 0])
    V = s.Matrix([x * x + 12 * a * h * x / 5, -3 * a * x * x / 5, 0, h * x + 2 * a * x * z / 5, h * h + z * z - 8 * a * x * w / 5])
    assert P == ey * V.T - V * ey.T

    # Derive the projective chart bracket by differentiating q_i / x.
    quotients = s.Matrix([h / x, y / x, z / x, w / x])
    chart = (quotients.jacobian(coords) * P * quotients.jacobian(coords).T).subs(x, 1).applyfunc(s.expand)
    D = [1 + 3 * a * h, h + a * z, h * h + z * z - a * w]
    ea = s.Matrix([0, 1, 0, 0])
    field = s.Matrix([D[0], 0, D[1], D[2]])
    assert chart == ea * field.T - field * ea.T
    derive = lambda f: s.expand(sum(d * s.diff(f, v) for d, v in zip(D, (h, z, w))))

    # Solve (D+a)Q = h^2+z^2 in original coordinates; no eigenmonomial formula.
    basis = [s.Integer(1), h, z, h * h, h * z, z * z]
    monoms = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
    operator = s.Matrix([[s.Poly(derive(b) + a * b, h, z).coeff_monomial(h ** i * z ** j) for b in basis] for i, j in monoms])
    assert s.factor(operator.det()) == 840 * a ** 6
    rhs = s.Matrix([0, 0, 0, 1, 0, 1])
    q = operator.inv() * rhs
    Q_original = s.expand(sum(c * b for c, b in zip(q, basis)))
    assert s.expand(derive(Q_original) + a * Q_original - h * h - z * z) == 0
    H, Z, W = s.symbols('H Z W')
    Hv = h + 1 / (3 * a)
    Zv = z - h / (2 * a) - 1 / (2 * a * a)
    h_inverse, z_inverse = H - 1 / (3 * a), Z + H / (2 * a) + 1 / (3 * a * a)
    Q = s.expand(Q_original.subs({h: h_inverse, z: z_inverse}, simultaneous=True))
    assert s.cancel(Q - parse(evidence['Q'], {'H': H, 'Z': Z})) == 0
    Wv = w - Q_original
    assert [s.cancel(derive(Hv) - 3 * a * Hv), s.cancel(derive(Zv) - a * Zv), s.cancel(derive(Wv) + a * Wv)] == [0, 0, 0]
    affine_forward = s.Matrix([Hv, y, Zv, Wv])
    affine_variables = s.Matrix([h, y, z, w])
    assert affine_forward.jacobian(affine_variables).det() == 1
    affine_inverse = s.Matrix([h_inverse, y, z_inverse, W + Q])
    assert (affine_forward.subs(dict(zip(affine_variables, affine_inverse)), simultaneous=True) - s.Matrix([H, y, Z, W])).applyfunc(s.cancel) == s.zeros(4, 1)

    # Homogenize our independently solved affine maps.
    Qh = lambda U, T, S: s.expand(S * S * Q.subs({H: U / S, Z: T / S}, simultaneous=True))
    A, B = h + x / (3 * a), z - h / (2 * a) - x / (2 * a * a)
    Cpoly = s.expand(x * w - Qh(A, B, x))
    forward = s.Matrix([x * A, x * x, x * y, x * B, Cpoly])
    HH, XX, YY, ZZ, WW = target = s.Matrix(s.symbols('HH XX YY ZZ WW'))
    inverse = s.Matrix([XX * (HH - XX / (3 * a)), XX ** 2, XX * YY, XX * (ZZ + HH / (2 * a) + XX / (3 * a * a)), XX * WW + Qh(HH, ZZ, XX)])
    assert (inverse.subs(dict(zip(target, forward)), simultaneous=True) - x ** 3 * coords).applyfunc(s.expand) == s.zeros(5, 1)
    assert (forward.subs(dict(zip(coords, inverse)), simultaneous=True) - XX ** 3 * target).applyfunc(s.expand) == s.zeros(5, 1)
    assert s.Poly(Cpoly, *coords).coeff_monomial(z ** 2) == -1 / (3 * a)
    assert s.Poly(inverse[4], *target).coeff_monomial(ZZ ** 2) == 1 / (3 * a)
    # Inverse homogeneous-degree-one map, checking the opposite tensor direction.
    psi = inverse / XX
    diagonal = s.Matrix([3 * a * XX * HH, 0, 0, a * XX * ZZ, -a * XX * WW])
    target_P = canonical(ey * diagonal.T - diagonal * ey.T, target)
    back = psi.jacobian(target) * target_P * psi.jacobian(target).T
    assert (back - P.subs(dict(zip(coords, psi)), simultaneous=True)).applyfunc(s.cancel) == s.zeros(5)
    assert s.cancel(psi.jacobian(target).det() - 1) == 0
    assert (psi.jacobian(target) * target - psi).applyfunc(s.cancel) == s.zeros(5, 1)

    # Two rational Casimirs and an actual inverse function-field chart.
    f, g = Hv / Zv ** 3, Zv * Wv
    assert s.cancel(derive(f)) == s.cancel(derive(g)) == 0
    c1, c2, t, Y = s.symbols('c1 c2 t Y')
    leaf = s.Matrix([c1 * t ** 3 - 1 / (3 * a), Y, t + c1 * t ** 3 / (2 * a) + 1 / (3 * a * a), c2 / t + Q.subs({H: c1 * t ** 3, Z: t}, simultaneous=True)])
    inverse_leaf = dict(zip(affine_variables, leaf))
    assert s.cancel(f.subs(inverse_leaf, simultaneous=True) - c1) == s.cancel(g.subs(inverse_leaf, simultaneous=True) - c2) == 0
    assert s.cancel(Zv.subs(inverse_leaf, simultaneous=True) - t) == 0
    assert s.cancel(derive(Zv) - a * Zv) == 0
    curve = [s.expand(t * leaf[0]), t, s.expand(t * leaf[2]), s.expand(t * leaf[3])]
    assert [s.degree(v, t) for v in curve] == [4, 1, 4, 7]
    highest7 = s.factor(s.Poly(curve[3], t).coeff_monomial(t ** 7))
    highest5 = s.factor(s.Poly(curve[3], t).coeff_monomial(t ** 5))
    assert highest7 == c1 ** 2 * (4 * a * a + 1) / (28 * a ** 3)
    assert highest5 == c1 / (5 * a * a)
    assert s.Poly(curve[3], t).coeff_monomial(t ** 6) == 0
    assert curve[3].subs(t, 0) == c2
    # Exact quadratic-extension specialization, not complex samples.
    special_degrees = []
    for polynomial in curve:
        degree = -1
        for (power,), coefficient in s.Poly(polynomial, t).terms():
            numerator, denominator = s.cancel(coefficient).as_numer_denom()
            assert s.gcd(s.Poly(denominator, a), s.Poly(4 * a * a + 1, a)).degree() == 0
            if s.rem(numerator, 4 * a * a + 1, a) != 0:
                degree = max(degree, power)
        special_degrees.append(degree)
    assert special_degrees == [4, 1, 4, 5]
    # t is recovered from B/x on the quotient curve, proving generic degree1.
    assert s.cancel((curve[2] - curve[0] / (2 * a) - curve[1] / (2 * a * a)) / curve[1] - t) == 0

    # Reconstruct exact legacy log coefficients from its source AST.
    legacylog = ROOT / 'legacy/v14/code/rank2_lower_bounds/verify_rank2_components.py'
    tree = ast.parse(legacylog.read_text())
    fun = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'log_bivector')
    assignment = next(n for n in fun.body if isinstance(n, ast.Assign) and any(isinstance(tgt, ast.Name) and tgt.id == 'coefficients' for tgt in n.targets))
    code = compile(ast.Expression(assignment.value), str(legacylog), 'eval')
    legacy_coefficients = eval(code, {'__builtins__': {}}, {'a': s.Rational(1, 4), 'b': -s.Rational(3, 4), 'c': s.Rational(7, 4)})
    r, rr = [1, -3, 2, 0], [0, 1, -3, 1]
    assert sum(i * j for i, j in zip(r, (1, 1, 1, 2))) == sum(i * j for i, j in zip(rr, (1, 1, 1, 2))) == 0
    assert legacy_coefficients == {(i, j): r[i] * rr[j] - r[j] * rr[i] for i, j in it.combinations(range(4), 2)}
    factors = [A, B, x, Cpoly]
    # Direct dlog f wedge dlog g, rather than the author's residue-pair sum.
    dlogf = [sum(r[k] * s.diff(factors[k], v) / factors[k] for k in range(4)) for v in coords]
    dlogg = [sum(rr[k] * s.diff(factors[k], v) / factors[k] for k in range(4)) for v in coords]
    product = s.prod(factors)
    omega = {(i, j): s.cancel(product * (dlogf[i] * dlogg[j] - dlogf[j] * dlogg[i])) for i, j in PAIRS}
    for (i, j), coefficient in omega.items():
        complement = [v for v in range(5) if v not in (i, j)]
        u, v, k = complement
        dual = parity(complement + [i, j]) * (coords[u] * P[v, k] - coords[v] * P[u, k] + coords[k] * P[u, v])
        assert s.cancel(coefficient + dual / a) == 0
        assert s.cancel(coefficient - parse(evidence['all10_log_form_coefficients'][str((i, j))])) == 0
    three_form = {(i, j, k): s.expand(s.diff(omega[j, k], coords[i]) - s.diff(omega[i, k], coords[j]) + s.diff(omega[i, j], coords[k])) for i, j, k in it.combinations(range(5), 3)}
    log_tensor = s.zeros(5)
    for i, j in PAIRS:
        complement = tuple(v for v in range(5) if v not in (i, j))
        log_tensor[i, j] = parity([i, j, *complement]) * three_form[complement]
        log_tensor[j, i] = -log_tensor[i, j]
    assert (log_tensor + 5 * P / a).applyfunc(s.cancel) == s.zeros(5)
    assert flatten(log_tensor) == s.Matrix([parse(v) for v in evidence['all150_log_tensor_coefficients']])

    # Exact LPB family: move constant vector e_y to the legacy projection center.
    U = s.Matrix(s.symbols('u0:5'))
    permutation = [0, 1, 4, 2, 3]
    E = s.zeros(5)
    for i, j in enumerate(permutation):
        E[j, i] = 1
    renamed = (E * P * E.T).subs(dict(zip(coords, E.T * U)), simultaneous=True)
    vnew = E * V.subs(dict(zip(coords, E.T * U)), simultaneous=True)
    e4 = s.Matrix([0, 0, 0, 0, 1])
    assert renamed == e4 * vnew.T - vnew * e4.T
    assert all(s.diff(v, U[4]) == 0 for v in vnew)
    assert canonical(renamed, U) == renamed
    oldpath = ROOT / 'legacy/v14/data/six_rank2_components.json'
    old = load(oldpath)
    cases = old if isinstance(old, list) else old.get('cases', old.get('components'))
    for name, dimension in [('LPB', 39), ('Log1112', 28)]:
        case = next(v for v in cases if v['name'] == name)
        assert case['projective_component_dimension'] == dimension
        assert case['ambient_affine_tangent_dimension'] == case['family_affine_differential_rank'] == dimension + 1
    # All affine and homogeneous polynomial maps are over QQ[a,a^-1].
    for expression in [*affine_forward, *affine_inverse, *forward, *inverse, Q]:
        denominator = s.denom(s.cancel(expression))
        assert denominator.free_symbols <= {a}
        assert len(s.Poly(denominator, a).terms()) == 1

    independently_derived = dict(operator_on_quadratic_polynomials_determinant=str(s.factor(operator.det())),
        Q_original_h_z=str(Q_original), Q_H_Z=str(Q),
        all150_raw_CA_and_canonical_parent_binding=True,
        full_inverse_tensor_pullback=True, quadratic_compositions_x3_X3=True,
        maps_no_common_coordinate_factor_for_all_a_nonzero=True,
        exact_rational_center='C(H/Z^3,ZW)', generic_leaf_chart='C times Gm',
        generic_curve_degrees=[4, 1, 4, 7], exact_special_curve_degrees=special_degrees,
        highest_degree7_coefficient=str(highest7), highest_degree5_coefficient=str(highest5),
        leaf_degree7_or5_basepoint_free_and_birational_proof_accepted=True,
        legacy_log_parameters=['1/4', '-3/4', '7/4'],
        direct_dlog_forms_all10_and_tensor_all150=True,
        exact_log_multiplier='-5/a', exact_LPB_source=True,
        closure_adds_a0_to_both_known_components=True,
        N_or_CB_membership_new_component_conductor_claim=False)
    (output / 'independent-evidence.json').write_text(json.dumps(independently_derived, indent=2) + '\n')
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md',
        ROOT / 'code/ambient_CA_self_pair_cremona_20261103_v63.py',
        ROOT / 'research/notes/AMBIENT_CA_SELF_PAIR_CREMONA_20261103_v63.md',
        producer / 'report.json', producer / 'certificate.json',
        parent / 'report.json', parent / 'certificate.json.gz', physicalpath,
        legacylog, ROOT / 'legacy/v14/code/rank2_lower_bounds/lpb_family.py', oldpath]
    report = dict(id=STEM, all_passed=True, producer_rerun=False,
        all_nonzero_a_exact_Poisson_Cremona_2_2_accepted=True,
        exact_rational_center_and_generic_leaf_degree7_drop5_accepted=True,
        actual_LPB39_Log1112_28_membership_including_a0_by_closure_accepted=True,
        CB_N_or_full_selfpair_conductor_claim=False,
        independent_method='Raw CA20; original-coordinate 6x6 operator inverse; inverse tensor pullback; exact quadratic-extension degree; direct dlog and legacy source AST.',
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-evidence.json': sha(output / 'independent-evidence.json')},
        seconds=time.monotonic() - start)
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

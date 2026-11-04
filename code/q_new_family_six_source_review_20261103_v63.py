#!/usr/bin/env python3
"""Independent actual source6/primitive circuit and exact 30x26 kernel audit."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import resource
import signal
import time

import sympy as s
from sympy.polys.rings import ring

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_NEW_FAMILY_SIX_SOURCE_20261103_v63_REVIEW'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for v in iter(lambda: f.read(1 << 20), b''):
            h.update(v)
    return h.hexdigest()


def pa(a, b, c=1):
    r = dict(a)
    for m, v in b.items():
        r[m] = r.get(m, 0)+c*v
        if not r[m]:
            del r[m]
    return r


def pm(a, b):
    r = defaultdict(int)
    for m, v in a.items():
        for n, w in b.items():
            r[m+n] += v*w
    return {m: v for m, v in r.items() if v}


def ps(values):
    out = {}
    for a in values:
        out = pa(out, a)
    return out


def main(out):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    folder = ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63'
    kernel_folder = ROOT/'reports/Q_NEW_FAMILY_SOURCE_KERNEL_20261103_v63'
    report, cert, dag = [read(folder/q) for q in ['report.json', 'certificate.json', 'source-six-circuit.json.gz']]
    kr, mc, kc = [read(kernel_folder/q) for q in ['report.json', 'matrix-certificate.json.gz', 'kernel-certificate.json.gz']]
    family_path = ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz'
    geo_path = ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json'
    spin_path = ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    inv_path = ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    family, geometry, spin, inv = map(read, [family_path, geo_path, spin_path, inv_path])
    gates = [ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63_REVIEW/report.json',
             ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63_REVIEW/report.json']
    assert report['all_passed'] and kr['all_passed'] and all(read(p)['all_passed'] for p in gates)
    u, v, c, X, W, B = s.symbols('u v c X W B')
    xx = s.Matrix(s.symbols('x y z w'))
    q = s.Matrix(s.symbols('q0:4'))
    variables = [u, v, c, X, W, B]
    pairs = list(itertools.combinations(range(4), 2))
    C0 = s.Matrix(family['full_affine_P']).subs(dict(zip(s.symbols('n0:3'), [u, 1, v])))
    Z = -s.Rational(3, 4)/u*(B+W/1536-W**2+X**2/64+s.Rational(2, 9))
    point = dict(zip(xx, [X, 0, Z, W]))
    value = C0.subs(point, simultaneous=True).applyfunc(s.cancel)
    Lcoeff = [C0.diff(z).subs(point, simultaneous=True).applyfunc(s.cancel) for z in xx]
    linear = sum((M*q[i] for i, M in enumerate(Lcoeff)), s.zeros(4))
    assert value == s.Matrix(cert['raw_value4'])
    assert (linear-s.Matrix(cert['raw_linear4'])).applyfunc(s.expand) == s.zeros(4)
    assert s.cancel(value[0, 1]*value[2, 3]-value[0, 2]*value[1, 3]+value[0, 3]*value[1, 2]-8*B/9) == 0
    R = s.Matrix(cert['rational_Darboux_frame4'])
    J = s.diag(s.Matrix([[0, 1], [-1, 0]]), s.Matrix([[0, 1], [-1, 0]]))
    assert (R*(c*value)*R.T-J).applyfunc(s.cancel) == s.zeros(4)
    assert (R-s.diag(1, 1/c, 1, 1/c)*R.subs(c, 1)).applyfunc(s.cancel) == s.zeros(4)
    Ri = R.inv()
    normcoeff = [(c*sum((R*Lcoeff[k]*R.T*Ri[k, j] for k in range(4)), s.zeros(4))).applyfunc(s.cancel) for j in range(4)]
    normalized = sum((M*q[j] for j, M in enumerate(normcoeff)), s.zeros(4))
    assert (normalized-s.Matrix(cert['normalized_linear4'])).applyfunc(s.cancel) == s.zeros(4)
    alpha = s.Matrix([(M[0, 1]+M[2, 3])/5 for M in normcoeff])
    assert (alpha-s.Matrix(cert['primitive_perspective_alpha'])).applyfunc(s.cancel) == s.zeros(4, 1)
    va = J*alpha
    prim = []
    for k, M in enumerate(normcoeff):
        ek = s.eye(4)[:, k]
        prim.append((-J*(M-2*alpha[k]*J+ek*va.T-va*ek.T)*J).applyfunc(s.cancel))
    assert all(s.cancel(M[0, 1]+M[2, 3]) == 0 for M in prim)
    assert all(s.cancel(prim[i][j, k]-prim[j][i, k]+prim[k][i, j]) == 0 for i, j, k in itertools.combinations(range(4), 3))
    beta24 = s.Matrix([prim[k][i, j] for i, j in pairs for k in range(4)])
    basis = s.Matrix(inv['primitive24_by16_basis'])
    ids = list(basis.T.rref()[1])
    beta16 = (basis[ids, :].inv()*beta24[ids, :]).applyfunc(s.cancel)
    assert (basis*beta16-beta24).applyfunc(s.cancel) == s.zeros(24, 1)
    assert (beta16-s.Matrix(cert['primitive_beta16'])).applyfunc(s.cancel) == s.zeros(16, 1)
    assert (beta24-s.Matrix(cert['primitive_beta24'])).applyfunc(s.cancel) == s.zeros(24, 1)
    weights = [(1-sum(w))//2 for w in inv['Cartan_weights']]
    assert weights == cert['exact_c_exponents_beta16']
    assert all(s.cancel(b-c**w*b.subs(c, 1)) == 0 for b, w in zip(beta16, weights))
    # The two finite gauges are bound to accepted whole-tensor formulas.
    Y, Z0, oldn, oldc = s.symbols('Y Z0 oldn oldc')
    gt = s.Matrix(geometry['modular_Ga_finite_action']).subs(dict(zip(s.symbols('n0:3'), [u, 1, v])))
    gt = gt.subs(s.Symbol('t'), -3*Y/(4*u))
    moved = gt*s.Matrix([1, X, Y, Z0, W])
    assert s.cancel(moved[2]) == 0
    assert s.cancel(moved[3]-(Z0-3*X*Y/4-3*v*Y**2/8-Y/(64*u))) == 0
    eq = geometry['parameter_equivalences'][0]
    assert eq['tensor_scalar'] == 'a**(-2)' and eq['parameter_image'] == ['n0', 'a*n1', 'n2/a']
    assert s.cancel(oldc*(1/oldn)**2-oldc/oldn**2) == 0
    # Independently evaluate only the ancestors of the sixteen beta roots in QQ(source6).
    ring6, *gens = ring('u,v,c,X,W,B', s.QQ)
    field6 = ring6.to_field()
    nodes, outputs = dag['nodes'], dag['outputs']
    assert len(nodes) == cert['circuit_nodes'] == 1025
    cache = {}
    def evaluate(i):
        if i in cache:
            return cache[i]
        op, arg = nodes[i]
        if op == 'q':
            z = field6(s.Rational(arg))
        elif op == 'v':
            z = field6(gens[arg])
        elif op == 'inv':
            z = 1/evaluate(arg)
        elif op == '+':
            z = sum((evaluate(j) for j in arg), field6.zero)
        else:
            assert op == '*'
            z = field6.one
            for j in arg:
                z *= evaluate(j)
        cache[i] = z
        return z
    assert all(evaluate(i) == field6.from_expr(z) for i, z in zip(outputs['beta16'], beta16))
    lookup = {}
    for i, (op, arg) in enumerate(nodes):
        if op in ('+', '*'):
            assert all(j < i for j in arg)
            key = op, tuple(arg)
        elif op == 'inv':
            assert arg < i
            key = op, arg
        else:
            key = op, arg
        assert key not in lookup
        lookup[key] = i
    zero, one = lookup['q', '0'], lookup['q', '1']
    def const(x):
        return lookup['q', str(s.Rational(x))]
    def plus(*values):
        z = sorted(i for i in values if i != zero)
        return zero if not z else z[0] if len(z) == 1 else lookup['+', tuple(z)]
    def times(*values):
        if zero in values:
            return zero
        z = sorted(i for i in values if i != one)
        return one if not z else z[0] if len(z) == 1 else lookup['*', tuple(z)]
    def mm(A, B):
        return [[plus(*(times(A[i][k], B[k][j]) for k in range(4))) for j in range(4)] for i in range(4)]
    tr = lambda A: plus(*(A[i][i] for i in range(4)))
    monnodes = [times(outputs['beta16'][i], outputs['beta16'][j]) for i, j in spin['quadratic_monomials']]
    def cov(key):
        return [[plus(*(times(const(a), monnodes[k]) for k, a in enumerate(spin[key][4*i+j]) if F(a))) for j in range(4)] for i in range(4)]
    C, D = map(cov, ['spin_adjoint_C_coefficients16x136', 'vector_adjoint_D_coefficients16x136'])
    assert C == outputs['C4'] and D == outputs['D4']
    CC, CD, DD = mm(C, C), mm(C, D), mm(D, D)
    S = times(const(s.Rational(1, 2)), tr(DD))
    PD = plus(times(const(s.Rational(1, 8)), tr(DD), tr(DD)), times(const(-s.Rational(1, 4)), tr(mm(DD, DD))))
    disc = plus(times(S, S), times(const(-4), PD))
    N = plus(tr(mm(mm(C, DD), mm(C, DD))), times(const(-1), S, tr(mm(CC, DD))), times(PD, tr(CC)))
    assert disc == outputs['spectral_discriminant8'] and N == outputs['N12']
    assert times(N, const(s.Rational(1, 2)), lookup['inv', disc]) == outputs['m']
    moment = [[plus(times(const(-1), C[i][j]), times(const(2), D[i][j])) for j in range(4)] for i in range(4)]
    assert tr(mm(moment, moment)) == outputs['F4']
    words = [tr(mm(A, B)) for A, B in [(CC, CC), (CC, CD), (CC, DD), (CD, CD), (CD, DD)]]
    assert words == outputs['trace8_basis']
    change = s.Matrix([[0, 0, -s.Rational(1, 2), 0, -1], [-s.Rational(1, 2), -1, 1, 0, 2],
        [s.Rational(1, 2), 0, -2, -s.Rational(1, 2), -2], [-s.Rational(1, 4), 0, s.Rational(1, 2), 0, 0],
        [s.Rational(1, 32), 0, s.Rational(1, 16), 0, s.Rational(1, 8)]])
    assert [plus(*(times(const(a), words[j]) for j, a in enumerate(row) if a)) for row in change.tolist()] == outputs['determinant_pencil5']
    # Prove that pencil-change matrix as full polynomials on the original W16,
    # by a separate denominator-cleared integer trace-word engine.
    entries = [spin[k] for k in ['spin_adjoint_C_coefficients16x136', 'vector_adjoint_D_coefficients16x136']]
    den = math.lcm(*(F(v).denominator for A in entries for row in A for v in row))
    powers = [9**i for i in range(16)]
    monkeys = [powers[i]+powers[j] for i, j in spin['quadratic_monomials']]
    def polys(A):
        return [[{monkeys[k]: int(F(a)*den) for k, a in enumerate(A[4*i+j]) if F(a)} for j in range(4)] for i in range(4)]
    def pmm(A, B):
        return [[ps(pm(A[i][k], B[k][j]) for k in range(4)) for j in range(4)] for i in range(4)]
    CP, DP = map(polys, entries)
    C2, CDP, D2 = pmm(CP, CP), pmm(CP, DP), pmm(DP, DP)
    tw = [ps(pm(A[i][j], B[j][i]) for i in range(4) for j in range(4)) for A, B in [(C2, C2), (C2, CDP), (C2, D2), (CDP, CDP), (CDP, D2)]]
    expected_pencil = [{int(m): F(v) for m, v in pp} for pp in spin['canonical_adjoint_pencil']['determinant_coefficients']]
    for row, pp in zip(change.tolist(), expected_pencil):
        ans = {}
        for a, pol in zip(row, tw):
            if a:
                ans = pa(ans, pol, F(a)/den**4)
        assert ans == pp
    # Original30x26 coefficient-action matrix, not author symbolic Lie formula.
    flat = lambda A, LC: [A[i, j] for i, j in pairs]+[LC[k][i, j] for i, j in pairs for k in range(4)]
    parameters = [u, v, X, W, B]
    columns = [s.Matrix(flat(value.diff(z), [M.diff(z) for M in Lcoeff])) for z in parameters]
    columns.append(s.Matrix(flat(value, Lcoeff)))
    for i, j in itertools.product(range(4), repeat=2):
        V = s.zeros(4)
        V[i, j] = 1
        lc = [V*M+M*V.T-sum((Lcoeff[a]*V[a, k] for a in range(4)), s.zeros(4)) for k, M in enumerate(Lcoeff)]
        columns.append(s.Matrix(flat(V*value+value*V.T, lc)))
    for j in range(4):
        alpha = s.eye(4)[:, j]
        va = value*alpha
        lc = [-2*alpha[k]*value+s.eye(4)[:, k]*va.T-va*s.eye(4)[k, :] for k in range(4)]
        columns.append(s.Matrix(flat(s.zeros(4), lc)))
    M = s.Matrix.hstack(*columns).applyfunc(s.expand)
    assert M == s.Matrix(mc['full30x26_polynomial_matrix'])
    PR, *pg = ring('u,v,X,W,B', s.QQ)
    matrix = [[PR.from_expr(z) for z in row] for row in M.tolist()]
    rowids, colids = list(range(30)), list(range(26))
    for step in mc['all_twenty_QQ_unit_pivots']:
        i, j = rowids.index(step['original_row']), colids.index(step['original_column'])
        assert colids == step['current_columns']
        row, pivot = matrix[i], PR.from_expr(s.Rational(step['pivot']))
        assert row[j] == pivot and [str(z) for z in row] == step['row_polynomials']
        matrix = [[f-r[j]/pivot*row[k] for k, f in enumerate(r) if k != j] for ii, r in enumerate(matrix) if ii != i]
        del rowids[i]
        del colids[j]
    assert rowids == mc['remaining_rows'] and colids == mc['remaining_columns']
    assert matrix == [[PR.from_expr(s.sympify(z)) for z in row] for row in mc['residual10x6_polynomial_matrix']]
    kernel = [PR.from_expr(s.sympify(z)) for z in kc['full26_polynomial_kernel']]
    original = [[PR.from_expr(z) for z in row] for row in M.tolist()]
    assert all(sum((a*b for a, b in zip(row, kernel)), PR.zero) == 0 for row in original)
    assert kernel[0] == 0 and any(kernel[1:6])
    # Rows1,2 each have one essential pivot, followed by one explicit3x3 block.
    resid = s.Matrix([[z.as_expr() for z in row] for row in matrix])
    minor = s.factor(resid[[1, 2, 3, 4, 7], [0, 1, 2, 3, 5]].det(method='domain-ge'))
    assert minor != 0
    assert any(kernel[j] for j in colids)
    evidence = dict(independent_raw_matrix30_by26=[[str(z) for z in row] for row in M.tolist()],
        exact_residual5minor_rows=[1, 2, 3, 4, 7], exact_residual5minor_columns=[0, 1, 2, 3, 5],
        exact_residual5minor=str(minor), actual_beta16=list(map(str, beta16)),
        all26_verified_kernel=[str(z) for z in kernel], universal_pencil_change_verified=True)
    with gzip.open(out/'independent-certificate.json.gz', 'wt') as f:
        json.dump(evidence, f)
        f.write('\n')
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', family_path, geo_path, spin_path, inv_path,
        *gates, *[folder/q for q in ['report.json', 'certificate.json', 'source-six-circuit.json.gz']],
        *[kernel_folder/q for q in ['report.json', 'matrix-certificate.json.gz', 'kernel-certificate.json.gz']],
        ROOT/'code/q_new_family_six_parameter_covariants_20261103_v63.py',
        ROOT/'code/q_new_family_source_kernel_20261103_v63.py']
    result = dict(id=STEM, all_passed=True,
        actual_six_parameter_raw_jet_full_coefficients_and_both_accepted_gauges=True,
        independent_rational_frame_primitive24_and_beta16=True,
        all16_scalar_exponents_and_denominators_domain_verified=True,
        all1025_circuit_nodes_bound_by_exact_beta_ancestors_and_C_D_word_recipes=True,
        universal_full_W16_pencil_change_by_independent_integer_coefficients=True,
        actual_N12_m_and_spectral_discriminant_circuit_bindings_accepted=True,
        independent_original30x26_matrix_and_all20_QQ_pivot_steps=True,
        independent_explicit_nonzero_residual5minor_and_full26_polynomial_kernel=True,
        actual_source_kernel_generic_rank25_unique_direction_and_du_zero=True,
        native_kernel_worker_not_replayed_or_assumed_for_rank_bound=True,
        no_new_point_high_degree_probe_or_full_actual_beta_quotient_field_claim=True,
        five_rational_first_integrals_not_in_this_receipt=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-certificate.json.gz': sha(out/'independent-certificate.json.gz')},
        seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'actual_source_kernel_generic_rank25_unique_direction_and_du_zero', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())

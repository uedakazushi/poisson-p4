#!/usr/bin/env python3
"""Exact Clifford model and two quadratic adjoint covariants; no new samples."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

from weight1_structure import COORDS, J, constraints, infinitesimal_action, sp_basis
from q_new_family_degree16_invariant_probe_20261103_v63 import character, echelon

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63'


def read(path):
    with path.open() as f:
        return json.load(f)


def matrix(rows):
    return s.Matrix([[s.Rational(x) for x in row] for row in rows])


def rows(M):
    return [[str(x) for x in row] for row in M.tolist()]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def rational_rank_pivots(M):
    _, pivots = DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref()
    return list(pivots)


def trace_product(A, B, monomials):
    """Full QQ quartic coefficient dictionary, not point evaluations."""
    result = defaultdict(lambda: s.S.Zero)
    for row in range(4):
        for col in range(4):
            aa = [(m, A[row * 4 + col, i]) for i, m in enumerate(monomials)
                  if A[row * 4 + col, i]]
            bb = [(m, B[col * 4 + row, i]) for i, m in enumerate(monomials)
                  if B[col * 4 + row, i]]
            for m, x in aa:
                for n, y in bb:
                    result[tuple(sorted(m + n))] += x * y
    return {k: v for k, v in result.items() if v}


def proportional(A, B):
    k = next(iter(B))
    ratio = A.get(k, 0) / B[k]
    assert ratio and all(A.get(m, 0) == ratio * B.get(m, 0)
                         for m in A.keys() | B.keys())
    return ratio


def determinant_pencil(C, D, monomials, F):
    """All QQ coefficients; base9 stores exponents <=8 without carrying."""
    powers = [9 ** i for i in range(16)]
    entries = []
    for rc, rd in zip(C.tolist(), D.tolist()):
        e = {}
        for j, (c, d) in enumerate(zip(rc, rd)):
            m = sum(powers[i] for i in monomials[j])
            if c:
                e[m, 0] = Fraction(str(c))
            if d:
                e[m, 1] = Fraction(str(d))
        entries.append(e)
    def multiply(a, b):
        result = defaultdict(Fraction)
        for (m, i), x in a.items():
            for (n, j), y in b.items():
                result[m + n, i + j] += x * y
        return {k: v for k, v in result.items() if v}
    determinant = defaultdict(Fraction)
    for perm in itertools.permutations(range(4)):
        sign = (-1) ** sum(perm[i] > perm[j] for i in range(4) for j in range(i + 1, 4))
        p = {(0, 0): Fraction(sign)}
        for i in range(4):
            p = multiply(p, entries[4 * i + perm[i]])
        for k, v in p.items():
            determinant[k] += v
    determinant = {k: v for k, v in determinant.items() if v}
    coefficients = [{m: x for (m, j), x in determinant.items() if j == i} for i in range(5)]
    ff = {(sum(powers[i] for i in m), 0): Fraction(str(v)) for m, v in F.items()}
    square = {m: v for (m, _), v in multiply(ff, ff).items()}
    assert all(square.get(m, 0) == 400 * coefficients[0].get(m, 0) + 1600 * coefficients[4].get(m, 0)
               for m in square.keys() | coefficients[0].keys() | coefficients[4].keys())
    basis, selected = {}, []
    for m in sorted(set().union(*(c.keys() for c in coefficients))):
        row = [c.get(m, Fraction(0)) for c in coefficients]
        for p, vector in sorted(basis.items()):
            z = row[p]
            if z:
                row = [x - z * y for x, y in zip(row, vector)]
        p = next((j for j, v in enumerate(row) if v), None)
        if p is not None:
            z = row[p]
            basis[p] = [v / z for v in row]
            selected.append(m)
    assert len(basis) == 5
    minor = s.Matrix([[s.Rational(c.get(m, 0)) for c in coefficients] for m in selected])
    assert minor.det() == -s.Rational(1, 16777216)
    def decode(key):
        exponents = []
        for _ in range(16):
            exponents.append(key % 9)
            key //= 9
        assert key == 0 and sum(exponents) == 8
        return exponents
    return dict(coefficient_order='p_i = coefficient t^(4-i) u^i in det(t C+u D)',
                polynomial_encoding='integer=sum exponent_i*9^i; exponent_i<=8',
                determinant_coefficients=[[[m, str(v)] for m, v in sorted(c.items())] for c in coefficients],
                coefficient_term_counts=list(map(len, coefficients)),
                selected_monomial_exponents=list(map(decode, selected)),
                selected_coefficient_minor=rows(minor), selected_coefficient_minor_determinant=str(minor.det()),
                full_polynomial_identity='F4^2 = 400 p_0 + 1600 p_4',
                all_degree8_invariants_spanned=True)


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    out.mkdir(parents=True, exist_ok=False)
    ip = ROOT / 'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    orbitpath = ROOT / 'reports/weight1_structure_20261030.json'
    oldpaths = {d: ROOT / f'reports/Q_NEW_FAMILY_DEGREE{d}_INVARIANT_PROBE_20261103_v63/report.json'
                for d in (16, 18, 20)}
    data = read(ip)
    B = matrix(data['primitive24_by16_basis'])
    rho = [matrix(A) for A in data['all_ten_representation_matrices']]
    omega = matrix(data['symplectic_form16'])
    gs = sp_basis()
    ps = list(itertools.combinations(range(4), 2))[:-1]
    alpha = []
    for e, (i, j) in enumerate(ps):
        a = s.zeros(4)
        a[i, j], a[j, i] = 1, -1
        if e == 0:
            a[2, 3], a[3, 2] = -1, 1
        alpha.append(a)
    gamma = [a * J.inv() for a in alpha]
    Q = s.Matrix(5, 5, lambda i, j: s.trace(gamma[i] * gamma[j]) / 4)
    assert Q.det() == s.Rational(1, 16)
    assert all(gamma[i] * gamma[j] + gamma[j] * gamma[i] == 2 * Q[i, j] * s.eye(4)
               for i in range(5) for j in range(5))
    CT = s.Matrix.hstack(*gamma)
    assert len(rational_rank_pivots(CT)) == 4
    P = [s.Matrix(5, 4, lambda e, k: B[COORDS.index((ps[e], k)), l]) for l in range(16)]
    embedding = s.Matrix.hstack(*(s.Matrix.vstack(*(p[e, :].T for e in range(5))) for p in P))
    assert CT * embedding == s.zeros(4, 16)
    assert len(rational_rank_pivots(embedding)) == 16
    d, tr = constraints()
    assert (d.col_join(tr)) * B == s.zeros(8, 16)
    assert len(rational_rank_pivots(B)) == 16
    R5 = [s.Matrix(5, 5, lambda i, j: (a.T * alpha[j] + alpha[j] * a)[ps[i][0], ps[i][1]])
          for a in gs]
    for a, rr, r in zip(gs, R5, rho):
        assert rr.T * Q + Q * rr == s.zeros(5)
        assert r.T * omega + omega * r == s.zeros(16)
        for e in range(5):
            assert a.T * gamma[e] - gamma[e] * a.T == sum(
                (rr[f, e] * gamma[f] for f in range(5)), s.zeros(4))
        for l in range(16):
            lhs = sum((r[k, l] * P[k] for k in range(16)), s.zeros(5, 4))
            assert lhs == rr * P[l] + P[l] * a
            assert B * r[:, l] == infinitesimal_action(B[:, l], a)
    rflat = s.Matrix.hstack(*(a.reshape(25, 1) for a in R5))
    selected = rational_rank_pivots(rflat.T)
    assert len(selected) == 10
    lift = rflat[selected, :].inv()
    monomials = list(itertools.combinations_with_replacement(range(16), 2))
    C, D, V = [], [], []
    for i, j in monomials:
        c = P[i].T * Q * P[j]
        v = P[i] * J.inv() * P[j].T
        if i != j:
            c += P[j].T * Q * P[i]
            v += P[j] * J.inv() * P[i].T
        c, v = c * J.inv(), v * Q
        assert c.T * J.inv() + J.inv() * c == s.zeros(4)
        assert v.T * Q + Q * v == s.zeros(5)
        co = lift * v.reshape(25, 1)[selected, :]
        assert rflat * co == v.reshape(25, 1)
        ds = sum((x * a.T for x, a in zip(co, gs)), s.zeros(4))
        C.append(c)
        D.append(ds)
        V.append(v)
    C = s.Matrix.hstack(*(c.reshape(16, 1) for c in C))
    D = s.Matrix.hstack(*(c.reshape(16, 1) for c in D))
    pair = C.reshape(2176, 1).row_join(D.reshape(2176, 1))
    witness_rows = rational_rank_pivots(pair.T)
    assert len(witness_rows) == 2
    witness = pair[witness_rows, :]
    gram = s.Matrix([[s.trace(a * b) for b in gs] for a in gs])
    assert all(s.denom(x) % 32003 for M in [omega, omega.inv(), gram.inv(), *rho]
               for x in M)
    moment_components = s.Matrix(10, 136, lambda k, m:
        (omega * rho[k])[monomials[m][0], monomials[m][1]] /
        (2 if monomials[m][0] == monomials[m][1] else 1))
    assert all(omega * r == (omega * r).T for r in rho)
    moment = s.Matrix.hstack(*(a.T.reshape(16, 1) for a in gs)) * gram.inv() * moment_components
    target = moment.reshape(2176, 1)
    combination = witness.inv() * target[witness_rows, :]
    assert pair * combination == target
    dynkin = s.Matrix([[s.trace(a * b) for b in rho] for a in rho])
    assert dynkin == 12 * gram
    # All quartic coefficients: this also binds the older unique-quartic normalization.
    F = trace_product(moment, moment, monomials)
    oldF = {tuple(m): s.Rational(v) for m, v in data['unique_quartic_sparse']}
    quartic_old_ratio = proportional(F, oldF)
    trace_ratios = [proportional(trace_product(A, E, monomials), F)
                    for A, E in ((C, C), (C, D), (D, D))]
    weights = list(map(tuple, data['Cartan_weights']))
    sym2 = Counter((weights[i][0] + weights[j][0], weights[i][1] + weights[j][1])
                   for i in range(16) for j in range(i, 16))
    weyl = []
    for perm in ((0, 1), (1, 0)):
        for signs in itertools.product((-1, 1), repeat=2):
            wr = (signs[0] * (2, 1)[perm[0]], signs[1] * (2, 1)[perm[1]])
            weyl.append((wr, signs[0] * signs[1] * (-1 if perm == (1, 0) else 1)))
    def multiplicity(ch, lam):
        return sum(sign * ch[(lam[0] + 2 - wr[0], lam[1] + 1 - wr[1])]
                   for wr, sign in weyl)
    def dimension(lam):
        x, y = lam
        return (x - y + 1) * (y + 1) * (x + 2) * (x + y + 3) // 6
    assert multiplicity(Counter(weights), (2, 1)) == 1 and dimension((2, 1)) == 16
    decomposition = []
    for lam, expected in (((4, 2), 1), ((3, 1), 1), ((2, 0), 2)):
        m = multiplicity(sym2, lam)
        assert m == expected
        decomposition.append(dict(highest_epsilon=list(lam), multiplicity=m, dimension=dimension(lam)))
    assert sum(x['multiplicity'] * x['dimension'] for x in decomposition) == 136
    chars = character(weights, 20)
    dimensions = [x['invariant_dimension'] for x in chars]
    assert {i: x for i, x in enumerate(dimensions) if x} == {0: 1, 4: 1, 8: 5, 12: 9, 14: 3, 16: 20, 18: 9, 20: 36}
    pencil = determinant_pencil(C, D, monomials, F)
    product_receipts = {}
    for degree, start, expected in ((16, 36, 19), (18, 60, 3), (20, 108, 36)):
        old = read(oldpaths[degree])
        ev = [row[start:] for row in old['ambient_evaluation']]
        pr, pc, _ = echelon(ev)
        assert len(pr) == expected
        product_receipts[degree] = dict(rank_mod32003=expected, input_start_column=start,
                                       selected_rows=pr, selected_product_columns=pc,
                                       no_new_evaluations=True)
    # Monomial upper counts using the already necessary minimal generators.
    generators = {4: 1, 8: 4, 12: 4, 14: 3, 16: 1, 18: 6}
    dp = [0] * 21
    dp[0] = 1
    for degree, number in generators.items():
        for _ in range(number):
            for k in range(degree, 21):
                dp[k] += dp[k - degree]
    assert dp == dimensions
    oldorbit = read(orbitpath)
    assert oldorbit['all_passed'] and oldorbit['generic_sp4_orbit_rank_witness']['rank'] == 10
    ow = oldorbit['generic_sp4_orbit_rank_witness']
    assert matrix(ow['matrix']).det() == s.Rational(ow['determinant']) != 0
    cp = out / 'certificate.json.gz'
    cert = dict(primitive_form_pairs=ps, primitive_form_matrices=list(map(rows, alpha)),
                gamma_matrices=list(map(rows, gamma)), quadratic_form5=rows(Q),
                gamma_trace4x20=rows(CT), vector_spinor_embedding20x16=rows(embedding),
                action5_matrices=list(map(rows, R5)), lie_lift_selected_rows=selected,
                lie_lift_inverse=rows(lift), quadratic_monomials=monomials,
                spin_adjoint_C_coefficients16x136=rows(C), vector_adjoint_D_coefficients16x136=rows(D),
                moment_transpose_coefficients16x136=rows(moment),
                independent_pair_rows=witness_rows, independent_pair_minor=rows(witness),
                independent_pair_minor_determinant=str(witness.det()),
                moment_transpose_combination_C_D=[str(x) for x in combination],
                quartic_trace_moment_square=[[list(m), str(v)] for m, v in sorted(F.items())],
                quartic_trace_ratios_C2_CD_D2=[str(x) for x in trace_ratios],
                moment_trace_square_over_old_unique_quartic=str(quartic_old_ratio),
                canonical_adjoint_pencil=pencil,
                sym2_decomposition=decomposition, character_records=chars,
                existing_lower_product_rank_receipts=product_receipts,
                reused_generic_orbit_rank=10, affine_invariant_quotient_dimension_upper_bound=6)
    with gzip.open(cp, 'wt') as f:
        json.dump(cert, f, indent=2)
        f.write('\n')
    inputs = [Path(__file__), ROOT / 'code/weight1_structure.py',
              ROOT / 'code/q_new_family_degree16_invariant_probe_20261103_v63.py', ip, orbitpath, *oldpaths.values()]
    result = dict(id=STEM, all_passed=True, exact_field='QQ',
                  all_25_Clifford_identities=True, gamma_kernel_dimension=16,
                  all_10_times_16_original_representation_intertwinings=True,
                  independent_quadratic_adjoint_covariants=2,
                  full_sym2_adjoint_multiplicity=2,
                  moment_transpose_combination_C_D=[str(x) for x in combination],
                  trace_index_rho_over_standard4=12,
                  quartic_trace_ratios_C2_CD_D2=[str(x) for x in trace_ratios],
                  moment_trace_square_over_old_unique_quartic=str(quartic_old_ratio),
                  determinant_pencil_coefficients_span_all_Inv8=True,
                  determinant_pencil_coefficient_minor=pencil['selected_coefficient_minor_determinant'],
                  determinant_pencil_identity='F4^2 = 400 p_0 + 1600 p_4',
                  invariant_dimensions_through20=dimensions,
                  minimal_generator_numbers_through20={**generators, 20: 0},
                  no_relations_between_minimal_generators_through20=True,
                  affine_invariant_quotient_dimension_upper_bound=6,
                  invariant_ring_not_polynomial=True,
                  negative_polar_argument='A finite-reflection Cartan restriction would force a polynomial invariant ring.',
                  no_new_point_evaluations=True, no_degree22_or_higher_probe=True,
                  NEW25_equation_or_discriminant_candidate_claimed=False,
                  certificate_sha256=sha(cp),
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs},
                  seconds=time.monotonic() - started)
    (out / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'moment_transpose_combination_C_D',
          'quartic_trace_ratios_C2_CD_D2', 'moment_trace_square_over_old_unique_quartic', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output.resolve())

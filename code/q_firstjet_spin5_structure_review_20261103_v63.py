#!/usr/bin/env python3
"""Independent Fraction/cleared-integer coefficient review of the Spin5 model."""
import argparse
from collections import Counter, defaultdict
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

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63_REVIEW'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def zero(n, m):
    return [[F(0) for _ in range(m)] for _ in range(n)]


def mat(x):
    return [[F(q) for q in row] for row in x]


def tr(A):
    return sum(A[i][i] for i in range(len(A)))


def tp(A):
    return list(map(list, zip(*A)))


def mm(A, B):
    BT = tp(B)
    return [[sum((x*y for x, y in zip(row, col) if x and y), F(0)) for col in BT] for row in A]


def add(A, B, c=1):
    return [[x+c*y for x, y in zip(a, b)] for a, b in zip(A, B)]


def scale(A, c):
    return [[c*x for x in row] for row in A]


def inv(A):
    return mat(s.Matrix(A).inv().tolist())


def plain(A):
    return [[str(q) for q in row] for row in A]


def polyadd(A, B, c=1):
    out = dict(A)
    for k, v in B.items():
        out[k] = out.get(k, 0)+c*v
        if not out[k]:
            del out[k]
    return out


def polymul(A, B):
    out = defaultdict(int)
    for m, a in A.items():
        for n, b in B.items():
            out[m+n] += a*b
    return {m: a for m, a in out.items() if a}


def polymm(A, B):
    n = len(A)
    return [[sum_polys([polymul(A[i][k], B[k][j]) for k in range(n)])
             for j in range(n)] for i in range(n)]


def sum_polys(values):
    out = {}
    for p in values:
        out = polyadd(out, p)
    return out


def modrank(A, p=32003):
    A = [[int(x) % p for x in row] for row in A]
    r = 0
    for j in range(len(A[0])):
        i = next((i for i in range(r, len(A)) if A[i][j]), None)
        if i is None:
            continue
        A[r], A[i] = A[i], A[r]
        a = pow(A[r][j], -1, p)
        A[r] = [a*x % p for x in A[r]]
        for i in range(r+1, len(A)):
            a = A[i][j]
            if a:
                A[i] = [(x-a*y) % p for x, y in zip(A[i], A[r])]
        r += 1
        if r == len(A):
            break
    return r


def main(out):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    author = ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/report.json'
    cp = author.parent/'certificate.json.gz'
    data_path = ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    data, cert, ar = map(read, [data_path, cp, author])
    assert ar['all_passed'] and ar['certificate_sha256'] == sha(cp)
    for name, digest in ar['input_sha256'].items():
        assert sha(ROOT/name) == digest, name
    pairs = list(itertools.combinations(range(4), 2))
    coords = [(ij, k) for ij in pairs for k in range(4)]
    B = mat(data['primitive24_by16_basis'])
    rho = list(map(mat, data['all_ten_representation_matrices']))
    omega = mat(data['symplectic_form16'])
    J = mat([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])
    Ji = scale(J, -1)
    generators = []
    for i in range(4):
        for j in range(i, 4):
            A = zero(4, 4)
            A[i][j] = A[j][i] = F(1)
            generators.append(mm(J, A))
    alpha = []
    for i, j in pairs[:-1]:
        a = zero(4, 4)
        a[i][j], a[j][i] = F(1), F(-1)
        if (i, j) == (0, 1):
            a[2][3], a[3][2] = F(-1), F(1)
        alpha.append(a)
    gamma = [mm(a, Ji) for a in alpha]
    Q = [[tr(mm(a, b))/4 for b in gamma] for a in gamma]
    assert s.Matrix(Q).det() == s.Rational(1, 16)
    I4 = [[F(i == j) for j in range(4)] for i in range(4)]
    for i in range(5):
        for j in range(5):
            assert add(mm(gamma[i], gamma[j]), mm(gamma[j], gamma[i])) == scale(I4, 2*Q[i][j])
    assert list(map(plain, alpha)) == cert['primitive_form_matrices']
    assert list(map(plain, gamma)) == cert['gamma_matrices'] and plain(Q) == cert['quadratic_form5']
    P = [[[B[coords.index((pairs[e], k))][l] for k in range(4)] for e in range(5)] for l in range(16)]
    trace_map = s.Matrix.hstack(*map(s.Matrix, gamma))
    emb = s.Matrix.hstack(*(s.Matrix([p[e][k] for e in range(5) for k in range(4)]) for p in P))
    assert trace_map.rank() == 4 and emb.rank() == 16 and trace_map*emb == s.zeros(4, 16)
    assert plain(emb.tolist()) == cert['vector_spinor_embedding20x16']
    assert s.Matrix(B).rank() == 16
    for l in range(16):
        forms = []
        for k in range(4):
            a = zero(4, 4)
            for ij in pairs:
                a[ij[0]][ij[1]] = B[coords.index((ij, k))][l]
                a[ij[1]][ij[0]] = -a[ij[0]][ij[1]]
            forms.append(a)
        assert all(forms[k][0][1]+forms[k][2][3] == 0 for k in range(4))
        for i, j, k in itertools.combinations(range(4), 3):
            assert forms[i][j][k]-forms[j][i][k]+forms[k][i][j] == 0
        for a, r in zip(generators, rho):
            physical = []
            for k in range(4):
                value = add(mm(tp(a), forms[k]), mm(forms[k], a))
                for i in range(4):
                    value = add(value, forms[i], a[i][k])
                physical.append(value)
            original = [physical[k][ij[0]][ij[1]] for ij, k in coords]
            assert original == [sum(B[i][j]*r[j][l] for j in range(16)) for i in range(24)]
    R5 = []
    for a, r in zip(generators, rho):
        rr = zero(5, 5)
        for j in range(5):
            acted = add(mm(tp(a), alpha[j]), mm(alpha[j], a))
            for i, ij in enumerate(pairs[:-1]):
                rr[i][j] = acted[ij[0]][ij[1]]
        assert add(mm(tp(rr), Q), mm(Q, rr)) == zero(5, 5)
        assert add(mm(tp(r), omega), mm(omega, r)) == zero(16, 16)
        for e in range(5):
            right = zero(4, 4)
            for f in range(5):
                right = add(right, gamma[f], rr[f][e])
            assert add(mm(tp(a), gamma[e]), mm(gamma[e], tp(a)), -1) == right
        for l in range(16):
            lhs = zero(5, 4)
            for k in range(16):
                lhs = add(lhs, P[k], r[k][l])
            assert lhs == add(mm(rr, P[l]), mm(P[l], a))
        R5.append(rr)
    assert list(map(plain, R5)) == cert['action5_matrices']
    flat = s.Matrix.hstack(*(s.Matrix([q for row in a for q in row]) for a in R5))
    selected = cert['lie_lift_selected_rows']
    inverse = inv(flat[selected, :].tolist())
    assert plain(inverse) == cert['lie_lift_inverse']
    monomials = list(itertools.combinations_with_replacement(range(16), 2))
    C, D = [], []
    for i, j in monomials:
        cc = mm(mm(tp(P[i]), Q), P[j])
        vv = mm(mm(P[i], Ji), tp(P[j]))
        if i != j:
            cc = add(cc, mm(mm(tp(P[j]), Q), P[i]))
            vv = add(vv, mm(mm(P[j], Ji), tp(P[i])))
        cc, vv = mm(cc, Ji), mm(vv, Q)
        # Explicit ordered selection avoids assuming pivot rows are sorted.
        co = mm(inverse, [[sum(vv, [])[k]] for k in selected])
        assert mm(mat(flat.tolist()), co) == [[q] for q in sum(vv, [])]
        dd = zero(4, 4)
        for q, a in zip(co, generators):
            dd = add(dd, tp(a), q[0])
        assert add(mm(tp(cc), Ji), mm(Ji, cc)) == zero(4, 4)
        assert add(mm(tp(dd), Ji), mm(Ji, dd)) == zero(4, 4)
        C.append(sum(cc, []))
        D.append(sum(dd, []))
    C, D = tp(C), tp(D)
    assert plain(C) == cert['spin_adjoint_C_coefficients16x136']
    assert plain(D) == cert['vector_adjoint_D_coefficients16x136']
    gram = [[tr(mm(a, b)) for b in generators] for a in generators]
    gi = inv(gram)
    components = []
    for r in rho:
        T = mm(omega, r)
        assert T == tp(T)
        components.append([T[i][j]/(2 if i == j else 1) for i, j in monomials])
    Aflat = [[tp(a)[i][j] for a in generators] for i in range(4) for j in range(4)]
    moment = mm(mm(Aflat, gi), components)
    assert moment == add(scale(C, -1), scale(D, 2))
    assert plain(moment) == cert['moment_transpose_coefficients16x136']
    assert [[tr(mm(a, b)) for b in rho] for a in rho] == scale(gram, 12)
    pairminor = [[sum(C, [])[i], sum(D, [])[i]] for i in cert['independent_pair_rows']]
    assert s.Matrix(pairminor).det() == s.Rational(cert['independent_pair_minor_determinant']) != 0
    # Entire determinant pencil by a DIFFERENT engine: Newton's formula for
    # traceless4x4 matrices, using common-denominator-cleared integer polys.
    den = math.lcm(*(q.denominator for A in [C, D] for row in A for q in row))
    powers = [9**i for i in range(16)]
    keys = [powers[i]+powers[j] for i, j in monomials]
    def polys(A):
        return [[{keys[k]: int(A[4*i+j][k]*den) for k in range(136) if A[4*i+j][k]}
                 for j in range(4)] for i in range(4)]
    CC, DD = polys(C), polys(D)
    C2, D2 = polymm(CC, CC), polymm(DD, DD)
    CD, DC = polymm(CC, DD), polymm(DD, CC)
    cross = [[polyadd(CD[i][j], DC[i][j]) for j in range(4)] for i in range(4)]
    matrices = [C2, cross, D2]
    traces = [sum_polys([A[i][i] for i in range(4)]) for A in matrices]
    Fnum = sum_polys([polyadd(polyadd(C2[i][i], cross[i][i], -2), D2[i][i], 4) for i in range(4)])
    # mu=-C+2D, so mu²=C²-2(CD+DC)+4D².
    Fdict = {m: F(v, den*den) for m, v in Fnum.items()}
    def encoded_quartic(records):
        return {sum(powers[i] for i in m): F(v) for m, v in records}
    assert Fdict == encoded_quartic(cert['quartic_trace_moment_square'])
    oldF = encoded_quartic(data['unique_quartic_sparse'])
    assert Fdict == {m: F(5, 2)*v for m, v in oldF.items()}
    for A, ratio in [(C2, F(1, 5)), (CD, F(-1, 10)), (D2, F(1, 10))]:
        tt = sum_polys([A[i][i] for i in range(4)])
        assert {m: F(v, den*den) for m, v in tt.items()} == {m: ratio*v for m, v in Fdict.items()}
    coeffs = []
    for degree in range(5):
        a, b = {}, {}
        for i in range(3):
            j = degree-i
            if not 0 <= j < 3:
                continue
            a = polyadd(a, polymul(traces[i], traces[j]))
            for row in range(4):
                for col in range(4):
                    b = polyadd(b, polymul(matrices[i][row][col], matrices[j][col][row]))
        numerator = polyadd(a, b, -2)
        coeffs.append({m: F(v, 8*den**4) for m, v in numerator.items()})
    pencil = cert['canonical_adjoint_pencil']
    saved = [{int(m): F(v) for m, v in pp} for pp in pencil['determinant_coefficients']]
    assert coeffs == saved
    Fsq = defaultdict(F)
    for m, a in Fdict.items():
        for n, b in Fdict.items():
            Fsq[m+n] += a*b
    assert {m: a for m, a in Fsq.items() if a} == polyadd(
        {m: 400*a for m, a in coeffs[0].items()}, {m: 1600*a for m, a in coeffs[4].items()})
    selected_monomials = [sum(e*p for e, p in zip(exps, powers)) for exps in pencil['selected_monomial_exponents']]
    pm = [[p.get(m, F(0)) for p in coeffs] for m in selected_monomials]
    assert s.Matrix(pm).det() == -s.Rational(1, 2**24)
    assert plain(pm) == pencil['selected_coefficient_minor']
    # Independent individual-weight symmetric-power DP and all8 Weyl terms.
    weights = [tuple(w) for w in data['Cartan_weights']]
    dp = [Counter() for _ in range(21)]
    dp[0][(0, 0)] = 1
    for wgt in weights:
        for degree in range(1, 21):
            for (a, b), n in dp[degree-1].items():
                dp[degree][a+wgt[0], b+wgt[1]] += n
    weyl = []
    for perm in [(0, 1), (1, 0)]:
        for signs in itertools.product([-1, 1], repeat=2):
            weyl.append(((signs[0]*(2, 1)[perm[0]], signs[1]*(2, 1)[perm[1]]),
                         signs[0]*signs[1]*(-1 if perm == (1, 0) else 1)))
    mult = lambda ch, lam: sum(sign*ch[(lam[0]+2-a, lam[1]+1-b)] for (a, b), sign in weyl)
    dimensions = []
    for degree, ch in enumerate(dp):
        assert sum(ch.values()) == math.comb(15+degree, degree)
        dimensions.append(mult(ch, (0, 0)))
    assert dimensions == ar['invariant_dimensions_through20']
    assert [mult(dp[2], a) for a in [(4, 2), (3, 1), (2, 0)]] == [1, 1, 2]
    assert mult(dp[1], (2, 1)) == 1 and 81+35+2*10 == 136
    product_ranks = {}
    product_paths = []
    for degree, col, expected in [(16, 36, 19), (18, 60, 3), (20, 108, 36)]:
        p = ROOT/f'reports/Q_NEW_FAMILY_DEGREE{degree}_INVARIANT_PROBE_20261103_v63/report.json'
        product_paths.append(p)
        old = read(p)
        rank = modrank([row[col:] for row in old['ambient_evaluation']])
        assert rank == expected
        product_ranks[degree] = rank
    counts = [0]*21
    counts[0] = 1
    gens = {4: 1, 8: 4, 12: 4, 14: 3, 16: 1, 18: 6}
    for degree, number in gens.items():
        for _ in range(number):
            for k in range(degree, 21):
                counts[k] += counts[k-degree]
    assert counts == dimensions
    orbitpath = ROOT/'reports/weight1_structure_20261030.json'
    ow = read(orbitpath)['generic_sp4_orbit_rank_witness']
    assert s.Matrix(ow['matrix']).det() == s.Rational(ow['determinant']) != 0
    prior = [ROOT/'reports/Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW.json',
             ROOT/'reports/Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in prior)
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', author, cp, data_path,
             ROOT/'research/notes/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63.md',
             orbitpath, *product_paths, *prior]
    result = dict(
        id=STEM, all_passed=True,
        independent_Fraction_Clifford25_and_original10_by16_coefficient_actions=True,
        full_gamma_kernel20_to4_and_original_closed_primitive_constraints=True,
        all_quadratic_C_D_moment_coefficients_independently_rebuilt=True,
        moment_transpose_equals_minus_C_plus2D=True, all100_trace_indices12=True,
        independent_full_QQ_quartic_trace_and_old_normalization_identities=True,
        independent_integer_Newton_pencil_all_coefficients_match=True,
        determinant_engine='det(P)=((tr(P²))²-2tr(P⁴))/8; all sparse integer coefficients, no point evaluation',
        common_quadratic_coefficient_denominator=den,
        pencil_term_counts=[len(p) for p in coeffs], independent_pencil_minor='-1/16777216',
        pencil_spans_all_degree8_invariants=True,
        sym2_adjoint_multiplicity2_and_full_C_D_covariance_accepted=True,
        all_character_dimensions_through20_independently_rebuilt=dimensions,
        reused_old_product_ranks_independently_reduced=product_ranks,
        minimal_generator_counts_and_no_relations_through20_accepted=True,
        invariant_ring_not_polynomial_by_atleast9_generators_and_dimension_atmost6=True,
        finite_reflection_polynomial_Cartan_restriction_cannot_apply_as_stated=True,
        no_new_samples_new_hypersurface_discriminant_or_degree22_search=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        seconds=time.monotonic()-start,
    )
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'pencil_term_counts',
        'independent_pencil_minor', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())

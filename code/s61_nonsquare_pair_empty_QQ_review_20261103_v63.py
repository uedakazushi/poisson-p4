#!/usr/bin/env python3
"""Independent coefficient-arithmetic audit of the whole QQ pair gate.

Only the previously independent actual-source/finite review is replayed.
All new identities use Fraction dictionaries and elementary univariate
long division.  No resultant, factorization, GB, finite field or field
structure on QQ[v]/Q23 is used.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import tarfile
import time

from s61_nonsquare_pair_finite_bad_locus_review_20261103_v63 import sha

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63_REVIEW'


def read(p):
    return json.loads(p.read_text())


def decode(terms, n):
    ans = {}
    for ex, c in terms:
        ex = tuple(ex)
        assert len(ex) == n and all(isinstance(k, int) and k >= 0 for k in ex)
        assert ex not in ans
        if F(c):
            ans[ex] = F(c)
    return ans


def add(p, q):
    r = p.copy()
    for ex, c in q.items():
        r[ex] = r.get(ex, F(0)) + c
        if not r[ex]:
            del r[ex]
    return r


def scale(p, c):
    return {ex: c*d for ex, d in p.items() if c*d}


def mul(p, q):
    r = {}
    for ex, c in p.items():
        for ey, d in q.items():
            mon = tuple(i+j for i, j in zip(ex, ey))
            r[mon] = r.get(mon, F(0)) + c*d
    return {mon: c for mon, c in r.items() if c}


def power(p, n, dim):
    r = {(0,)*dim: F(1)}
    for _ in range(n):
        r = mul(r, p)
    return r


def rem(p, q):
    """QQ[v] long division; q need not be irreducible or squarefree."""
    assert q and all(len(ex) == 1 for ex in q)
    n = max(ex[0] for ex in q)
    assert n > 0
    r = p.copy()
    while r:
        k = max(ex[0] for ex in r)
        if k < n:
            break
        c = r[(k,)] / q[(n,)]
        for (i,), d in q.items():
            key = (i+k-n,)
            r[key] = r.get(key, F(0)) - c*d
            if not r[key]:
                del r[key]
    return r


def in_v(p):
    assert all(ex[0] == 0 for ex in p)
    return {(ex[1],): c for ex, c in p.items()}


def in_Hv(p):
    return {(0, ex[0]): c for ex, c in p.items()}


def subtract(p, q):
    return add(p, scale(q, F(-1)))


def main(out):
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(120)
    out.mkdir(parents=True, exist_ok=False)
    parent = ROOT/'reports/S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63'
    author = read(parent/'report.json')
    assert author['frozen'] and author['all_exact_QQ_checks_passed']
    archive = ROOT/author['archive']
    assert sha(archive) == author['archive_sha256']
    assert archive.stat().st_size == author['archive_bytes'] == 96138
    snapshot = out/'immutable-input'
    snapshot.mkdir()
    with tarfile.open(archive) as tar:
        for member in tar:
            assert member.isfile()
            dst = snapshot/member.name
            assert dst.resolve().is_relative_to(snapshot.resolve())
            dst.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as src, dst.open('wb') as target:
                shutil.copyfileobj(src, target)
    for rel, expected in author['input_sha256'].items():
        p = snapshot/rel
        if not p.exists():
            p = ROOT/rel
        assert sha(p) == expected, rel

    # Reuse our own independent 804-term/source proof, not an author's helper.
    finite_code = ROOT/'code/s61_nonsquare_pair_finite_bad_locus_review_20261103_v63.py'
    finite_receipt = ROOT/'reports/S61_NONSQUARE_PAIR_FINITE_BAD_LOCUS_20261103_v63_REVIEW/report.json'
    accepted = read(finite_receipt)
    assert accepted['all_passed']
    for rel, expected in accepted['input_sha256'].items():
        assert sha(ROOT/rel) == expected, rel
    assert sha(ROOT/'reports/S61_NONSQUARE_PAIR_FINITE_BAD_LOCUS_20261103_v63/actual_pair_finite_QQ_exact.tar.gz') == author['finite_theorem_frozen_archive_sha256']
    with (out/'actual804-independent-replay.log').open('w') as log:
        subprocess.run([sys.executable, str(finite_code), '--output', str(out/'actual804-independent-replay')], check=True, timeout=60, stdout=log)
    assert read(out/'actual804-independent-replay/report.json')['all_passed']

    source = read(snapshot/'work/nonsquare-v62/pair-finite-proof-v76/certificate.json')
    data = read(snapshot/'work/nonsquare-v62/pair-finite-subresultant-v78/whole-pair-cleanup-certificate.json')
    eqs = [decode(q, 3) for q in source['linear_t_equations']]  # exponents t,v,H
    ab = []
    for q in eqs:
        assert all(ex[0] <= 1 for ex in q)
        ab.append(tuple({(ex[2], ex[1]): c for ex, c in q.items() if ex[0] == k} for k in [1, 0]))
    N01, A, B = [decode(q, 2) for q in source['elimination_polynomials']]
    R = decode(source['projection_certificate']['R'], 2)
    f = {(2,): F(3), (1,): F(13), (0,): F(12)}
    cubic, large = data['cubic_branch'], data['degree23_branch']
    C, Q23 = decode(cubic['C'], 1), decode(large['Q'], 1)
    assert C == {(3,): F(8), (2,): F(30), (1,): F(71), (0,): F(58)}
    assert max(ex[0] for ex in Q23) == 23
    product = mul(mul(power(f, 2, 1), power(C, 2, 1)), Q23)
    Rv = in_v(R)
    assert Rv and product and max(ex[0] for ex in Rv) == 33
    k = max(product)
    factor = Rv[k]/product[k]
    assert factor and scale(product, factor) == Rv

    # C=0: beta3=-f*C and t is a unit; alpha3=N01=0 is impossible.
    alpha, beta = ab[2]
    assert beta == scale(in_Hv(mul(f, C)), F(-1))
    U, V = [decode(q, 2) for q in cubic['resultant_Bezout_alpha4_N01']]
    RC = decode(cubic['alpha4_N01_resultant'], 2)
    assert add(mul(U, alpha), mul(V, N01)) == RC
    RCv = in_v(RC)
    RCinv = decode(cubic['resultant_inverse_mod_C'], 1)
    assert rem(mul(RCv, RCinv), C) == {(0,): F(1)}

    # Q23=0: all following equalities are in a finite algebra, not a field.
    S1 = decode(data['linear_subresultant'], 2)
    UA, UB = [decode(q, 2) for q in data['subresultant_A_B_Bezout']]
    assert add(mul(UA, A), mul(UB, B)) == S1
    assert max(ex[0] for ex in S1) == 1
    LC = {(ex[1],): c for ex, c in S1.items() if ex[0] == 1}
    cc = {(ex[1],): c for ex, c in S1.items() if ex[0] == 0}
    LCinv = decode(large['linear_leading_inverse'], 1)
    Hb = decode(large['Hbar'], 1)
    assert rem(mul(LC, LCinv), Q23) == {(0,): F(1)}
    assert not rem(add(Hb, mul(cc, LCinv)), Q23)

    def compose(poly, images):
        ps = []
        for i, image in enumerate(images):
            n = max((ex[i] for ex in poly), default=0)
            powers = [{(0,): F(1)}]
            for _ in range(n):
                powers.append(rem(mul(powers[-1], image), Q23))
            ps.append(powers)
        result = {}
        for ex, c in poly.items():
            term = {(0,): c}
            for powers, i in zip(ps, ex):
                term = rem(mul(term, powers[i]), Q23)
            result = add(result, term)
        return rem(result, Q23)

    v = {(1,): F(1)}
    ae = compose(alpha, [Hb, v])
    be = compose(beta, [Hb, v])
    ainv = decode(large['alpha4_inverse'], 1)
    Tb = decode(large['Tbar'], 1)
    assert rem(mul(ae, ainv), Q23) == {(0,): F(1)}
    assert not rem(add(mul(ae, Tb), be), Q23)
    assert not compose(A, [Hb, v]) and not compose(B, [Hb, v])
    assert not compose(N01, [Hb, v])
    for al, be in ab:
        assert not rem(add(mul(compose(al, [Hb, v]), Tb), compose(be, [Hb, v])), Q23)

    # This actual minor is rebuilt from all 12 certified physical coefficients.
    rows = [[decode(q, 3) for q in row] for row in source['natural_rows']]
    assert len(rows) == 2 and all(len(row) == 6 for row in rows)
    minor = subtract(mul(rows[0][5], rows[1][0]), mul(rows[1][5], rows[0][0]))
    actual = compose(minor, [Tb, v, Hb])
    assert actual == decode(large['remaining_pair_minor_col0'], 1)
    inverse = decode(data['degree23_remaining_minor_inverse'], 1)
    assert rem(mul(actual, inverse), Q23) == {(0,): F(1)}

    deps = [ROOT/f'reports/{stem}/report.json' for stem in [
        'S61_NONSQUARE_GCD_GATE_PROPORTIONAL_EXCLUSION_20261103_v63_REVIEW',
        'S61_NONSQUARE_FINITE_GCD_EXCEPTION_NAKAYAMA_20261103_v63_REVIEW',
        'S61_NONSQUARE_FOUR_QUADRIC_BINARY_REDUCTION_20261103_v63_REVIEW']]
    assert all(read(p)['all_passed'] for p in deps)
    evidence = dict(
        exact_factorization={'degree': 33, 'nonzero_QQ_scalar': str(factor), 'powers': {'f': 2, 'C3': 2, 'Q23': 1}},
        C_branch_literally_checked=['beta3=-f*C', 'U*alpha3+V*N01=Rc(v)', 'Rc*inverse=1 mod C'],
        Q23_branch_literally_checked=['UA*A+UB*B=L*H+K', 'L*inverse=1 mod Q23', 'Hbar=-K/L',
            'alpha3(Hbar)*inverse=1 mod Q23', 'tbar=-beta3/alpha3', 'all3 linear-t equations and all3 eliminants',
            'full actual col0/u2 minor evaluation', 'actual minor * inverse=1 mod Q23'],
        Hbar_degree=max(ex[0] for ex in Hb), tbar_degree=max(ex[0] for ex in Tb),
        actual_minor_term_count=len(minor), actual_minor_residue_degree=max(ex[0] for ex in actual),
        no_irreducibility_squarefreeness_field_flatness_or_modp_assumption=True,
        pointwise_to_module='Maximal-minor rank failure is empty over QQbar on aPJ; finite-type/local Nakayama gives the zero cokernel. Powers C^2 cause no omitted geometric points.',
        source_binding='Fresh independent 804-term reconstruction and all12 natural coefficients, known-unit row changes, f=0 exclusions and necessary-minor equations.'
    )
    (out/'independent-certificate.json').write_text(json.dumps(evidence, indent=2)+'\n')
    shutil.rmtree(snapshot)
    result = dict(id=STEM, all_passed=True,
        fresh_independent_actual804_source_and_finite_proof_passed=True,
        all_new_QQ_memberships_inverses_and_actual_minor_checked_by_Fraction_arithmetic=True,
        whole_pair_rank2_on_a_P_J_open=True, auxiliary_five_rank5_on_a_P_J_via_accepted_binary_review=True,
        raw24_module_surjective_on_FULLU_via_accepted_gcd_binary_finiteSigma_reviews=True,
        FULLU='a*P*F*B*C*J!=0', original_product41_bridge_not_adopted_here=True,
        full57_quartic_or_ambient_exhaustion_claim=False, no_new_GB_search_prime_or_finite_field_promotion=True,
        seconds=time.monotonic()-started)
    note = ROOT/f'research/notes/{STEM}.md'
    paths = [Path(__file__), note, parent/'report.json', archive,
        ROOT/'code/s61_nonsquare_pair_empty_QQ_20261103_v63.py',
        ROOT/'research/notes/S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63.md', finite_code, finite_receipt]+deps
    result['input_sha256'] = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    result['evidence_sha256'] = {'independent-certificate.json': sha(out/'independent-certificate.json'),
        'actual804-independent-replay/report.json': sha(out/'actual804-independent-replay/report.json')}
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: v for k, v in result.items() if k not in ['input_sha256', 'evidence_sha256']}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

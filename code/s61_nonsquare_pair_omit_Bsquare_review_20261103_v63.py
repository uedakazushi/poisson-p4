#!/usr/bin/env python3
"""Audit only the strong2x5 arithmetic theorem; reject unsupported syzygy step."""
import argparse
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import resource
import shutil
import signal
import tarfile
import time
from s61_nonsquare_pair_empty_QQ_review_20261103_v63 import read, sha, decode, add, mul, rem, subtract, power

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_NONSQUARE_PAIR_OMIT_BSQUARE_20261103_v63_REVIEW'


def main(out):
    started = time.monotonic(); resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000)); signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    parent = ROOT/'reports/S61_NONSQUARE_PAIR_OMIT_BSQUARE_20261103_v63'; author = read(parent/'report.json')
    assert author['frozen'] and author['all_literal_QQ_checks_passed']
    archive = ROOT/author['archive']; assert archive.stat().st_size == 123510 and sha(archive) == author['archive_sha256']
    snap = out/'immutable-input'; snap.mkdir()
    with tarfile.open(archive) as tar:
        for member in tar:
            assert member.isfile(); dst = snap/member.name; assert dst.resolve().is_relative_to(snap.resolve())
            dst.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as src, dst.open('wb') as target: shutil.copyfileobj(src, target)
    for rel, h in author['input_sha256'].items(): assert sha(snap/rel if (snap/rel).exists() else ROOT/rel) == h, rel
    old_review = ROOT/'reports/S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63_REVIEW/report.json'
    accepted = read(old_review); assert accepted['all_passed']
    for rel, h in accepted['input_sha256'].items(): assert sha(ROOT/rel) == h, rel
    source = read(snap/'work/nonsquare-v62/pair-finite-proof-v76/certificate.json')
    whole = read(snap/'work/nonsquare-v62/pair-finite-subresultant-v78/whole-pair-cleanup-certificate.json')
    extra = read(snap/'work/nonsquare-v62/quad78-small-gates-v83/omit-B2/certificate.json')
    assert extra['retained_pair_columns'] == [1, 2, 3, 4, 5]
    rows = [[decode(p, 3) for p in row] for row in source['natural_rows']]
    # All three necessary equations and hence the C/Q23 split remain available.
    assert [r['pair_coefficient_index'] for r in source['linear_equation_normalizations']] == [2, 3, 4]
    assert all(r['pair_coefficient_index'] in extra['retained_pair_columns'] for r in source['linear_equation_normalizations'])
    big = whole['degree23_branch']; modulus = decode(big['Q'], 1)
    assert modulus == decode(extra['degree23_Q'], 1)
    Hb, Tb = decode(big['Hbar'], 1), decode(big['Tbar'], 1)
    def evaluate(p):
        images = [Tb, {(1,): Q(1)}, Hb]; powers = []
        for i, image in enumerate(images):
            array = [{(0,): Q(1)}]
            for _ in range(max((ex[i] for ex in p), default=0)): array.append(rem(mul(array[-1], image), modulus))
            powers.append(array)
        result = {}
        for ex, c in p.items():
            term = {(0,): c}
            for array, k in zip(powers, ex): term = rem(mul(term, array[k]), modulus)
            result = add(result, term)
        return rem(result, modulus)
    minor = subtract(mul(rows[0][5], rows[1][1]), mul(rows[1][5], rows[0][1]))
    value = evaluate(minor); assert value == decode(extra['remaining_actual_pair_minor_col1'], 1)
    inverse = decode(extra['remaining_minor_inverse_mod_Q'], 1)
    assert rem(mul(value, inverse), modulus) == {(0,): Q(1)}
    # Special branches use the surviving first-row coefficients BEFORE H is fixed.
    def specialize(p, vv, hh=None):
        result = {}
        for (it, iv, ih), c in p.items():
            ex = (it, ih) if hh is None else (it,)
            coefficient = c*vv**iv*(Q(1) if hh is None else hh**ih)
            result = add(result, {ex: coefficient})
        return result
    branches = []
    for record, (vv, hh, expected_unit) in zip(extra['special_branches'], [(Q(-3), Q(9), Q(3, 625)), (Q(-4, 3), Q(7, 3), Q(41, 20000))]):
        assert Q(record['v']) == vv and Q(record['H']) == hh and Q(record['q0_cross_unit_coefficient']) == expected_unit
        cross = subtract(rows[0][3], {ex: Q(4, 15)*c for ex, c in rows[0][4].items()})
        assert specialize(cross, vv) == {(1, 0): expected_unit}
        assert not specialize(rows[0][5], vv)
        # row1_u2 = -3*(3H(v+4)+5v-12)/32000 fixes this H uniquely.
        assert 3*(vv+4) and 3*hh*(vv+4)+5*vv-12 == 0
        rr = [[specialize(p, vv, hh) for p in row] for row in rows]
        expected_cols = list(itertools.combinations(range(1, 5), 2)); primitive = []
        assert len(record['normalization']) == len(record['raw_minors']) == len(record['primitive_minors']) == 6
        for norm, encoded_raw, encoded_primitive, cols in zip(record['normalization'], record['raw_minors'], record['primitive_minors'], expected_cols):
            assert norm['columns'] == list(cols); i, j = cols
            actual = subtract(mul(rr[0][i], rr[1][j]), mul(rr[0][j], rr[1][i]))
            assert actual == decode(encoded_raw, 1)
            pr = decode(encoded_primitive, 1); pt, pP = norm['unitpowers_t_P']; assert min(pt, pP) >= 0
            factor = mul({(pt,): Q(norm['QQcontent'])}, power({(0,): hh, (1,): Q(-8)}, pP, 1))
            assert mul(pr, factor) == actual; primitive.append(pr)
        cofs = [decode(p, 1) for p in record['Bezout_coefficients']]; assert len(cofs) == 6
        unit = {}
        for co, pr in zip(cofs, primitive): unit = add(unit, mul(co, pr))
        assert unit == {(0,): Q(1)}
        branches.append(dict(v=str(vv), H=str(hh), before_H_fixed_first_row_unit=str(expected_unit)+'*t', six_actual_retained_minors_and_literal_QQ_Bezout_verified=True))
    assert len(extra['special_branches']) == 2
    evidence = dict(accepted_source_and_C_Q23_candidate_proof=str(old_review.relative_to(ROOT)),
        all_three_candidate_minors_use_only_retained_columns=[2, 3, 4, 5],
        actual_col1_u2_minor_rebuilt=True, residue_times_saved_inverse_is1_mod_Q23=True,
        no_irreducibility_or_field_assumption=True, special_branches=branches,
        coefficient_matrix_rank2_everywhere_on_aPJ_after_omitting_Bsquare=True,
        pointwise_to_module='All geometric rank-failure points excluded; maximal-minor ideal is unit on aPJ and local Nakayama gives zero cokernel.',
        conditional_five_linear_syzygy_and_raw78_conclusions_NOT_accepted=True,
        gap='S*Q1 contains A*B*S/6075, so the B*S coefficient also receives a scroll-cofactor contribution; similarly for B*U.',
        existing_product41_uses_constant_five_quadric_criterion_and_is_unaffected=True)
    (out/'independent-certificate.json').write_text(json.dumps(evidence, indent=2)+'\n'); shutil.rmtree(snap)
    result = dict(id=STEM, all_passed_within_explicit_arithmetic_scope=True,
        actual_omitted_Bsquare_2_by5_rank2_on_entire_aPJ_accepted=True,
        independent_Fraction_col1_inverse_and_two_special_sixminor_Bezout_checks=True,
        conditional_five_linear_syzygy_claim_accepted=False, conditional_raw_quartic_rank78_claim_accepted=False,
        raw78_G1_G2_sufficiency_needs_scroll_module_elimination=True,
        source_packet_original_product41_and_prior_CI_open_not_invalidated=True,
        no_GB_new_field_prime_or_heavy_source_replay=True, seconds=time.monotonic()-started)
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md',
        ROOT/'research/notes/S61_NONSQUARE_PAIR_OMIT_BSQUARE_CONDITIONAL_SYZYGY_GAP_20261103_v63.md',
        parent/'report.json', archive, old_review]
    result['input_sha256'] = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    result['evidence_sha256'] = {'independent-certificate.json': sha(out/'independent-certificate.json')}
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n'); signal.alarm(0)
    print(json.dumps({k: v for k, v in result.items() if k not in ['input_sha256', 'evidence_sha256']}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

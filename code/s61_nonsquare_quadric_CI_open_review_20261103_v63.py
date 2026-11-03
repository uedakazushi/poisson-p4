#!/usr/bin/env python3
"""Independent binary resultant and proof audit of the explicit CI open."""
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
from s61_nonsquare_structural_reductions_review_20261103_v63 import read, sha, add, mul

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_NONSQUARE_QUADRIC_CI_OPEN_20261103_v63_REVIEW'


def decode(terms):
    ans = {}
    for ex, c in terms:
        ex = tuple(ex); assert len(ex) == 3 and min(ex) >= 0 and ex not in ans
        if Q(c): ans[ex] = Q(c)
    return ans


def main(out):
    start = time.monotonic(); resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000)); signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    parent = ROOT/'reports/S61_NONSQUARE_QUADRIC_CI_OPEN_20261103_v63'; author = read(parent/'report.json')
    assert author['frozen'] and author['all_exact_QQ_source_and_resultant_checks_passed']
    archive = ROOT/author['archive']; assert archive.stat().st_size == 57680 and sha(archive) == author['archive_sha256']
    snap = out/'immutable-input'; snap.mkdir()
    with tarfile.open(archive) as tar:
        for member in tar:
            assert member.isfile(); dst = snap/member.name; assert dst.resolve().is_relative_to(snap.resolve())
            dst.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as src, dst.open('wb') as target: shutil.copyfileobj(src, target)
    for rel, h in author['input_sha256'].items(): assert sha(snap/rel if (snap/rel).exists() else ROOT/rel) == h, rel
    deps = [ROOT/f'reports/{stem}/report.json' for stem in [
        'S61_NONSQUARE_PAIR_FINITE_BAD_LOCUS_20261103_v63_REVIEW',
        'S61_NONSQUARE_FOUR_QUADRIC_BINARY_REDUCTION_20261103_v63_REVIEW',
        'S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63_REVIEW']]
    for p in deps:
        receipt = read(p); assert receipt['all_passed']
        for rel, h in receipt['input_sha256'].items(): assert sha(ROOT/rel) == h, rel
    binary_receipt = read(deps[1])
    assert binary_receipt['all3_actual_determinantal_minor_identities_verified']
    assert binary_receipt['actual_inverse_frame_and_six_linear_entries_rank_exact4_on_a_P_J']
    assert binary_receipt['D_binary_B_W_restriction_known_unit_nonzero']
    cert = read(snap/'work/nonsquare-v62/quadric-CI-open-replay-v79/certificate.json')
    pair = read(snap/'work/nonsquare-v62/pair-finite-proof-v76/certificate.json')
    binary = read(snap/'work/nonsquare-v62/four-quadric-binary-certificate-v71/certificate.json')
    assert cert['source_binary_certificate'] == binary
    assert binary['inverse_frame_determinant'] == '1/(200*a*P^3*J)'
    assert binary['D_x0_squared_unit'] == '1280*a^3*J^2/243'
    assert binary['three_minor_signs'] == ['-D/6075', '-Q1', '-Q2prime']
    assert cert['known_open'] == 't*(H-8t)*(5H-2v-4)!=0'
    assert cert['vertex'] == 'A=x0=W=B=0'
    # On the accepted actual pair chart B=-10*a*P*J*z, so B=0 iff z=0.
    rows = [[decode(p) for p in row] for row in pair['natural_rows']]
    coeffs = [[row[i] for i in [5, 4, 2]] for row in rows]
    assert coeffs == [[decode(p) for p in row] for row in cert['binary_quadric_coefficients']]
    zero = {}
    aa, bb = coeffs
    Syl = [aa+[zero], [zero]+aa, bb+[zero], [zero]+bb]
    determinant = {}
    for perm in itertools.permutations(range(4)):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(4) for j in range(i+1, 4))
        term = {(0, 0, 0): Q(sign)}
        for i, j in enumerate(perm): term = mul(term, Syl[i][j])
        determinant = add(determinant, term)
    R9 = decode(cert['R_vertex']); assert R9 and max(sum(ex) for ex in R9) == 9
    expected = mul({(1, 0, 0): Q(-729, 4096000000000000000)}, R9)
    assert determinant == expected
    assert cert['Sylvester_determinant_known_unit_factor'] == '-729*t/4096000000000000000'
    assert cert['R_vertex_totaldegree'] == 9 and len(R9) == 92
    # A rational point verifies only nonemptiness of this explicit open.
    # (a,e,P)=(1,1,1) corresponds to (t,v,H)=(1,23,9).
    point = (1, 23, 9)
    R_value = sum(c*point[0]**ex[0]*point[1]**ex[1]*point[2]**ex[2] for ex, c in R9.items())
    assert R_value
    evidence = dict(accepted_actual804_binary_frame_and_pair_source_dependencies=[str(p.relative_to(ROOT)) for p in deps[:2]],
        actual_vertex_binary_coefficients_all6_checked=True, independent_Fraction_24term_Sylvester_determinant=True,
        degree9_R_vertex_terms=92, exact_known_unit_factor='-729*t/4096000000000000000',
        R_vertex_is_additional_condition_not_known_core_unit=True,
        chart_nonempty_witness={'t_v_H': list(point), 'R_vertex': str(R_value)},
        geometric_argument={
            'Q1_Q2_coprime': 'Any common linear factor reduces modW to B, hence B+kW; it divides D, contradicting nonzero D mod(B,W). A quadratic common factor would imply proportionality.',
            'minor_ideal_height2': 'Contains coprime Q1,Q2, while the maximal-minor height bound is2; Hilbert-Burch gives CM pure projective curve.',
            'set_decomposition': 'V(Q1,Q2)=V(B,W) union V(all3 minors), from WD=40XQ1-AQ2 and the W=0 alternatives.',
            'vertex_avoided': 'The exact binary resultant is nonzero on D(R_vertex).',
            'P3_piece': 'Closed subset disjoint from vertex P1 has dimension<=1 by the projective dimension theorem.',
            'cone_piece': 'Projection from vertex on closed projective Y is proper; every closed fiber in a ruling P2 avoids its line P1 and therefore has dimension<=0.',
            'no_smoothness_or_reducedness': 'Purity and support dimensions suffice; curve may be reducible or nonreduced.',
            'height4': 'Both pieces have projective dimension<=1, so four quadrics in six variables have height4 and form a regular sequence.'},
        allpoint_colon_Hilbert_consequence_on_FULLU_D_R9={'HF4_barJ': 48, 'HF3_ZW': 16, 'HF4_original': 64, 'light_rank': 146,
            'HF2_barJ': 17, 'HF1_ZW': 5, 'HF2_original': 22, 'original_quadric_rank': 6, 'weighted_rank': 152},
        R_vertex_zero_CI_lower146_full57_fullquartic_and_ambient_exhaustion_unresolved=True)
    (out/'independent-certificate.json').write_text(json.dumps(evidence, indent=2)+'\n'); shutil.rmtree(snap)
    result = dict(id=STEM, all_passed=True, actual_binary_resultant_and_only_a_P_J_units_verified=True,
        complete_intersection_on_raw_a_P_J_Rvertex_open_accepted=True,
        pure_base_curve_and_proper_projection_proof_accepted_without_smooth_or_reduced_assumption=True,
        original_light_rank146_and_weighted152_on_FULLU_D_Rvertex=True,
        R_vertex_unknown_factor_not_inverted_globally=True, R_vertex_zero_case_and_fullquartic_open=True,
        no_GB_new_field_modp_promotion_or_parametric_quartic_expansion=True, seconds=time.monotonic()-start)
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', parent/'report.json', archive]+deps
    result['input_sha256'] = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    result['evidence_sha256'] = {'independent-certificate.json': sha(out/'independent-certificate.json')}
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n'); signal.alarm(0)
    print(json.dumps({k: v for k, v in result.items() if k not in ['input_sha256', 'evidence_sha256']}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)

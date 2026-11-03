#!/usr/bin/env python3
"""Compose independently accepted raw finiteness and the physical product bridge."""
import argparse
import hashlib
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_NONSQUARE_PRODUCT41_FINITE_BAD_LOCUS_20261103_v63'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output):
    start = time.monotonic()
    paths = [
        ROOT/'reports/S61_NONSQUARE_PAIR_FINITE_BAD_LOCUS_20261103_v63_REVIEW/report.json',
        ROOT/'reports/S61_NONSQUARE_RAW24_PRODUCT41_BRIDGE_20261103_v63_REVIEW/report.json',
    ]
    finite, bridge = [json.loads(p.read_text()) for p in paths]
    assert finite['all_passed'] and bridge['all_passed']
    assert finite['actual_raw24_defect_dimension_at_most_zero_on_FULLU_using_prior_reviews']
    assert bridge['pointwise_RAW24_rank24_implies_original_product_rank41']
    assert not bridge['original_product_rank41_outside_finite_FULLU_subset_accepted']
    checked = {}
    for path in paths:
        report = json.loads(path.read_text())
        for key in ['input_sha256', 'evidence_sha256']:
            for name, expected in report.get(key, {}).items():
                found = next((p for p in [ROOT/name, path.parent/name]
                              if p.is_file() and sha(p) == expected), None)
                assert found is not None, (path, name)
                checked[str(found.relative_to(ROOT))] = expected
    output.mkdir(parents=True, exist_ok=False)
    note = ROOT/'research/notes'/(ID+'.md')
    result = {
        'id': ID, 'all_passed': True,
        'field': 'QQ, every geometric point on FULLU=a*P*F*B0*C*J!=0',
        'actual_original_six_quadric_times_seven_light_rank41_outside_finite_subset': True,
        'original_product_rank_failure_dimension_at_most_zero': True,
        'proof': 'Original product defect is contained in the independently proved finite actual RAW24 defect.',
        'whole_open_rank41_or_empty_bad_locus_not_claimed': True,
        'full57_quartic_Jzero_ambient_removal_not_claimed': True,
        'constructive_pivot_count_unchanged': True,
        'no_new_CAS_or_field_computation': True,
        'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), note, *paths]},
        'checked_dependency_sha256': checked,
        'seconds': time.monotonic()-start,
    }
    (output/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'all_passed': True, 'seconds': result['seconds']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())

#!/usr/bin/env python3
"""Bounded proof/source audit of the good-prime seventh-jet upper bound."""
import argparse
import hashlib
import json
from math import comb
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]


def verify():
    start = time.time()
    primary = 'reports/S61_MINIMAL_FREE_FLAG_SEVENTH_UPPER_BOUND_20261101.json'
    data = json.loads((ROOT/primary).read_text())
    assert data['all_passed']
    for path, digest in data['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
    assert data['good_prime'] == 32003 and data['specialization'] == [1,2,3,4,5]
    assert comb(18,7) == data['truncated_ring_columns'] == 31824
    assert 31824-data['exact_modular_ideal_rank'] == 4335
    assert sum(data['generic_hilbert_function_through_six']) == 2568
    assert 4335-2568 == data['generic_seventh_hilbert_upper_bound'] == 1767
    assert data['known_seventh_quotient_before_new_relations']-1767 == 30
    assert data['generic_new_seventh_initial_rank_lower_bound'] == 30
    prime = data['good_prime']
    # Independent arithmetic check of the stated native/numpy safety bounds.
    assert 45*(prime-1)**2 < 2**63
    maximum_pairs = max(comb(10+d,d)*comb(17-d,7-d) for d in range(1,4))
    assert maximum_pairs*(prime-1) < 2**53
    paths = list(data['input_sha256'])+[primary,
        'research/notes/S61_MINIMAL_FREE_FLAG_SEVENTH_UPPER_BOUND_20261101.md',
        str(Path(__file__).relative_to(ROOT)),
        'research/notes/S61_MINIMAL_FREE_FLAG_SEVENTH_UPPER_BOUND_REVIEW_20261101.md']
    return dict(id='S61_MINIMAL_FREE_FLAG_SEVENTH_UPPER_BOUND_REVIEW_20261101',
                all_passed=True,proof_and_source_review_accepted=True,
                primary_input_hashes_verified=True,independent_arithmetic_checked=True,
                independent_raw7_or_closure_reimplementation=False,
                primary_producer_replayed=False,
                accepted_generic_seventh_upper_bound=1767,
                accepted_new_seventh_directions_lower_bound=30,
                limitations=['This review checks the finite-matrix proof and source, not a second raw7 contraction or closure implementation.',
                             'The known homogeneous ideal HF7=1797 is a primary computational dependency.',
                             'No generic HF7 equality, full initial ideal or ambient t-nilpotence follows.'],
                input_sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
                seconds=time.time()-start)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print('FREE FLAG SEVENTH UPPER REVIEW PASS',result['seconds'])

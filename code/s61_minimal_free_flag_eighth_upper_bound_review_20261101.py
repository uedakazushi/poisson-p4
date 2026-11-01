#!/usr/bin/env python3
"""Proof/source audit of a finite-row-space eighth-order upper bound."""
import argparse
import hashlib
import json
from math import comb
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]


def verify():
    start=time.time()
    primary='reports/S61_MINIMAL_FREE_FLAG_EIGHTH_UPPER_BOUND_20261101.json'
    data=json.loads((ROOT/primary).read_text())
    assert data['all_passed']
    for path,digest in data['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    assert comb(19,8)==data['truncated_ring_columns']==75582
    assert 75582-data['exact_modular_ideal_rank']==6745
    assert sum(data['generic_hilbert_function_through_six'])==2568
    assert 6745-2568==data['generic_sum_hilbert_seventh_eighth_upper_bound']==4177
    assert data['generic_quotient_by_m9_length_upper_bound']==6745
    prime=data['good_prime'];assert prime==32003
    assert 45*(prime-1)**2<2**63
    assert max(comb(10+d,d)*comb(18-d,8-d) for d in range(1,5))*(prime-1)<2**53
    paths=list(data['input_sha256'])+[primary,
        'research/notes/S61_MINIMAL_FREE_FLAG_EIGHTH_UPPER_BOUND_20261101.md',
        str(Path(__file__).relative_to(ROOT)),
        'research/notes/S61_MINIMAL_FREE_FLAG_EIGHTH_UPPER_BOUND_REVIEW_20261101.md']
    return dict(id='S61_MINIMAL_FREE_FLAG_EIGHTH_UPPER_BOUND_REVIEW_20261101',all_passed=True,
        proof_and_source_review_accepted=True,primary_input_hashes_verified=True,
        finite_row_membership_suffices_without_complete_standard_basis=True,
        accepted_generic_length_m9_upper_bound=6745,accepted_generic_hf7_plus_hf8_upper_bound=4177,
        independent_raw8_or_basis_reimplementation=False,primary_producer_replayed=False,
        limitations=['The review uses the certified relation count only as a rank lower bound; it does not independently certify the exact modular ideal rank.',
                     'No separate generic HF8 bound 2410 or any higher formal-ring conclusion is inferred.'],
        input_sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print('FREE FLAG EIGHTH UPPER REVIEW PASS',result['seconds'])

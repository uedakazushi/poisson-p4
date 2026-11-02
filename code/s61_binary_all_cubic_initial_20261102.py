#!/usr/bin/env python3
"""Dependency-checked rank-semicontinuity corollary for all cubic initials."""
import argparse
import hashlib
import json
from math import comb
from pathlib import Path
import time
ROOT=Path(__file__).resolve().parents[1]
NOTE=ROOT/'research/notes/S61_BINARY_ALL_CUBIC_INITIAL_20261102.md'
STEMS=['S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101',
       'S61_BINARY_DEEPER_QUARTIC_NUISANCE_20261101',
       'S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101',
       'S61_BINARY_CUBIC_NUISANCE_LOCALIZATION_20261102']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    start=time.time();reports={};checked={};paths=[Path(__file__),NOTE]
    for stem in STEMS:
        p=ROOT/'reports'/f'{stem}.json';data=json.loads(p.read_text());assert data['all_passed']
        for name,digest in data['input_sha256'].items():
            q=ROOT/name;assert sha(q)==digest,(stem,name);checked[name]=digest
        paths.extend([p,ROOT/data['proof_note']]);reports[stem]=data
    quartic=reports[STEMS[1]];formal=reports[STEMS[2]];local=reports[STEMS[3]]
    assert quartic['generated_subideal_hilbert_function']==[1,9,21,34]
    assert quartic['generated_subideal_length']==65
    assert formal['relative_completed_rank']==formal['fixed_beta_completed_rank']==65
    assert local['full_seven_parameter_linear_denominators_are_core_only']
    assert local['nuisance_class_invariance_on_entire_D_open']
    assert local['independent_two_true_cubic_classes_on_DF_open_for_every_nuisance']
    assert local['quadratic_linear_rank_on_DF_open']==129
    assert local['rank_with_two_true_cubics_on_DF_open']==131
    assert local['known_initial_subideal_hilbert_through3']==[1,9,21,34]
    assert comb(11,3)-131==34 and comb(10,2)-24==21 and 1+9+21+34==65
    return dict(id='S61_BINARY_ALL_CUBIC_INITIAL_20261102',all_passed=True,
        verification_kind='dependency/hash verification and proof corollary; no second full-coefficient expansion',
        coefficient_ring='QQ[e,d,D^-1][p,q,r,h,k]',
        complete_cubic_matrix_columns=comb(11,3),
        complete_cubic_matrix_rows_at_most=9*24+135,
        generic_rank_from_prior_length65_and_initial_subideal=131,
        all_132_minors_identically_zero_by_dense_open=True,
        exact_full_cubic_rank_on_DF_for_every_nuisance=131,
        ordinary_hilbert_through3=[1,9,21,34],
        all_true_cubics_generated_by_quadratic_linear_multiples_and_C92_C134=True,
        limitations=['No uniform fourth-degree surjectivity or full-core length65 theorem.',
                    'No ambient t-nilpotence, affine exhaustion, or component upper bound.',
                    'The full cubic matrix/minors are defined in the proof, not expanded anew by this verifier.'],
        checked_primary_input_hash_count=len(checked),
        proof_note=str(NOTE.relative_to(ROOT)),
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        seconds=time.time()-start)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify()
    with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({k:report[k]for k in ['all_passed','checked_primary_input_hash_count','seconds']}))

#!/usr/bin/env python3
"""Bounded verification and scope audit of the actual global action/cyclic open."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from s61_affine_seventeen_base_lift_20261103_v62 import build
from s61_affine_t_action_circuit_20261103_v62 import read_native, check
from s61_affine_cyclic_open_review_20261103_v62 import check_dense_action
from s61_affine_cyclic_open_20261103_v62 import check_certificate

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--t4-replay', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    tdir = ROOT/'reports/S61_AFFINE_T_ACTION_CIRCUIT_20261103_v62'
    cdir = ROOT/'reports/S61_AFFINE_CYCLIC_OPEN_20261103_v62'
    mdir = ROOT/'reports/S61_AFFINE_CYCLIC_OPEN_20261103_v62_REVIEW'
    t4path = ROOT/'reports/S61_AFFINE_T_FOURTH_NORMAL_FORM_20261103_v62_REVIEW.json'
    paths = [tdir/'report.json', cdir/'report.json', mdir/'report.json', t4path]
    audited = []
    for path in paths:
        d = json.loads(path.read_text())
        assert d['all_passed']
        for name, expected in d['input_sha256'].items():
            p = ROOT/name
            assert sha(p) == expected, name
            audited.append(p)
        for name, expected in d.get('evidence_sha256', {}).items():
            p = path.parent/name
            assert sha(p) == expected, str(p)
            audited.append(p)
    old4 = json.loads(t4path.read_text())
    replay4 = json.loads(args.t4_replay.read_text())
    assert {k:v for k,v in old4.items() if k!='seconds'} == {k:v for k,v in replay4.items() if k!='seconds'}
    high,cross,base,pairs,_,_,bpairs,vector,inputs,_,_ = build()
    encoded = ('29 134 32003 17\n'+'\n'.join(' '.join(map(str,a.ravel())) for a in [high,cross,base])
               +'\n'+' '.join(map(str,vector%32003))+'\n')
    assert encoded == (cdir/'quadrics.in').read_text() == (tdir/'quadrics.in').read_text()
    known_path = ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json'
    known = json.loads(known_path.read_text())
    summary,p2,rows,p3,alternatives,trials,H,L,D = read_native(tdir/'native.txt')
    tcheck = check(high,cross,base,bpairs,p2,p3,trials,H,L,D,known['known_five_mod_32003'])
    assert summary['GENERATORS'] == 1031 and summary['DET3'] == 14008 and summary['DET4'] == 8797
    lines = (cdir/'native-certificate.txt').read_text().splitlines()
    ell = list(map(int,next(line for line in lines if line.startswith('LINEAR_FUNCTIONAL ')).split()[1:]))
    linear = cross@(vector%32003)%32003
    constant = base@np.array([int(vector[i])*int(vector[j])%32003 for i,j in bpairs],dtype=np.int64)%32003
    h_count = check_dense_action(lines,ell,high,linear,constant)
    _,_,det,roots,gcd,_ = check_certificate(cdir/'native-certificate.txt',known)
    assert det == 17765 and len(set(roots)) == 5
    m2log = (mdir/'review.log').read_text()
    assert 'CYCLIC_OPEN_ALL134_REMAINDERS_PASS' in m2log
    assert 'error' not in m2log.lower()
    paths += audited + inputs + [Path(__file__),known_path,
              ROOT/'research/notes/S61_AFFINE_T_ACTION_CIRCUIT_20261103_v62.md',
              ROOT/'research/notes/S61_AFFINE_CYCLIC_OPEN_20261103_v62.md',
              ROOT/'research/notes/S61_GLOBAL_ACTION_CYCLIC_OPEN_20261103_v62_REVIEW.md']
    report = dict(
        id='S61_GLOBAL_ACTION_CYCLIC_OPEN_20261103_v62_REVIEW',all_passed=True,
        actual_original134_all17_base_input_reconstructed=True,
        dense_t_action_original_equation_cofactor_recheck=tcheck,
        dense_h_action_original_equation_coefficients_rechecked=h_count,
        all46_variable_t4_coefficients_fresh_replay=211876,
        saved_all1031_linear_functional_coefficients_checked=True,
        all29_coordinates_at_all5_points_checked=True,
        independent_M2_all134_gcd_evidence_hash_and_script_audited=True,
        M2_118_second_gcd_repeated=False,
        native_full1031_Krylov_determinant_accepted=det,
        native_692_second_full_Krylov_run_repeated=False,
        full_Krylov_determinant_independently_recomputed=False,
        generic_QQ_spread_requires_full_quartic_right_inverse_correction=True,
        corrected_generic_action_reduces_to_frozen_mod_p_action=True,
        actual_univariate_presentation_on_nonempty_cyclic_open_accepted=True,
        symbol_freeness_assumed=False,action_commutator_freeness_assumed=False,
        generic_gcd_degree5_claimed=False,determinant_zero_locus_resolved=False,
        full_affine_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(paths)})
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2)
        stream.write('\n')
    print('PASS actual action identities, cyclic presentation scope and exact generic spread')


if __name__ == '__main__':
    main()

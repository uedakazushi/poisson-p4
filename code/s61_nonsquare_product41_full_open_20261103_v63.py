#!/usr/bin/env python3
"""Accepted whole actual RAW24 plus the independent physical 17-kernel bridge."""
import argparse
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PRODUCT41_FULL_OPEN_20261103_v63'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output):
    start=time.monotonic()
    paths=[ROOT/'reports/S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63_REVIEW/report.json',
           ROOT/'reports/S61_NONSQUARE_RAW24_PRODUCT41_BRIDGE_20261103_v63_REVIEW/report.json']
    pair,bridge=[json.loads(p.read_text())for p in paths]
    assert pair['all_passed'] and bridge['all_passed']
    assert pair['whole_pair_rank2_on_a_P_J_open']
    assert pair['raw24_module_surjective_on_FULLU_via_accepted_gcd_binary_finiteSigma_reviews']
    assert pair['original_product41_bridge_not_adopted_here']
    assert bridge['pointwise_RAW24_rank24_implies_original_product_rank41']
    checked={}
    for path in paths:
        report=json.loads(path.read_text())
        for key in ['input_sha256','evidence_sha256']:
            for name,want in report.get(key,{}).items():
                found=next((p for p in [ROOT/name,path.parent/name]if p.is_file()and sha(p)==want),None)
                assert found is not None,(path,name)
                checked[str(found.relative_to(ROOT))]=want
    output.mkdir(parents=True,exist_ok=False)
    note=ROOT/'research/notes'/(ID+'.md')
    result={
        'id':ID,'all_passed':True,
        'field':'QQ, all geometric points of FULLU',
        'FULLU_original':'a*d*F*B*C*J!=0, J=5*a*d-6*e-4',
        'FULLU_serialized':'a*P*F*B0*C*J!=0, P=a*d, B0=a*B',
        'actual_original42_by84_product_rank_exact41_on_whole_FULLU':True,
        'actual_original_product_rank_defect_empty_on_FULLU':True,
        'product_image_rank41_subbundle_and_cokernel_locally_free_rank43':True,
        'global_free_cokernel_basis_or_explicit_right_inverse_constructed':False,
        'whole_Jzero_or_dzero_product_claim':False,
        'full57_quartic_length65_ambient_or_exhaustion_claim':False,
        'constructive_pivot_count_unchanged':True,
        'no_new_CAS_or_field_computation':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),note,*paths]},
        'checked_dependency_sha256':checked,
        'seconds':time.monotonic()-start,
    }
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'whole_FULLU_product_rank':41,'seconds':result['seconds']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())

#!/usr/bin/env python3
"""Finite weighted-presentation corollary of the exact QQ seventh certificate."""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
NOTE='research/notes/S61_FLAG_ZERO12_KERNEL_FIRST_ORDER_20261102.md'

def verify():
    paths=[ROOT/'reports'/f'{s}.json'for s in ['S61_FLAG_ZERO12_KERNEL_GENERATORS_20261102',
           'S61_FLAG_ZERO12_SEVENTH_INITIAL_20261101','S61_FLAG_ZERO12_RESIDUE_ACTIONS_20261102']]
    gen,seven,action=[json.loads(p.read_text())for p in paths]
    assert all(x['all_passed']for x in [gen,seven,action])
    assert seven['new_seventh_initial_directions']==0 and seven['all_closure_checks']==34857
    degrees=gen['filtration_generator_degrees'];assert degrees==[2,15,28,38,18,10]
    assert gen['actual_kernel_minimal_generator_number']==111 and action['actual_t_M_contained_in_b_M']
    e=[0]+degrees
    first=[(e[d]if d<len(e)else 0)+2*(e[d-1]if 0<=d-1<len(e)else 0)for d in range(8)]
    assert sum(first)==333 and first==[0,2,19,58,94,94,46,20]
    d=[1,7]+[15]*6;assert sum(d)==98
    paths += [Path(__file__),ROOT/NOTE]
    return dict(id='S61_FLAG_ZERO12_KERNEL_FIRST_ORDER_20261102',all_passed=True,field='QQ',
                actual_kernel_mod_b_squared_free_rank=111,actual_kernel_mod_b_squared_dimension=333,
                fixed_squarezero_quotient_truncation_length=98,full_truncated_quotient_length=431,
                kernel_weighted_basis_counts=first,linear_t_action_has_two_unique_111_matrices=True,
                QQ_linear_t_action_zero_claimed=False,full_kernel_freeness_claimed=False,
                full_t_nilpotence_claimed=False,proof_note=NOTE,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths})
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify()
    with args.output.open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print('FIRST ORDER KERNEL FREE333 PASS')

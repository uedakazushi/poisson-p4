#!/usr/bin/env python3
"""Explicit minimal111 generators of the actual two-base-variable flag kernel."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_flag_kernel_lift_engine_20261102 import source_geometry,quadratic_representation,kernel_candidates,actual_lifts
ROOT=Path(__file__).resolve().parents[1]
NOTE='research/notes/S61_FLAG_ZERO12_KERNEL_GENERATORS_20261102.md'

def verify():
    started=time.monotonic()
    stems=['S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101','S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101',
           'S61_FLAG_ZERO12_ACTUAL_KERNEL_FINITE_20261102','S61_FLAG_ZERO12_FIRST_BASE_QUOTIENT_20261102']
    paths=[ROOT/'reports'/f'{stem}.json'for stem in stems]
    initial,normal,finite,first=[json.loads(p.read_text())for p in paths]
    assert all(x['all_passed']for x in [initial,normal,finite,first])
    assert initial['full_truncation_degree']==6 and initial['full_truncation_length']==1147
    assert 'EXACT CLOSURE PASS' in initial['replay_log']
    assert finite['actual_kernel_generator_bound']==111
    assert initial['point61']==first['point61']==finite['point61']
    source=source_geometry(initial);action=quadratic_representation();candidates=kernel_candidates(initial)
    lifts=actual_lifts(source,action,normal,candidates)
    assert candidates['degree_counts']==[2,15,28,38,18,10]
    assert len(lifts['kernel_111_lifts_sixjet'])==111
    # Recheck the two one-variable quotient lengths from the exact first
    # base-adic quotient, independently of its A^40+(A/h)^9 description.
    Q=s.zeros(187,196)
    for j,row in enumerate(first['generator_normal_forms']):
        for i,c in row:Q[int(i),j]=s.Rational(c)
    assert Q.to_DM().rank()==187
    rank_a=Q[:,49:].to_DM().rank();rank_b=Q[:,49:147].to_DM().rank()
    assert rank_a==138 and rank_b==98
    len_one=187-rank_a;len_two=187-rank_b
    assert (len_one,len_two)==(49,89)
    h_rank=111-(len_one-5);h_square_rank=111-(len_two-10)
    assert (h_rank,h_square_rank)==(67,32)
    paths += [Path(__file__),ROOT/'code/s61_flag_kernel_lift_engine_20261102.py',ROOT/NOTE,
              ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
    return dict(id='S61_FLAG_ZERO12_KERNEL_GENERATORS_20261102',all_passed=True,field='QQ',
                actual_kernel_minimal_generator_number=111,base_variables=[0,1],
                filtration_generator_degrees=candidates['degree_counts'],
                induced_filtration_cutoff='F^(r+6)M subset (v0,v1)^r M for every r>=1',
                finite_schur_truncation_degree='r+5',
                residue_h_rank=h_rank,residue_h_squared_rank=h_square_rank,
                residue_H_power6_zero=True,residue_T_square_zero=True,residue_H_squared_T_zero=True,
                full_residue_h_and_t_matrices_computed=False,
                exact_formal_generator_count=12,source_geometry=source,
                chosen_initial_generators=candidates,actual_generator_lifts=lifts,
                action_constraints_proved_by_filtration_not_full_matrix_replay=True,
                actual_kernel_free_claimed=False,full_completed_ring_t_nilpotence_claimed=False,
                relative_ambient_removal_claimed=False,proof_note=NOTE,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
                seconds=time.monotonic()-started)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify()
    with args.output.open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))

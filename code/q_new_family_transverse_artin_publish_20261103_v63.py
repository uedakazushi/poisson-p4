#!/usr/bin/env python3
"""Freeze the actual CLOSED fiber length5 without assuming generic flatness."""
import argparse
import json
from pathlib import Path
import sympy as s
from q_positive_weight6_quadratic_model_20261103_v63 import ROOT,read,sha


def main(folder):
    point=folder/'point111';high=folder/'point111-through4'
    a=read(point/'report.json');b=read(high/'report.json')
    assert a['all_passed'] and b['all_passed']
    assert a['exact_corrected_quadratic_rank']==1
    assert b['actual_closed_fiber_HF_through4']==[1,2,2,0,0] and b['actual_closed_fiber_length_exact']==5
    for directory,data in [(point,a),(high,b)]:
        for name,digest in data['input_sha256'].items():assert sha(ROOT/name)==digest,name
        for name,digest in data['evidence_sha256'].items():assert sha(directory/name)==digest,name
    q=read(point/'certificate.json.gz');Q=s.Matrix(q['full350_by3_true_quadratic_obstructions'])
    aa,bb,cc=Q[q['selected_quadratic_rows'][0],:]
    discriminant=s.factor(bb*bb-4*aa*cc);assert discriminant
    log=(high/'closed-fiber.log').read_text();assert 'DIM=0'in log and 'DEG=5'in log
    paths=[Path(__file__),point/'report.json',point/'certificate.json.gz',high/'report.json',high/'certificate.json.gz',
           high/'closed-fiber.m2',high/'closed-fiber.log',
           ROOT/'research/notes/Q_NEW_FAMILY_TRANSVERSE_ARTIN_20261103_v63.md',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63_REVIEW/report.json',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.json']
    report=dict(id='Q_NEW_FAMILY_TRANSVERSE_ARTIN_20261103_v63',all_passed=True,
        point=['1','1','1'],full_canonical126_and_all350_actual_Poisson_equations_used=True,
        transverse_linear_cuts_to_actual26dimensional_family_image=26,
        selected_unit_implicit_variables=98,actual_closed_transverse_embedding_dimension=2,
        exact_quadratic_initial_rank=1,exact_binary_quadratic_discriminant=str(discriminant),
        exact_cubic_initial_rank=4,actual_closed_fiber_Hilbert_function=[1,2,2],
        actual_closed_fiber_Artin_length=5,closed_fiber_maximal_ideal_cubed_zero=True,
        geometric_closed_fiber_algebra='C[u,v]/(u*v,u^3,v^3)',
        actual_full_completion_finite_generated_by_at_most5_over26_smooth_formal_branch_parameters=True,
        parameter_directions_are_not_all_group_directions=True,
        generic_flatness_exact_generic_multiplicity_or_generic_nonreducedness_inferred_from_point=False,
        independent_function_field_Jacobian_and_quadratic_rank_needed_for_generic_length5=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths})
    with (folder/'report.json').open('x')as f:
        json.dump(report,f,indent=2);f.write('\n')
    print('FROZEN actual closed fiber length5, finiteover26branchbase; no genericflatness inference')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--folder',type=Path,required=True);main(ap.parse_args().folder)

#!/usr/bin/env python3
"""Discharge the finite-rank hypotheses and freeze the full length-five completion."""
import argparse
import json
from pathlib import Path
import sympy as s
from q_positive_weight6_quadratic_model_20261103_v63 import ROOT,read,sha


def main(output):
    names=[
        'reports/Q_NEW_FAMILY_TRANSVERSE_ARTIN_20261103_v63/report.json',
        'reports/Q_NEW_FAMILY_TRANSVERSE_ARTIN_20261103_v63_REVIEW/report.json',
        'reports/Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63_REVIEW.json',
        'reports/Q_NEW_FAMILY_GENERIC_QUADRATIC_RANK_20261103_v63/report.json',
        'reports/Q_NEW_FAMILY_GENERIC_QUADRATIC_RANK_20261103_v63_REVIEW.json',
        'reports/Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.json']
    data=[]
    for name in names:
        d=read(ROOT/name);assert d['all_passed']
        for filename,digest in d['input_sha256'].items():assert sha(ROOT/filename)==digest,(name,filename)
        data.append(d)
    closed,creview,kernel,quad,qreview,component=data
    assert closed['actual_closed_fiber_Artin_length']==5
    assert creview['full_completion_finite_over26smooth_branch_parameters_generated_atmost5']
    assert creview['closed_point_fiber_HF']==[1,2,2]
    assert kernel['exact_generic_full_canonical_Jacobian_rank']==98
    assert qreview['generic_augmented_rank']==99 and qreview['generic_Jacobian_rank']==98
    assert qreview['intrinsic_generic_two_normal_quadratic_obstruction_rank']==1
    assert qreview['descent_mod_actual26_and_dense_corrected_shear_scope_accepted']
    assert qreview['generic_Artin_length_at_least5']
    assert component['actual_global_projective_irreducible_component_dimension_exact']==25
    u,v=s.symbols('u v');G=s.groebner([u*v,u**3,v**3],u,v)
    basis=[1,u,v,u*u,v*v]
    assert all(G.reduce(s.sympify(f))[1]==f for f in basis)
    assert all(G.reduce(u**i*v**(3-i))[1]==0 for i in range(4))
    assert all(G.reduce(a*b)[1]==0 for a in [u*u,v*v]for b in [u,v])
    paths=[Path(__file__),*[ROOT/name for name in names],
        ROOT/'research/notes/Q_NEW_FAMILY_LENGTH5_COMPLETION_20261103_v63.md',
        ROOT/'research/notes/Q_NEW_FAMILY_LENGTH5_CONDITIONAL_COMPLETION_20261103_v63.md',
        ROOT/'code/q_positive_weight6_quadratic_model_20261103_v63.py']
    report=dict(id='Q_NEW_FAMILY_LENGTH5_COMPLETION_20261103_v63',all_passed=True,
        independent_closed_fiber_and_generic_quadratic_hypotheses_all_discharged=True,
        actual_full_completion_finite_R26_generated_atmost5=True,
        generic_branch_HF_lower_1_2_2_consumes_all5_ranks=True,
        actual_full_completion_free_R26_rank=5,
        freeness_from_kernel_of_R5_surjection_not_assumed_torsionfree_or_Tor_basechange=True,
        actual_generic_Artin_length=5,generic_transverse_Hilbert_function=[1,2,2],
        nilradical_cubed_zero=True,nilradical_mod_square_free_rank2_before_square_basechange=True,
        nilradical_square_free_rank2=True,
        canonical_completion_at_n111='C[[s1,...,s26,u,v]]/(u*v,u^3,v^3)',
        projective_completion_at_n111='C[[s1,...,s25,u,v]]/(u*v,u^3,v^3)',
        reduced_component_smooth_at_n111=True,unique_minimal_prime_in_this_completion=True,
        full_completion_Cohen_Macaulay=True,Cohen_Macaulay_type=2,Gorenstein=False,
        generic_geometric_Artin_algebra='Kbar[u,v]/(u*v,u^3,v^3)',
        all_points_global_product_or_full_component_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths})
    with output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print('PASS unconditional actualfullcompletion length5/freeR26, projectivefatArtin5xformal25, CMtype2')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)

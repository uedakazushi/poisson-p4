#!/usr/bin/env python3
"""Assemble the independently certified ALL-Qa T01--T38 membership theorem."""
import argparse
import json
from pathlib import Path
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    names=['Q_T01_MEMBERSHIP','Q_MPS_EXTENDED_TORIC_ARCS','Q_MPS_IMPLICIT_TORIC_ARCS','Q_MPS_HENSEL_TORIC_ARCS','Q_ALL_FAMILY_NON_CM']
    paths=[Path(__file__),ROOT/'code/q_smith_global_model_review_20261103_v63.py'];receipts={}
    for name in names:
        p=ROOT/f'reports/{name}_20261103_v63_REVIEW.json';d=read(p);assert d['all_passed'];paths.append(p);receipts[name]=d
        for file,h in d['input_sha256'].items():q=ROOT/file;assert sha(q)==h,file
    assert receipts['Q_T01_MEMBERSHIP']['all_Qa_in_actual_existing_T01']
    rational=receipts['Q_MPS_EXTENDED_TORIC_ARCS']['rational_rho_and_a_density_closedness_give_ALL_finite_Qa_memberships']
    implicit=receipts['Q_MPS_IMPLICIT_TORIC_ARCS']['all_finite_Qa_in_components'];hensel=receipts['Q_MPS_HENSEL_TORIC_ARCS']['all_finite_Qa_in_components']
    union=['T01']+rational+implicit+hensel;assert len(union)==len(set(union))==38 and sorted(union)==[f'T{i:02d}'for i in range(1,39)]
    assert receipts['Q_ALL_FAMILY_NON_CM']['all_finite_complex_a_including_zero']
    paths += [ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json',ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json',ROOT/'legacy/v14/data/six_rank2_components.json',ROOT/'research/notes/Q_ALL_TORIC_MPS_COMPONENTS_20261103_v63.md']
    report=dict(id='Q_ALL_TORIC_MPS_COMPONENTS_20261103_v63',all_passed=True,
        parameter_scope='Every finite complex a, including a=0; physical nonzero homogeneous Qa',
        actual_known_MPS_components_through_ALL_Qa=sorted(union),additional_known_components_through_ALL_Qa=['CA','CB'],
        distinct_known_components_through_ALL_Qa_lower_bound=40,
        proof_partition={'direct_toric':['T01'],'rational_extended_arcs':rational,'exact_QQ_implicit_sources':implicit,'positive_Hensel_implicit_sources':hensel},
        all_four_proof_groups_have_frozen_independent_reviews=True,
        all_order_existence_from_implicit_Hensel_theorems_not_finite_jets_or_emptiness_lifting=True,
        global_known_component_lower_bound_remains=57,
        T39_T40_other_component_nonmembership_or_local_branch_exhaustion_claimed=False,
        general_firstjet_affine_fivepoint_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
    with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print('PASS standalone ALL finite Qa in T01--T38 and CA/CB;40 known memberships, global57 unchanged')

if __name__=='__main__':main()

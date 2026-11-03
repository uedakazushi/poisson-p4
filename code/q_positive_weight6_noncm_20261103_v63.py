#!/usr/bin/env python3
"""Full positive-orbit non-CM from exactdimension3, edim8 and multiplicity<=5."""
import argparse
import json
from pathlib import Path
from q_positive_weight6_quadratic_model_20261103_v63 import ROOT,read,sha


def main(output):
    paths=[ROOT/'reports/Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.json',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63/report.json',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_SMALL_MODEL_20261103_v63_REVIEW/report.json',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63/report.json',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63_REVIEW/report.json',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW.json']
    records=[]
    for path in paths:
        d=read(path);assert d['all_passed']
        for filename,digest in d['input_sha256'].items():assert sha(ROOT/filename)==digest,(str(path),filename)
        records.append(d)
    component,quadratic,small,fifth,freview,orbit=records
    dimension=component['actual_full_transverse_local_dimension_exact'];assert dimension==3
    embedding=quadratic['formal_transverse_embedding_dimension'];assert embedding==8
    assert small['exact_implicit95_and_formal8_presentation']
    assert fifth['through5_ideal_dimension']==dimension and fifth['through5_ideal_degree']==5
    assert freview['exact_fifth_cone_dimension3_degree5']
    assert orbit['nonzero_eta_one_complex_projective_orbit_and_closure_Qa8over9']
    cm_lower=embedding-dimension+1;assert cm_lower==6>5
    paths+=[Path(__file__),ROOT/'code/q_positive_weight6_quadratic_model_20261103_v63.py',
            ROOT/'research/notes/Q_POSITIVE_WEIGHT6_NONCM_20261103_v63.md']
    result=dict(id='Q_POSITIVE_WEIGHT6_NONCM_20261103_v63',all_passed=True,
        actual_full_transverse_dimension=3,actual_full_transverse_embedding_dimension=8,
        Hilbert_Samuel_multiplicity_upper_bound=5,Cohen_Macaulay_multiplicity_lower_bound_would_be=6,
        full_actual_Poisson_local_and_completed_rings_non_Cohen_Macaulay_at_positive_point=True,
        all_nonzero_positive_axis_points_full_Poisson_non_Cohen_Macaulay=True,
        projective_local_dimension=25,projective_embedding_dimension=30,
        independent_of_generic_length5_or_generic_quadratic_rank=True,
        reduced_component_nonCM_or_all_other_points_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths})
    with output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print('PASS fullpositiveorbit NON-CM: e<=5 <edim8-dim3+1=6, no genericlength5 used')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)

#!/usr/bin/env python3
"""Pin the independently audited pure commutative-algebra argument."""
import argparse
import json
from pathlib import Path
from s61_nonsquare_structural_reductions_review_20261103_v63 import read,sha

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_QUADRATIC_QUARTIC_QUOTIENT_REFLEXIVE_20261103_v63_REVIEW'


def main(output):
    primary=ROOT/'reports/S61_NONSQUARE_QUADRATIC_QUARTIC_QUOTIENT_REFLEXIVE_20261103_v63.json';author=read(primary)
    for rel,h in author['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    q22=ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63_REVIEW/report.json'
    finite=ROOT/'reports/S61_NONSQUARE_VERTEX_SCROLL_FINITE_20261103_v63_REVIEW/report.json'
    a,b=read(q22),read(finite);assert a['all_passed']and b['all_passed']
    assert a['entire_kernel_module_basis_accepted']and a['product_image_module_locally_free_rank']==146
    assert b['original_light146_weighted152_outside_finite_FULLU_via_degree3_colon']
    assert author['coefficient_regular_ring_dimension']==3 and author['quadratic_product_cokernel_reflexive_rank']==64
    note=ROOT/f'research/notes/{STEM}.md'
    result=dict(id=STEM,independent_mathematical_review_passed=True,
        proof_kind='Pure depth/associated-prime/S2 argument; no new arithmetic',
        exact_projective_resolution_length_at_most1=True,abstract_image_locally_free146_not_ambient_subbundle=True,
        closed_point_depth_at_least2=True,finite_support_torsion_excluded=True,S2_and_torsion_free_imply_reflexive64=True,
        Ext_ge2_zero_and_Ext1_finite_support=True,
        free_flat_ambient146_or_fullquartic_not_implied=True,new_finite_Sigma_rank_proof_not_a_dependency=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),note,primary,q22,finite,ROOT/'research/notes/S61_NONSQUARE_QUADRATIC_QUARTIC_QUOTIENT_REFLEXIVE_20261103_v63.md']})
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in result.items()if k!='input_sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)

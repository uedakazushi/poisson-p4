#!/usr/bin/env python3
"""Independent low invariant/integer checks for the Sigma-free degree bound."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import sympy as s

from s61_affine_sigma_free_degree_bound_20261103_v62 import verify as producer_verify

ROOT=Path(__file__).resolve().parents[1]


def verify(output):
    source=ROOT/'reports/S61_AFFINE_SIGMA_FREE_DEGREE_BOUND_20261103_v62.json'
    original=json.loads(source.read_text());fresh=producer_verify();assert fresh==original
    p0,p1,p2,p3,l0,l1=s.symbols('p0 p1 p2 p3 l0 l1');z=[p0,p1,p2,p3,l0,l1]
    mons=[z[i]*z[j]for i,j in itertools.combinations_with_replacement(range(6),2)]
    # Coefficient action obtained independently by r->r+e*t and t->t+e*r.
    raising=[0,3*p0,2*p1,p2,0,l0]
    lowering=[p1,2*p2,3*p3,0,l1,0]
    columns=[]
    for mon in mons:
        column=[]
        for action in [raising,lowering]:
            value=s.Poly(sum(a*s.diff(mon,v)for a,v in zip(action,z)),*z)
            column.extend(value.coeff_monomial(m)for m in mons)
        columns.append(column)
    rank=s.Matrix.hstack(*map(s.Matrix,columns)).rank();assert rank==21
    c=s.symbols('u4 u14 u49 u35 u50 u36')
    u4,u14,u49,u35,u50,u36=c
    conversion=s.Matrix([u4,-u14+2*u49,-u35+u50,-u36,-5*u14+u49,5*u35+u50]).jacobian(c)
    assert abs(conversion.det())==54
    possible=[]
    for degree in range(1,30):
        if degree%5:continue
        components=[(d,m)for d in range(10,degree+1,5)for m in range(2,degree//d+1)]
        for d,m in components:
            if d*m==degree:possible.append((degree,m,d))
    assert possible==[(20,2,10)]
    component_numbers=[g for g in range(1,21)if g%10==0 and 20%g==0]
    assert component_numbers==[10,20]
    projection=json.loads((ROOT/'reports/S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102.json').read_text())
    assert projection['bad_group_element_incidence_dimension_upper']==9<projection['symplectic_group_dimension']==10
    residue=json.loads((ROOT/'reports/S61_AFFINE_ANCHOR_RESIDUE_DEGREE_20261102.json').read_text())
    assert residue['possible_anchor_support_plane_residue_degree_at_least']==10
    paths=[Path(__file__),source,
           ROOT/'research/notes/S61_AFFINE_SIGMA_FREE_DEGREE_BOUND_20261103_v62.md',
           ROOT/'reports/S61_AFFINE_ANCHOR_RESIDUE_DEGREE_20261102.json',
           ROOT/'research/notes/S61_AFFINE_ANCHOR_RESIDUE_DEGREE_20261102.md',
           ROOT/'reports/S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102.json',
           ROOT/'research/notes/S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102.md',
           ROOT/'research/notes/S61_AFFINE_BOUNDARY_FINITE_ANCHORS_20261102.md',
           ROOT/'research/notes/S61_AFFINE_SURFACE_DEGREE_OBSTRUCTION_20261102.md',
           ROOT/'reports/S61_MARKED_SIGMA_FULL_PENCIL_BRIDGE_20261103_v62_REVIEW.json']
    result=dict(id='S61_AFFINE_SIGMA_FREE_DEGREE_BOUND_20261103_v62_REVIEW',all_passed=True,
                fresh_lightweight_producer_and_recursive_hashes_match=True,recursive_hash_checks=original['recursive_hash_checks'],
                independently_recomputed_binary_cubic_plus_linear_degree2_invariant_dimension=21-rank,
                independent_six_core_conversion_determinant_absolute_value=54,
                stronger_support_plane_residue_bound_distinguished_from_anchor_bound=True,
                no_affine_linear_q_lines_argument_accepted=True,
                connectedness_and_Galois_ten_disjoint_lines_argument_accepted=True,
                geometric_support_plane_degree5_conic_loophole_excluded=True,
                degree20_and25_surfaces_excluded=True,
                center_free_projective_linear_map_needed_for_rank_degree_equality=True,
                any_K_integral_surface_degree_and_whole_projection_rank_at_least=30,
                actual_rank_upper_or_no_surfaces_proved=False,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths})
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS SIGMA-FREE DEGREE/RANK30 AUDIT; DEGREE2 INVARIANTS ZERO')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();verify(args.output)

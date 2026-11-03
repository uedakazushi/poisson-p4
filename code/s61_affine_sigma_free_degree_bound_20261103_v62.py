#!/usr/bin/env python3
"""Check prerequisites for no affine lines and the degree30 surface bound.

Geometric implications are proved in the companion note. In particular the
support-plane residue bound, stronger than the anchor residue bound, is used.
No numerical upper bound on an actual affine projection rank is computed.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def verify():
    dependencies={};visited=set();count=0
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    def audit(path):
        nonlocal count
        path=path.resolve()
        if path in visited:return
        visited.add(path)
        if path.suffix!='.json':return
        data=json.loads(path.read_text())
        if not isinstance(data,dict):return
        for key,base in [('input_sha256',ROOT),('evidence_sha256',path.parent),('output_sha256',path.parent)]:
            for name,expected in data.get(key,{}).items():
                target=base/name;assert sha(target)==expected,(str(path),key,name)
                count+=1
                if key=='input_sha256':audit(target)
    def load(name):
        p=ROOT/'reports'/name;d=json.loads(p.read_text());assert d['all_passed']
        audit(p);dependencies[str(p.relative_to(ROOT))]=sha(p);return d
    whole=load('S61_MARKED_SIGMA_FULL_PENCIL_BRIDGE_20261103_v62.json')
    assert whole['whole_marked_square_stratum_removed_for_general_beta']
    review=load('S61_MARKED_SIGMA_FULL_PENCIL_BRIDGE_20261103_v62_REVIEW.json')
    assert review['whole_marked_square_removed_for_general_beta_accepted']
    anchors=load('S61_AFFINE_BOUNDARY_FINITE_ANCHORS_20261102.json')
    assert anchors['all_generic_true_boundary_points_in_binary_bundle']and anchors['tautological_divisor_ample']
    assert anchors['every_true_boundary_curve_meets_finite_intersection']
    assert anchors['generic_true_boundary_intersection_with_divisor_dimension_upper']==0
    contact=load('S61_AFFINE_SURFACE_DEGREE_OBSTRUCTION_20261102.json')
    assert contact['binary_radical']==['5*u35 + u50','5*u14 - u49']
    assert contact['flag_radical']==['u5','5*u14 - 4*u49']
    assert contact['simple_contact_boundary_curve_only_possible']=='marked square section Sigma'
    assert contact['surface_degree_divisor']==5
    residue=load('S61_AFFINE_ANCHOR_RESIDUE_DEGREE_20261102.json')
    assert residue['possible_anchor_support_plane_residue_degree_at_least']==10
    assert residue['binary_pair_SL2_invariant_polynomials']['2']['kernel_dim']==0
    assert residue['degree_of_divisor_for_core_degree_n']=='5n/2'
    surface=load('S61_AFFINE_SURFACE_BOUND_20261102.json')
    assert surface['generic_affine_dimension_upper']==2 and surface['common_zero_beta_image_dimension_upper']==15
    projection=load('S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102.json')
    assert projection['generic_affine_projection_is_finite']and projection['linear_coordinates']==['u6','u42']
    assert projection['bad_group_element_incidence_dimension_upper']==9
    assert projection['symplectic_group_dimension']==10
    assert projection['common_center_dimension_upper']=={'binary_orbit':3,'flag_orbit':2}
    assert projection['orbit_dimensions']=={'binary_orbit':5,'flag_orbit':4}
    length=load('S61_AFFINE_ANCHOR_LENGTH_BOUND_20261102.json')
    assert length['marked_intersection_degree_divisible_by']==5
    # The geometric proof gives each K-boundary-curve degree >=10 and each
    # hyperplane coefficient >=2. Enumerate only the resulting integer cases.
    below30=[(coefficient*degree,coefficient,degree)
             for degree in range(10,30,5)for coefficient in range(2,30//degree+1)
             if coefficient*degree<30]
    assert below30==[(20,2,10)]
    plane_degrees=[r for r in range(10,11)if 10%r==0];assert plane_degrees==[10]
    geometric_components=[g for g in range(1,21)if g%10==0 and 20%g==0]
    assert geometric_components==[10,20]
    assert [20//g for g in geometric_components]==[2,1]
    c,e,d=20,2,30
    assert d==c+5*e and c%5==0 and c>=5*e and c>=20
    note=ROOT/'research/notes/S61_AFFINE_SIGMA_FREE_DEGREE_BOUND_20261103_v62.md'
    for p in [Path(__file__),note,
              ROOT/'research/notes/S61_AFFINE_ANCHOR_RESIDUE_DEGREE_20261102.md',
              ROOT/'research/notes/S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102.md',
              ROOT/'research/notes/S61_AFFINE_DEGREE_TEN_EXCLUSION_20261102.md']:
        dependencies[str(p.relative_to(ROOT))]=sha(p)
    return dict(id='S61_AFFINE_SIGMA_FREE_DEGREE_BOUND_20261103_v62',all_passed=True,
                base='K=C(W16), with geometric statements after algebraic closure',
                all_recorded_dependency_hashes_current=True,recursive_hash_checks=count,
                general_true_boundary_disjoint_from_marked_square=True,
                no_first_order_dt_direction_at_any_true_boundary_point=True,
                no_affine_linear_q_line_in_general_fiber=True,
                anchor_support_plane_residue_degree_at_least=10,
                each_K_boundary_curve_q_degree_at_least=10,
                degree10_K_boundary_curve='Ten conjugate vertical projective lines over one degree10 support plane',
                degree5_support_plane_conic_case_excluded_by_stronger_support_plane_bound=True,
                sole_numerical_surface_boundary_below30=[{'surface_degree':a,'coefficient':b,'curve_degree':c}for a,b,c in below30],
                degree20_geometric_component_counts=geometric_components,
                degree20_geometric_component_degrees=[20//g for g in geometric_components],
                degree20_excluded_by_absence_of_affine_lines=True,
                affine_surface_degree_over_K_at_least=30,
                exact_fixed_pair_projective_center_avoids_entire_true_boundary=True,
                rank_equals_degree_uses_center_free_linear_projective_map=True,
                two_coordinate_generic_module_rank_if_surface_at_least=30,
                sufficient_rank_upper_bound_for_no_surface=29,
                sufficient_actual_finite_anchor_Artin_length_for_no_surface=19,
                numerically_unexcluded_cycle_data={'marked_intersection_c':20,'degree_over_C_e':2,'surface_degree_d':30},
                no_affine_surfaces_proved=False,general_affine_dimension_upper_remains=2,
                all_true_boundary_curves_avoid_marked_divisor_D_claimed=False,
                actual_rank_or_Artin_length_upper_bound_computed=False,
                whole_affine_fiber_finite_or_five_point_exhaustion_claimed=False,
                input_sha256=dependencies)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify()
    with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print('PASS',report['id'],report['recursive_hash_checks'],'recursive input hashes')

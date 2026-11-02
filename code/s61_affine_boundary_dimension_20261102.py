#!/usr/bin/env python3
"""Audit all reduced infinity types and their contribution to the affine boundary."""
import argparse
import hashlib
import json
import time
import itertools
import sympy as sp
from s61_model import S61, X
from s61_nil22_infinity_20261101 import beta_matrix
from weight1_structure import sp_basis
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_AFFINE_BOUNDARY_DIMENSION_20261102'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_four_orbits():
    x,y,z,w=X[1:];v=sp.Matrix([x,y,z,w]);model=S61()
    out={}
    for name,f,g in [('binary',x*w,0),('flag',0,-x*x),
                     ('canonical',x*z,-x*w),('double',z*z,-2*z*w)]:
        q=sp.zeros(4);q[1,2]=f;q[1,3]=g;q-=q.T
        H=16-len(beta_matrix(model,q).to_DM().rref()[1])
        def coefficients(tensor):
            polynomials=[sp.Poly(tensor[i,j],*v) for i,j in itertools.combinations(range(4),2)]
            return sp.Matrix([poly.coeff_monomial(v[i]*v[j]) for poly in polynomials
                              for i,j in itertools.combinations_with_replacement(range(4),2)])
        columns=[coefficients(q)]
        for A in sp_basis():
            dq=sum(((A*v)[i]*q.diff(v[i]) for i in range(4)),sp.zeros(4))-A*q-q*A.T
            columns.append(coefficients(dq))
        orbit=len(sp.Matrix.hstack(*columns).to_DM().rref()[1])-1
        out[name]={'projective_orbit_dimension':orbit,'beta_linear_kernel_dimension':H,
                   'joint_dimension_upper':H+orbit}
    return out


def certify():
    started=time.monotonic()
    hashes={}
    def report(name):
        path=Path('reports')/(name+'.json')
        record=json.loads((ROOT/path).read_text())
        assert record['all_passed']
        hashes[str(path)]=sha(ROOT/path)
        return record
    four=check_four_orbits()
    assert [v['joint_dimension_upper'] for v in four.values()]==[19,19,17,17],four
    full=report('S61_GENERIC_REDUCED_INFINITY_REVIEW_20261101')
    assert len(full['modular_cases'])==5 and full['corrected_boundary_included']
    mixed=report('S61_R0_MIXED_ORBIT_REDUCTION_20261101')
    orbit_stats={name:d['stats'] for name,d in mixed['orbits'].items()}
    assert set(orbit_stats)=={f'T{i}' for i in range(1,8)}
    for d in orbit_stats.values():assert d['Hdim']+d['orbitdim']==d['joint_rank']<=18
    assert orbit_stats['T3']['beta_rank']==15
    flag=report('S61_FLAG_UNIFORM_OPEN_20261101')
    assert flag['flag']['stats']==dict(Hdim=15,orbitdim=4,beta_rank=16,joint_rank=19)
    flag0=report('S61_FLAG_BINARY_CONNECTION_FORMAL_20261101')
    flag1=report('S61_FLAG_BOUNDARY_RELATIVE_REMOVED_OPEN_20261101')
    assert flag0['nonempty_ambient_open_on_flag_b0_divisor'] and flag1['ambient_remaining_boundary_dimension_upper']==17
    binary=report('S61_R0_INFINITY_20261101')
    assert binary['binary_plane']['G6_identically_zero'] and len(binary['binary_plane']['four_equations'])==4
    rank=report('S61_BINARY_RANK_DROP_FAMILY_20261101')
    assert rank['projective_Q_stabilizer_dimension']==5 and rank['first_jet_image_dimension']==16
    ends=report('S61_BINARY_ENDPOINTS_20261101')
    assert all(d['generic_fixed_beta_stratum_removed'] and d['homogenizer_square_zero'] for d in ends['cases'])
    reduced=report('S61_BINARY_H11_RANK_DROP_RESIDUAL_FINITE_20261102')
    assert reduced['residual_joint_dimension_upper']==16 and reduced['full135_projective_jacobian_rank_upper_in_scope']==34
    local=report('S61_BINARY_LOCAL_DIMENSION_20261101')
    binary65=report('S61_BINARY_FORMAL_RING_20261101')
    assert local['quartic_rank']==210 and local['quartic_minor_determinant_mod_prime']==11420
    assert local['normal_hilbert_function_upper_bound']==[1,7,23,34] and binary65['uniform_normal_length_lower_bound']==65
    vertex=report('NULL_CORRELATION_DUALITY_20261101')['generic_zero_curve']
    assert vertex['smooth'] and vertex['connected'] and vertex['genus']==1 and vertex['degree']==5
    minimal=report('S61_MINIMAL_VERTEX_EXHAUSTION_20261101')
    assert [f['label'] for f in minimal['families']]==['M0','M1']
    chart_counts={f['label']:{k:len(f[k]) for k in ['baseline_free_lines','rational_tails','isolated_axis_charts']} for f in minimal['families']}
    assert sum(v['rational_tails'] for v in chart_counts.values())==18
    assert sum(v['baseline_free_lines'] for v in chart_counts.values())==4
    for name in ['S61_R0_GENERIC_EXHAUSTION_20261101','S61_R0_UNIFORM_VERTEX_20261101',
            'S61_NONBINARY_CANONICAL_EXHAUSTION_20261101','S61_NONBINARY_DOUBLE_ROOT_EXHAUSTION_20261101',
            'S61_R0_MIXED_CURVE_EXHAUSTION_20261101','S61_R0_MIXED_SURFACE_T7_EXHAUSTION_20261101',
            'S61_R0_T5_ZERO_L_EXHAUSTION_20261101','S61_AFFINE_TWENTYFIVE_PLANE_20261102',
            'S61_AFFINE_FINITE_NORMAL_MODULE_20261102']:
        report(name)
    for name in ['S61_GENERIC_REDUCED_INFINITY_EXHAUSTION_20261101','S61_MINIMAL_VERTEX_EXHAUSTION_20261101',
            'S61_R0_MIXED_ORBIT_REDUCTION_20261101','S61_FLAG_UNIFORM_OPEN_20261101',
            'S61_FLAG_BINARY_CONNECTION_FORMAL_20261101','S61_FLAG_BOUNDARY_RELATIVE_REMOVED_OPEN_20261101',
            'S61_BINARY_ENDPOINTS_20261101','S61_BINARY_H11_RANK_DROP_RESIDUAL_FINITE_20261102',STEM]:
        path=Path('research/notes')/(name+'.md')
        hashes[str(path)]=sha(ROOT/path)
    for name in ['s61_model.py','s61_nil22_infinity_20261101.py','s61_infinity_flags.py','weight1_structure.py']:
        path=Path('code')/name;hashes[str(path)]=sha(ROOT/path)
    hashes[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    return dict(id=STEM,all_passed=True,raw_generic_infinity_exact_dimension=3,four_recomputed_orbits=four,
                nonsquare_binary_joint_dimension=14+(10-rank['projective_Q_stabilizer_dimension']),
                flag_joint_dimension=flag['flag']['stats']['joint_rank'],mixed_R0_orbit_data=orbit_stats,
                minimal_chart_counts=chart_counts,minimal_generic_fixed_beta_dimension_upper=1,
                fixed_beta_genuine_affine_boundary_dimension_upper=2,
                binary_rank_at_most_34_residual_joint_dimension_upper=16,flag_residual_joint_dimension_upper=17,
                ordinary_binary_transverse_initial_rank_targets=[5,50,210],
                only_possible_dimension_two_boundary_type='nonsquare binary, full projective Jacobian rank 35',
                generic_affine_fiber_dimension_upper=3,finite_affine_projection_variable_count=3,
                actual_three_linear_forms_exhibited=False,finite_over_beta_proved=False,all_five_exhausted=False,
                input_sha256=hashes,seconds=time.monotonic()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    out=parser.parse_args().output;assert not out.exists()
    result=certify()
    with out.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(f"PASS: raw infinity 3, affine boundary at most 2, affine fiber at most 3; {len(result['input_sha256'])} hashes")

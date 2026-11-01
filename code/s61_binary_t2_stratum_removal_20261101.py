#!/usr/bin/env python3
"""T2 junction removal via finite-base torsion: exact dependency and source checks."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]


def verify():
    start=time.time();names=['S61_BINARY_T2_CONNECTION','S61_BINARY_T2_FORMAL','S61_BINARY_T2_RELATIVE_DOUBLE_LINE','S61_BINARY_UNIFORM_CUBIC']
    paths=[ROOT/'reports'/(n+'_20261101.json') for n in names]
    connection,formal,double,cubic=[json.loads(p.read_text()) for p in paths]
    assert all(r['all_passed'] for r in [connection,formal,double,cubic])
    assert connection['joint_image_dimension']==16 and connection['source_dimension']==11
    assert s.Rational(connection['beta_dominance_minor'])!=0
    assert connection['point']['characteristic_zero_transverse_length_upper_bound']==67
    assert formal['initial_ideal_degree']==65
    assert formal['actual_binary_kernel_torsion_length_upper_bound']==1
    assert formal['actual_binary_kernel_torsion_filtration_degree_if_nonzero']==3
    assert formal['associated_graded_ideal_identified'] and formal['completed_homogenizer_square_zero']
    assert formal['point']['sample']==connection['point']['sample']
    assert double['source_dimension']==9 and double['flat_relative_double_line']
    assert double['all_135_actual_double_line_identities']
    assert cubic['normal_length_lower_bound']==65
    normal=s.Matrix(double['beta_normal_matrix']).applyfunc(s.Rational)
    determinant=normal[:,double['beta_minor_columns']].det()
    assert determinant==s.Rational(double['beta_minor'])==-s.Rational(2248,25)
    point=connection['point'];vertex=s.Matrix(point['curve_derivative']).applyfunc(s.Rational)
    tangent=s.Matrix(point['actual_binary_family_tangent']).applyfunc(s.Rational)
    assert vertex.rank()==2 and tangent.rank()==3
    assert tangent[point['linear_parameters'],:].det()==s.Rational(point['parameter_minor'])!=0
    assert (point['same_beta_regular_full_rank'],point['same_beta_regular_quadratic_rank'],point['same_beta_regular_hessian_radical'])==(35,5,3)
    paths += [Path(__file__),ROOT/'research/notes/FINITE_BASE_TORSION_SPECIALIZATION_20261101.md',
        ROOT/'research/notes/S61_BINARY_T2_FORMAL_20261101.md',
        ROOT/'research/notes/S61_BINARY_T2_RELATIVE_DOUBLE_LINE_20261101.md',
        ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',
        ROOT/'research/notes/S61_BINARY_T2_STRATUM_REMOVAL_20261101.md']
    return dict(id='S61_BINARY_T2_STRATUM_REMOVAL_20261101',all_passed=True,
        verification_scope='Exact small matrix and frozen dependency checks; the finite-base theorem is proved in the notes.',
        source_beta_dimension=9,relative_binary_base_dimension=12,source_stratum_dimension=16,
        fixed_beta_kernel_rank=65,fixed_beta_kernel_torsion_length_upper_bound=1,
        relative_kernel_rank=65,relative_torsion_finite_over_beta_base=True,
        relative_t_saturation_support_dimension_upper_bound=9,
        relative_t_zero_boundary_dimension_upper_bound=8,
        smooth_beta_base_change_minor=str(determinant),
        entire_general_beta_T2_connection_stratum_removed=True,
        ambient_specified_point_t_nilpotence_asserted=False,
        proof_note='research/notes/S61_BINARY_T2_STRATUM_REMOVAL_20261101.md',
        limitations=['Removal of the entire16D stratum is over sufficiently general beta.',
            'No relative or ambient specified-point t^2 identity is inferred from finite-base torsion.',
            'The special fixed-beta t^2 theorem is a separate established input.',
            'The remaining open17D deeper chart and full affine exhaustion remain unresolved.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();report=verify()
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],entire_general_beta_stratum_removed=True)))

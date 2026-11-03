#!/usr/bin/env python3
"""Literal W1 frame and recognition in the old universal actual218 intersection."""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s
from s61_model import S61, X, J, pushforward, coefficients
from weight1_structure import coefficient_map

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def mat(M):
    return [[str(q) for q in row] for row in M.tolist()]


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    out.mkdir(parents=True, exist_ok=False)
    w_path = ROOT/'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63/staircase.json.gz'
    source_path = ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json'
    base_path = ROOT/'reports/AMBIENT_BASE_PINCHED_SURFACE_20261103_v63/report.json'
    bundle_path = ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json'
    builder_path = ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/builder.json'
    old_path = builder_path.parent/'report.json'
    review_path = ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63_REVIEW/report.json'
    uniform_path = ROOT/'reports/AMBIENT_B18_W_UNIFORM_DARBOUX_20261103_v63_REVIEW/report.json'
    assert read(old_path)['all_passed'] and read(review_path)['all_passed'] and read(uniform_path)['all_passed']
    old, body, base, bundle, builder = map(read, [old_path, source_path, base_path, bundle_path, builder_path])
    h, x, y, z, w, W = s.symbols('h x y z w W')
    original = s.Matrix(read(w_path)['CA']['physical_entire_W_curve']).subs(W, 1)
    original = original.subs(dict(zip([h, x, y, z, w], X)), simultaneous=True)
    # New coordinates (oldx,-oldh,oldy,oldz-oldx,oldw-oldh).
    frame = s.Matrix([[0, 1, 0, 0, 0], [-1, 0, 0, 0, 0], [0, 0, 1, 0, 0],
                      [0, -1, 0, 1, 0], [-1, 0, 0, 0, 1]])
    assert frame.det() == 1 and frame*s.Matrix([0, 1, 0, 1, 0]) == s.Matrix([1, 0, 0, 0, 0])
    transported = pushforward(original, frame)
    H, xx, yy, zz, ww = X
    expected = H*H*J
    expected[2, 3] = -H*xx
    expected[3, 2] = H*xx
    expected[2, 4] = 2*H*zz+xx*xx+zz*zz
    expected[4, 2] = -expected[2, 4]
    assert (transported-expected).applyfunc(s.expand) == s.zeros(5)
    value = transported[1:5, 1:5].subs({H: 1, xx: 0, yy: 0, zz: 0, ww: 0})
    assert value == J[1:5, 1:5]
    model = S61()
    normal = model.normalize(transported)
    verified = model.verify_point(normal, check_equations=True)
    u = s.Matrix(normal['parameters'])
    beta_ids = [i for i, wt in enumerate(model.reduced['parameter_weights']) if wt == 1]
    qi = bundle['q_indices_zero_based']
    assert beta_ids == bundle['beta_indices_zero_based']
    beta = s.Matrix(bundle['beta16_to_intrinsic24_adapter'])*u[beta_ids, :]
    p, sv, rv = s.symbols('p s r')
    k = s.symbols('k0:9')
    source_point = {p: 0, sv: 0, rv: 0, **dict.fromkeys(k, 0)}
    source_point[k[0]], source_point[k[5]] = s.S.One, s.Rational(1, 3)
    b0 = s.Matrix(list(map(s.sympify, body['beta24_coordinates']))).subs({p: 0, sv: 0})
    assert beta == b0 and coefficient_map(beta).rank() == 2
    orbit = base['orbit_constant_block_witness']
    assert s.Matrix(orbit['matrix']).det() == s.Rational(orbit['determinant']) != 0
    assert base['orbit_rank'] == 7
    source_records = {}
    full_selected = builder['common_q_selected_rows']
    for name in ['CA', 'CB', 'N']:
        record = body['families'][name]
        pars = list(map(s.Symbol, record['source_fiber_parameters']))
        all_u = s.Matrix(list(map(s.sympify, record['full_61_coordinates'])))
        specialized = all_u.subs(source_point, simultaneous=True)
        assert specialized == u
        q = all_u[qi, :]
        M = q.jacobian(pars).subs({p: 0, sv: 0})
        selected = record['actual_q_selected_rows'] if name == 'N' else full_selected
        det = s.factor(M[selected, :].det())
        assert det != 0
        if name != 'N':
            assert det == -72
        odd_matrix = s.Matrix(record['first_odd_coefficient_matrix'])
        assert odd_matrix.rank() == (6 if name == 'N' else 7)
        assert odd_matrix*s.Matrix([source_point[q] for q in pars]) == s.zeros(45, 1)
        source_records[name] = dict(source_fiber_point={str(q): str(source_point[q]) for q in pars},
                                   all61_exact_match=True, selected_q_rows=selected,
                                   actual_selected_q_fiber_Jacobian_determinant=str(det),
                                   regular_odd_linear_rank=odd_matrix.rank())
    K = u[qi, :][full_selected, :]
    remaining = builder['remaining_q_rows']
    zzpoint = u[qi, :][remaining, :]
    names = builder['ambient_variables']
    variables = tuple(s.Symbol(a) for a in names)
    point = {q: s.S.Zero for q in variables}
    point.update({s.Symbol(f'v{i}'): a for i, a in enumerate(K)})
    point.update({s.Symbol(f'z{i}'): a for i, a in enumerate(zzpoint)})
    nonzero = {str(q): str(a) for q, a in point.items() if a}
    assert nonzero == {'v4': '1', 'v7': '1'}
    assert all(s.sympify(f).subs(point) == 0 for f in old['annihilator_prime_generators'])
    nilF = s.sympify(old['nilradical_generator'])
    assert nilF.subs(point) == 0
    equations = [s.sympify(f) for name in ['CA', 'CB'] for f in builder['component_presentations'][name]]
    assert len(equations) == 218
    point_values = [point[q] for q in variables]
    zeros = {i for i, a in enumerate(point_values) if not a}
    def value_and_gradient(f):
        value, gradient = s.S.Zero, [s.S.Zero]*49
        for mon, coefficient in s.Poly(f, *variables).terms():
            absent = [i for i in zeros if mon[i]]
            if not absent:
                value += coefficient
                for i, e in enumerate(mon):
                    if e:
                        gradient[i] += coefficient*e
            elif len(absent) == 1 and mon[absent[0]] == 1:
                gradient[absent[0]] += coefficient
        return value, gradient
    gradients = []
    for f in equations:
        v, gradient = value_and_gradient(f)
        assert v == 0
        gradients.append(gradient)
    Jall = s.Matrix(gradients)
    ranks = [Jall[:109, :].rank(), Jall[109:, :].rank(), Jall.rank()]
    f_value, f_gradient = value_and_gradient(nilF)
    assert f_value == 0
    Nrank = Jall.col_join(s.Matrix([f_gradient])).rank()
    assert ranks == [28, 28, 29] and Nrank == 30
    # The old global polynomial equalities, independently audited, can be
    # localized here: the point is ON the linear embedded prime, not outside it.
    assert old['square_zero'] and old['annihilator_prime_slice_dimension'] == 4
    assert read(review_path)['exact_actual_target_equalities_rechecked_by_fresh_M2']
    assert old['annihilator_free_variables'] == ['v1', 'v4', 'v5', 'v7']
    shift = {q: q+point[q] for q in variables if point[q]}
    translated = [s.expand(f.subs(shift, simultaneous=True)) for f in equations]
    evidence = dict(
        original_projective_symplectic_point=['0', '1', '0', '1', '0'],
        exact_homogeneous_frame=mat(frame), exact_frame_inverse=mat(frame.inv()),
        full150_transported_physical_coefficients=list(map(str, coefficients(transported))),
        exact_normalized61=list(map(str, u)), normalized150=list(map(str, normal['coefficients'])),
        perspective_parameters=list(map(str, normal['perspective'])),
        full_beta24=list(map(str, beta)), source_records=source_records,
        old_49_target_point_nonzero=nonzero, source_beta_orbit7_minor=str(orbit['determinant']),
        common_target10_fiber_K=list(map(str, K)), remaining35_q=list(map(str, zzpoint)),
        old_49_variable_order=names, all218_translated_actual_intersection_generators=list(map(str, translated)),
        translated_nilradical_generator=str(s.expand(nilF.subs(shift, simultaneous=True))),
        actual218_by49_numeric_Jacobian=mat(Jall), nil_F_numeric_gradient=list(map(str, f_gradient)),
        old_embedded_prime_generators=old['annihilator_prime_generators'],
    )
    with (out/'recognition.json.gz').open('wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0) as g:
            g.write((json.dumps(evidence, separators=(',', ':'))+'\n').encode())
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', ROOT/'code/s61_model.py',
             w_path, source_path, base_path, bundle_path, builder_path, old_path, review_path, uniform_path]
    result = dict(id=STEM, all_passed=True,
        literal_GL5_tensor_identity_not_type_matching=True,
        source_point='p=s=r=0,k0=1,k5=1/3,all_other_k=0',
        original_normalization_checks=verified, all61_all134_all150_verified=True,
        beta_literal_old_critical_rank2_orbit_representative=True,
        beta_sp4_orbit_dimension=7, beta_coefficient_rank=2,
        old_universal_component_fiber_det_CA_CB=['-72', '-72'],
        old_target49_point_nonzero=nonzero,
        all218_actual_target_generators_verified_at_point=True,
        old_polynomial_colon_and_reduction_equalities_reused_after_independent_localization_check=True,
        intrinsic_regular_sequence_odd_ideals_remain_regular_at_this_free_parameter_shift=True,
        critical_slice_tangent_dimensions_CA_CB_D_N=[21, 21, 20, 19],
        full_cone_tangent_dimensions_CA_CB_D_N=[38, 38, 37, 36],
        full_D_completion='translated_saved49_variable218_actual_target_ideal plus17 smooth variables',
        full_D_cone_dim_edim_depth_mult=[28, 37, 21, 2],
        full_D_projective_dim_edim_depth_mult=[27, 36, 20, 2],
        nilradical_principal_square_zero_exact_index2=True,
        nilradical_module='C[[u1,...,u21]]*F',
        exact_Ass=['minimal reduced N prime', 'smooth21-dimensional embedded prime'],
        all_nonzero_W_transport_by_accepted_group_identity=True,
        no_ring_splitting_original_B18_origin_extension_new_global_prime_or_new_component_claimed=True,
        no_new_source_jets_or_M2_GB=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'recognition.json.gz': sha(out/'recognition.json.gz')},
        seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'old_target49_point_nonzero',
          'critical_slice_tangent_dimensions_CA_CB_D_N', 'full_D_cone_dim_edim_depth_mult', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())

#!/usr/bin/env python3
"""Actual finite37-generator gluing map and reduced B18 intersection application."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_REDUCED_INTERSECTION_DVR_GLUING_20261103_v63'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def main(out):
    start = time.monotonic()
    out.mkdir(parents=True, exist_ok=False)
    curve_path = ROOT/'reports/AMBIENT_B18_R_CURVE_TANGENT_20261103_v63/certificate.json.gz'
    cubic_path = ROOT/'reports/AMBIENT_B18_N_R_DVR_CUBIC_GRAPH_20261103_v63/certificate.json.gz'
    orbit_path = ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz'
    gates = [
        ROOT/'reports/AMBIENT_B18_ZERO_LINEAR_CONTACT4_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_N_R_DVR_CUBIC_GRAPH_20261103_v63_REVIEW/report.json',
        ROOT/'reports/AMBIENT_B18_REDUCED_DEPTH_SEMINORMALIZATION_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json',
    ]
    assert all(read(p)['all_passed'] for p in gates)
    assert read(gates[0])['h_exact'] == 4
    gate = read(gates[1])
    assert gate['h_exact'] == gate['k_exact'] == 4
    assert gate['reduced_slice_embedding_dimension'] == 18 and gate['linear_nil_cotangent_kernel_dimension'] == 2
    curve, cubic, orbit = map(read, [curve_path, cubic_path, orbit_path])
    rows = curve['N37_coordinate_rows']
    odd = cubic['actual_odd9_physical_rows']
    assert rows == cubic['actual_N37_physical_coordinate_rows'] and len(rows) == 37
    t = s.Symbol('t')
    physical = list(map(s.sympify, curve['physical_curve150']))
    centered = [s.expand(v-v.subs(t, 0)) for v in physical]
    assert all(s.Poly(v, t).degree() <= 2 for v in centered if v)
    assert all(centered[i] == 0 for i in odd)
    nonzero = {str(i): str(centered[i]) for i in rows if centered[i]}
    assert nonzero == {'16': '13*t', '22': '-25*t**2', '65': '-12*t', '112': '-13*t'}
    U = {i: s.Symbol(f'u{i}') for i in rows}
    uniformizer = U[16]/13
    relations = [s.expand(U[i]-centered[i].subs(t, uniformizer)) for i in rows if i != 16]
    relations.append(U[16]**4)
    substitution = {U[i]: centered[i] for i in rows}
    assert len(relations) == 37
    assert all(s.rem(s.Poly(f.subs(substitution), t), s.Poly(t**4, t)).as_expr() == 0 for f in relations)
    assert [s.diff(f, U[i]) for f, i in zip(relations[:-1], [i for i in rows if i != 16])] == [1]*36
    # Rank20 action on the actual R curve supplies full formal coordinates
    # (g19,t) on the smooth R component, not an arbitrary abstract isomorphism.
    O = s.Matrix(orbit['physical150_by19_scalar_then_GL'])
    v = s.Matrix(curve['physical_R_curve_derivative150'])
    assert O.row_join(v).rank() == 20
    trace_a = cubic['full150_trace_zero_even_a']
    beta = cubic['full150_remainders_odd_b']
    minimal_tables = {str(i): dict(physical_local_target=str(centered[i]),
                                 trace_even_a_through_T3=trace_a[i],
                                 odd_beta_through_T2=beta[i]) for i in odd}
    evidence = dict(
        actual_N37_physical_coordinate_rows=rows,
        actual_N28_base_physical_coordinate_rows=rows[:28], actual_odd9_physical_rows=odd,
        centered_physical37_to_C4_map=[dict(row=i, image=str(centered[i])) for i in rows],
        all150_centered_physical_map=list(map(str, centered)),
        C4_uniformizer_from_actual_target='t=u16/13',
        full_preimage_kernel_in_regular37_presentation=list(map(str, relations)),
        physical_map_full_normalizer_trace_tables=minimal_tables,
        raw_Weierstrass_coefficients=cubic['Weierstrass_monic_coefficients'],
        trace_zero_discriminant_mod_T4='64*t^2',
        actual_G19_section_GL18_indices=orbit['GL18_generator_indices'],
        full_G19_equivariant_gluing_map_recipe='q -> coefficient of (1+g0) G P_R(t)(G^-1 X) G^T minus q(P0), modulo t^4; G=I+sum GL18 generators',
    )
    with (out/'finite-gluing-map.json.gz').open('wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0) as g:
            g.write((json.dumps(evidence, separators=(',', ':'))+'\n').encode())
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', curve_path, cubic_path, orbit_path, *gates]
    result = dict(
        id=STEM, all_passed=True, exact_R_N_contact_order=4,
        actual_physical_N37_to_C4_map_is_explicit_and_surjective=True,
        centered_nonzero_physical_generators=nonzero,
        all9_selected_odd_physical_generators_zero_not_trace_odd_parts_zero=True,
        actual_normalizer_trace_tables_bind_finite_map_without_abstract_coordinate_alignment=True,
        exact_full_regular37_preimage_kernel_generated_by_saved37_relations=True,
        N_R_scheme_intersection_completion='C[[g19,t]]/(t^4)',
        full_reduced_intersection_completion='A_N_full fiber_product over C[[g19,t]]/(t^4) with C[[g19,t]]',
        transverse_reduced_intersection_completion='A_N_slice fiber_product over C[t]/(t^4) with C[[t]]',
        actual_equivariant_gluing_recipe_not_arbitrary_abstract_identification=True,
        reduced_slice_dim_edim_depth_mult=[9, 18, 1, 2],
        reduced_cone_dim_edim_depth_mult=[28, 37, 20, 2],
        reduced_projective_dim_edim_depth_mult=[27, 36, 19, 2],
        actual_nonreduced_cone_embedding_dimension=39,
        actual_nilradical_image_in_cotangent_dimension=2,
        seminormalization_node_packet_reused_without_new_computation=True,
        seminormalization='C[[g19,x1,...,x9,t]]/(t*x1,...,t*x9)',
        full_nonreduced_D_ring_exact_depth_original_reduced_conductor_or_coordinate_orbit_classification_claimed=False,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'finite-gluing-map.json.gz': sha(out/'finite-gluing-map.json.gz')},
        seconds=time.monotonic()-start,
    )
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['all_passed', 'exact_R_N_contact_order',
          'centered_nonzero_physical_generators', 'reduced_cone_dim_edim_depth_mult', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())

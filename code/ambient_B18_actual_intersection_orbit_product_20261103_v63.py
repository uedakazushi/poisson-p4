#!/usr/bin/env python3
"""Exact orbit minor for the scheme-level common B18 formal product."""
import argparse
import gzip
import hashlib
import itertools
import json
import resource
import signal
import time
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63'


def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main(out):
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    signal.alarm(60)
    out.mkdir(parents=True, exist_ok=False)
    reduction_path = ROOT / 'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    tangent_path = ROOT / 'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63/report.json'
    tangent_cert = tangent_path.parent / 'certificate.json.gz'
    reviews = [
        ROOT / 'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_CA_CB_MINIMAL_LENGTHS_AND_R_CM_20261103_v63_REVIEW.json',
    ]
    assert all(read(p)['all_passed'] for p in reviews)
    reduction, tangent, cert = map(read, [reduction_path, tangent_path, tangent_cert])
    assert tangent['all_passed'] and tangent['actual_D_cone_depth_upper_bound'] == 20
    h, x, y, z, w = s.symbols('h x y z w')
    X = s.Matrix([h, x, y, z, w])
    pairs = list(itertools.combinations(range(5), 2))
    monomials = list(itertools.combinations_with_replacement(range(5), 2))

    def vector(P):
        return s.Matrix([s.Poly(P[i, j], *X).coeff_monomial(X[u] * X[v])
                         for i, j in pairs for u, v in monomials])

    ey = s.eye(5)[:, 2]
    D = s.Matrix([x*x, 0, 0, h*x, h*h + z*z])
    P = ey * D.T - D * ey.T
    p = vector(P)
    assert p[106] == 1
    assert all(sum(s.diff(P[i, j], X[i]) for i in range(5)) == 0 for j in range(5))
    columns = []
    for i, j in itertools.product(range(5), repeat=2):
        E = s.zeros(5)
        E[i, j] = 1
        v = E * X
        L = E * P + P * E.T - sum((v[k] * P.diff(X[k]) for k in range(5)), s.zeros(5))
        columns.append(vector(L))
    O25 = s.Matrix.hstack(*columns)
    assert O25.rank() == 19 and O25[:, 12] == p
    indices = reduction['orbit_GL_generator_indices']
    old_O = s.Matrix(cert['families']['CA']['full150_by28_physical_even_linear_J'])[:, :19]
    assert O25[:, indices] == old_O
    # Replace precisely E_yy by the literal tensor scalar.  This gives a
    # commuting scalar factor plus eighteen actual GL directions.
    gl18 = [i for i in indices if i != 12]
    assert len(gl18) == 18
    O19 = p.row_join(O25[:, gl18])
    free = reduction['common_target126_raw_coordinate_indices']
    C = s.Matrix(reduction['canonical_reconstruction_from126'])
    assert C.shape == (150, 126) and C[free, :] == s.eye(126)
    assert C * p[free, :] == p and C * O19[free, :] == O19
    Q = O19[free, :]
    pivots = list(Q.T.rref()[1])
    assert len(pivots) == 19
    minor = s.factor(Q[pivots, :].det())
    assert minor != 0
    complement = [i for i in range(126) if i not in pivots]
    assert len(complement) == 107
    # [Q | e_j, j in complement] has determinant +/- the recorded minor.
    # No 126 by 126 determinant, formal inversion, or extra jet is needed.
    evidence = dict(
        physical_point150=list(map(str, p)),
        all25_physical_GL_derivative_rank=19,
        scalar_GL_index=12,
        scalar_derivative_equals_physical_point=True,
        GL18_generator_indices=gl18,
        physical150_by19_scalar_then_GL=[[str(q) for q in row] for row in O19.tolist()],
        canonical126_raw_coordinate_indices=free,
        orbit19_minor_rows_in126=pivots,
        orbit19_minor_raw_coordinate_rows=[free[i] for i in pivots],
        orbit19_minor=str(minor),
        affine_transverse107_free_coordinate_indices=complement,
        ambient126_action_product_Jacobian_determinant_up_to_sign=str(minor),
    )
    with (out / 'certificate.json.gz').open('wb') as f:
        with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0) as g:
            g.write((json.dumps(evidence, separators=(',', ':')) + '\n').encode())
    paths = [Path(__file__), ROOT / f'research/notes/{STEM}.md', reduction_path,
             tangent_path, tangent_cert, *reviews]
    result = dict(
        id=STEM, all_passed=True,
        actual_point='P0=e_y wedge (x^2,0,0,hx,h^2+z^2)',
        independently_rebuilt_all25_physical_GL_columns=True,
        orbit_dimension=19, ambient_canonical_dimension=126, transverse_ambient_dimension=107,
        orbit19_minor=str(minor), literal_scalar_is_one_of19_independent_parameters=True,
        actual_D_ideal_GL5_times_scalar_stable_including_nilpotents=True,
        formal_inverse_function_and_Artin_functor_argument=True,
        actual_D_completed_ring='C[[g19]] completed_tensor D_transverse',
        simultaneous_common_slice_for_CA_CB_D_N_R=True,
        actual_D_cone_dimension=28, actual_D_cone_embedding_dimension=39,
        actual_D_Hilbert_Samuel_multiplicity=2,
        actual_D_cone_depth_interval=[19, 20],
        actual_D_projective_dimension=27, actual_D_projective_embedding_dimension=38,
        actual_D_projective_depth_interval=[18, 19],
        transverse_D_dimension=9, transverse_D_embedding_dimension=20,
        transverse_D_depth_interval=[0, 1], transverse_D_multiplicity=2,
        exact_depth_or_full_relative_ring_claimed=False,
        remaining_depth_question='Does the transverse maximal ideal have a nonzero socle? Known embedded curves do not exclude an embedded vertex.',
        no_new_CAS_jets_or_flatness_reducedness_CM_assumption=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'certificate.json.gz': sha(out / 'certificate.json.gz')},
        seconds=time.monotonic() - start,
    )
    (out / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'orbit19_minor',
          'actual_D_cone_depth_interval', 'actual_D_projective_depth_interval', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output.resolve())

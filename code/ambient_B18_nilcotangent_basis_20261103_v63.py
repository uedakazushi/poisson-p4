#!/usr/bin/env python3
"""Exact physical cotangent basis of the nilradical at generic B18."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_NILCOTANGENT_BASIS_20261103_v63'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def main(out):
    start = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    out.mkdir(parents=True, exist_ok=False)
    np = ROOT/'reports/AMBIENT_B18_R_CURVE_TANGENT_20261103_v63/certificate.json.gz'
    dp = ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63/certificate.json.gz'
    kp = ROOT/'reports/AMBIENT_B18_N_R_DVR_CUBIC_GRAPH_20261103_v63_REVIEW/report.json'
    n, d, k = map(read, [np, dp, kp])
    assert k['all_passed'] and k['k_exact'] == 4
    assert k['reduced_cone_embedding_dimension'] == 37
    assert k['linear_nil_cotangent_kernel_dimension'] == 2
    N = s.Matrix(n['full150_by37_actual_N_normalized_tangent'])
    D = s.Matrix(d['common150_by39_normalized_tangent'])
    rows = d['common39_physical_coordinate_rows']
    assert D[rows, :] == s.eye(39)
    assert N[n['N37_coordinate_rows'], :] == s.eye(37)
    coordinates = N[rows, :]
    assert D*coordinates == N
    vectors = coordinates.T.nullspace()
    assert len(vectors) == 2
    physical = [{str(rows[i]): str(a) for i, a in enumerate(v) if a}
                for v in vectors]
    assert physical == [{'46': '1', '95': '1'}, {'63': '1'}]
    assert s.Matrix.hstack(*vectors).rank() == 2
    # The radial direction has zero value under both linear forms.
    rp = ROOT/'reports/AMBIENT_B18_R_CURVE_TANGENT_20261103_v63/report.json'
    assert read(rp)['all_passed']
    q = s.symbols('t')
    P0 = s.Matrix([s.sympify(a).subs(q, 0) for a in n['physical_curve150']])
    assert P0[106] == 1 and all((v.T*P0[rows, :])[0] == 0 for v in vectors)
    evidence = dict(common39_rows=rows, N_coordinates_in_D39=coordinates.tolist(),
                    physical_linear_forms=physical)
    (out/'certificate.json').write_text(json.dumps(evidence, default=str, indent=2)+'\n')
    inputs = [Path(__file__), ROOT/f'research/notes/{STEM}.md', np, dp, kp, rp]
    report = dict(id=STEM, all_passed=True, actual_contact_k=4,
                  all150_N_tangent_columns_in_actual_D_tangent_verified=True,
                  actual_reduced_tangent_equals_T_N37_by_exact_edim=True,
                  nil_cotangent_kernel_dimension=2, physical_linear_forms=physical,
                  basis_means_classes_with_formal_nilpotent_lifts=True,
                  raw_linear_functions_themselves_claimed_nilpotent=False,
                  same_basis_after_projective_chart_q106_equals1=True,
                  full_nilradical_generators_nil_index_or_depth_claimed=False,
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs},
                  evidence_sha256={'certificate.json': sha(out/'certificate.json')},
                  seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({x: report[x] for x in ['all_passed', 'nil_cotangent_kernel_dimension',
                                          'physical_linear_forms', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args().output.resolve())

#!/usr/bin/env python3
"""Strict stabilizer matrix review; the smooth-curve finiteness proof is geometric."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from s61_model import S61, X, quotient_coefficients

ROOT = Path(__file__).resolve().parents[1]


def verify():
    started = time.monotonic(); x0, x, y, z, w = X
    P = s.zeros(5)
    P[1, 2] = P[3, 4] = x0*x0; P[2, 3] = -x*z; P[2, 4] = x*w; P -= P.T
    model = S61(); normalized = model.normalize(P); model.verify_point(normalized)
    assert all(normalized['parameters'][i] == 0 for i, degree in enumerate(model.reduced['parameter_weights']) if degree == 1)
    columns = []
    for i in range(5):
        for j in range(5):
            A = s.zeros(5); A[i, j] = 1
            columns.append(s.Matrix(quotient_coefficients(A*P+P*A.T-X[j]*P.diff(X[i]))))
    strict = s.Matrix.hstack(*columns)
    assert strict.rank() == 21
    assert strict.row_join(s.Matrix(quotient_coefficients(P))).rank() == 21
    # Explicit kernel: identity, d_y, x*d_y, z*d_z-w*d_w.
    e0 = s.eye(5); e1 = s.zeros(5); e1[2, 0] = 1
    e2 = s.zeros(5); e2[2, 1] = 1; e3 = s.diag(0, 0, 0, 1, -1)
    kernel = s.Matrix.hstack(*[M.reshape(25, 1) for M in [e0, e1, e2, e3]])
    assert kernel.rank() == 4 and strict*kernel == s.zeros(150, 4)
    assert all(A*B == B*A for A in [e1, e2, e3] for B in [e1, e2, e3])
    paths = [Path(__file__), ROOT/'code/s61_model.py',
        ROOT/'reports/TORIC_STRICT_STABILIZER_20261101.json',
        ROOT/'research/notes/TORIC_STRICT_STABILIZER_20261101.md',
        ROOT/'research/notes/SMOOTH_CURVE_FINITE_STABILIZER_20261101.md']
    return {'id': 'STABILIZER_OBSTRUCTIONS_REVIEW_20261101', 'all_passed': True,
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'strict_action_rank': 21, 'projectivized_action_with_scalar_rank': 21,
        'strict_stabilizer_dimension': 3, 'projective_tensor_stabilizer_dimension': 4,
        'explicit_strict_kernel_mod_identity': ['d_y', 'x*d_y', 'z*d_z-w*d_w'],
        'accepted_geometric_claims': ['Every T01 point has a dimension4 abelian Cartan-limit subalgebra in its strict stabilizer.',
            'On Delta!=0 both vector and projective Sp4 stabilizers of beta are finite.'],
        'scope': 'No limiting wedge-image membership assertion; no assertion that an arbitrary modular matrix stabilizes beta.',
        'proof_note': 'research/notes/STABILIZER_OBSTRUCTIONS_REVIEW_20261101.md',
        'elapsed_seconds': round(time.monotonic()-started, 3),
        'input_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); report = verify(); args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out: json.dump(report, out, indent=2); out.write('\n')
    print(json.dumps({'all_passed': True, 'seconds': report['elapsed_seconds']}))


if __name__ == '__main__': main()

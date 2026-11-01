#!/usr/bin/env python3
"""A strict-stabilizer obstruction beyond the projective dimension condition."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s
from s61_model import S61, X, quotient_coefficients
from toric_zero_jet_dimension_20261101 import witness, encode

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    x0, x, y, z, w = X
    P = s.zeros(5)
    P[1, 2] = P[3, 4] = x0*x0
    P[2, 3], P[2, 4] = -x*z, x*w
    P -= P.T
    model = S61()
    normalized = model.normalize(P)
    model.verify_point(normalized)
    assert all(normalized['parameters'][i] == 0 for i, degree in enumerate(model.reduced['parameter_weights']) if degree == 1)
    columns = []
    for i in range(5):
        for j in range(5):
            a = s.zeros(5)
            a[i, j] = 1
            columns.append(s.Matrix(quotient_coefficients(a*P+P*a.T-X[j]*P.diff(X[i]))))
    strict = s.Matrix.hstack(*columns)
    projective = strict.row_join(s.Matrix(quotient_coefficients(P)))
    strict_witness, projective_witness = witness(strict), witness(projective)
    assert strict_witness['rank'] == projective_witness['rank'] == 21
    kernel = strict.nullspace()
    assert len(kernel) == 4  # The scalar identity is the zero projective vector field.
    report = dict(id='TORIC_STRICT_STABILIZER_20261101', all_passed=True, field='QQ',
                  potential='G=x*z*w', normalized_S61_coordinates=list(map(str, normalized['parameters'])),
                  strict_action_rank_witness=strict_witness,
                  projective_action_rank_witness=projective_witness,
                  strict_kernel_matrices=[encode(s.Matrix(5, 5, list(v))) for v in kernel],
                  strict_projective_vector_field_stabilizer_dimension=3,
                  projectivized_tensor_stabilizer_dimension=4,
                  consequence='This genuine zero-first-jet Poisson point is outside T01, although it satisfies the weaker projective-stabilizer-dimension>=4 condition.',
                  universal_necessary_condition='Every T01 point has a four-dimensional abelian Cartan-limit subalgebra in its strict Poisson stabilizer.',
                  limitations=['No wedge-square image-membership assertion at the limiting Lie subalgebra: evaluation can lose rank.',
                               'No classification of the remaining lower-dimensional zero-first-jet strata.'],
                  proof_note='research/notes/TORIC_STRICT_STABILIZER_20261101.md',
                  source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in (Path(__file__), ROOT/'code/s61_model.py')})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report, out, indent=2)
        out.write('\n')
    print(json.dumps({'all_passed': True, 'strict_stabilizer_dimension': 3, 'outside_T01': True}))


if __name__ == '__main__':
    main()

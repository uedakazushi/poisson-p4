#!/usr/bin/env python3
"""Independent saved polynomial-basis and active-algebra certificate checks."""
import argparse
import hashlib
import json
from itertools import combinations_with_replacement as cwr
from pathlib import Path
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def verify():
    start = time.time()
    primary = 'reports/S61_FLAG_TRUE_CORE_QUADRATIC_20261101.json'
    data = json.loads((ROOT/primary).read_text())
    assert data['all_passed']
    assert all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
               for path, digest in data['input_sha256'].items())
    Q = s.Matrix(data['polynomial_quadratic_basis']).applyfunc(s.sympify)
    columns = data['constant_minor']['columns']
    assert Q.shape == (15, 55) and Q[:, columns] == s.eye(15)
    assert s.Rational(data['constant_minor']['value']) != 0
    assert all(s.denom(value).is_Number for value in Q)
    variables = list(map(s.Symbol, data['variables']))
    assert max(s.Poly(value, *variables).total_degree() for value in Q) == 2
    C = s.Matrix(data['constant_active_quadrics']).applyfunc(s.Rational)
    active = data['active_variables']
    pairs = list(cwr(active, 2))
    all_pairs = list(cwr(range(10), 2))
    embedded = s.zeros(6, 55)
    for j, mon in enumerate(pairs):
        embedded[:, all_pairs.index(mon)] = C[:, j]
    assert (embedded[:, columns]*Q-embedded).applyfunc(s.expand) == s.zeros(6, 55)
    assert C.rank() == 6
    triples = list(cwr(active, 3))
    products = s.zeros(24, 20)
    for i in range(6):
        for j, variable in enumerate(active):
            for k, mon in enumerate(pairs):
                products[4*i+j, triples.index(tuple(sorted((*mon, variable))))] += C[i, k]
    selected = data['active_cubic_product_selected_rows']
    inverse = s.Matrix(data['active_cubic_product_inverse']).applyfunc(s.Rational)
    assert inverse*products[selected, :] == s.eye(20)
    assert data['raw_t2_nonzero'] and data['t2_in_quadratic_ideal']
    paths = list(data['input_sha256'])+[primary,
        'research/notes/S61_FLAG_TRUE_CORE_QUADRATIC_20261101.md',
        str(Path(__file__).relative_to(ROOT)),
        'research/notes/S61_FLAG_TRUE_CORE_QUADRATIC_REVIEW_20261101.md']
    return dict(id='S61_FLAG_TRUE_CORE_QUADRATIC_REVIEW_20261101', all_passed=True,
                independent_polynomial_basis_pivot_check=True,
                independent_active_quadrics_membership=True,
                independent_active_cubic_inverse=True,
                proof_and_source_review_accepted=True,
                accepted_hilbert_function_through2=[1, 10, 40], t_order_at_least=3,
                active_algebra_hilbert_function=[1, 4, 4], active_algebra_length=9,
                independent_raw135_quadratic_or_t2_reimplementation=False,
                limitations=['Raw A2 reconstruction, the constant minor and t2 membership use the source-audited primary producer.',
                             'The active algebra length9 is not an ambient transverse length.',
                             'No universal higher initial ideal, ambient nilpotence or saturation removal.'],
                input_sha256={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
                seconds=time.time()-start)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print('TRUE FLAG CORE QUADRATIC REVIEW PASS', result['seconds'])

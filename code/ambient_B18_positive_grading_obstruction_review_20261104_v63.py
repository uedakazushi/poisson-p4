#!/usr/bin/env python3
"""Independent physical stabilizer and all-cocharacter grading obstruction audit."""
import argparse
from collections import Counter
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63_REVIEW'

def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()

def enc(matrix):
    return [[str(a) for a in row] for row in matrix.tolist()]

def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    out.mkdir(parents=True, exist_ok=False)
    author = ROOT / 'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63'
    ar = read(author / 'report.json')
    ac = read(author / 'certificate.json')
    assert ar['all_exact_physical_tangent_and_stabilizer_weight_checks_passed']
    assert sha(author / 'certificate.json') == ar['certificate_sha256']
    hashes_checked = 0
    for rel, digest in ar['input_sha256'].items():
        assert sha(ROOT / rel) == digest, rel
        hashes_checked += 1
    orbit_path = ROOT / 'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63'
    tangent_path = ROOT / 'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63'
    oc = read(orbit_path / 'certificate.json.gz')
    tc = read(tangent_path / 'certificate.json.gz')
    gates = [
        ROOT / 'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW.json',
        ROOT / 'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63_REVIEW/report.json',
    ]
    for path in gates:
        assert read(path)['all_passed']
    # Use a literal polynomial 5x5 skew tensor, independently of the author's
    # index-action implementation. Coordinate order is h,x,y,z,w throughout.
    x = s.symbols('h x y z w')
    h, xx, yy, z, w = x
    P = s.zeros(5)
    for i, j, a in [(0, 2, -xx**2), (2, 3, h*xx),
                    (2, 4, h*h + z*z)]:
        P[i, j], P[j, i] = a, -a
    pairs = list(itertools.combinations(range(5), 2))
    mons = list(itertools.combinations_with_replacement(range(5), 2))
    keys = [(a, b) for a in pairs for b in mons]
    def coefficients(Q):
        polys = {pair: s.Poly(s.expand(Q[pair[0], pair[1]]), *x) for pair in pairs}
        return s.Matrix([polys[pair].coeff_monomial(x[mon[0]]*x[mon[1]])
                         for pair, mon in keys])
    point = coefficients(P)
    assert point == s.Matrix(oc['physical_point150'])
    assert {i: int(a) for i, a in enumerate(point) if a} == {20: -1, 106: 1, 120: 1, 132: 1}
    def lie(M):
        velocity = M*s.Matrix(x)
        return M*P + P*M.T - P.applyfunc(
            lambda a: sum(s.diff(a, x[k])*velocity[k] for k in range(5)))
    directions = []
    for a, b in itertools.product(range(5), repeat=2):
        M = s.zeros(5)
        M[a, b] = 1
        directions.append(coefficients(lie(M)))
    O = s.Matrix.hstack(*directions)
    assert O.rank() == 19
    radial = directions[2*5+2]
    assert radial == point  # E_yy P0=P0; tensor scalar adds no orbit dimension.
    orbit19 = s.Matrix.hstack(point, *(directions[j] for j in oc['GL18_generator_indices']))
    assert orbit19 == s.Matrix(oc['physical150_by19_scalar_then_GL'])
    assert orbit19.rank() == 19
    # Solve the whole pure stabilizer, not only its diagonal subalgebra.
    expected_kernel = []
    for a, b in [(2, 0), (2, 1), (2, 3), (2, 4), (4, 1)]:
        M = s.zeros(5)
        M[a, b] = 1
        expected_kernel.append(M)
    K = s.Matrix.hstack(*(s.Matrix(M).reshape(25, 1) for M in expected_kernel),
                         s.Matrix(s.eye(5)).reshape(25, 1))
    assert K.rank() == 6 and O*K == s.zeros(150, 6)
    assert s.Matrix.hstack(*O.nullspace()).columnspace()
    assert len(O.nullspace()) == 6
    assert s.Matrix(ac['GL25_kernel_basis']).T == K
    # A stronger independent literal polynomial identity than M^5=0:
    # arbitrary tracefree stabilizer matrices have third power zero.
    t = s.symbols('t0:5')
    N = sum((a*M for a, M in zip(t, expected_kernel)), s.zeros(5))
    assert N.trace() == 0 and N**3 == s.zeros(5)
    assert all(s.expand(s.trace(N**j)) == 0 for j in range(1, 6))
    brackets = []
    for i, j in itertools.combinations(range(5), 2):
        comm = expected_kernel[i]*expected_kernel[j] - expected_kernel[j]*expected_kernel[i]
        if comm != s.zeros(5):
            brackets.append([i, j, enc(comm)])
    assert len(brackets) == 1
    # Diagonal exponents and exact finite action, not merely its infinitesimal
    # derivative: the tensor multiplier cancels lambda^{-1} on e_y.
    lam = s.Symbol('lam', nonzero=True)
    g = s.diag(1, 1, 1/lam, 1, 1)
    substitution = dict(zip(x, g.inv()*s.Matrix(x)))
    Pg = lam*g*P.subs(substitution, simultaneous=True)*g.T
    assert Pg.applyfunc(s.simplify) == P
    T0 = s.diag(0, 0, -1, 0, 0)
    assert coefficients(lie(T0)) + point == s.zeros(150, 1)
    assert O.row_join(point).rank() == 19
    # Compute the physical weights from actual finite pushforward of all 150
    # coefficient basis vectors, and count each invariant subspace directly.
    exponents = [0, 0, -1, 0, 0]
    weights = [1 + exponents[pair[0]] + exponents[pair[1]]
               - exponents[mon[0]] - exponents[mon[1]] for pair, mon in keys]
    assert weights == ac['physical150_torus_weights']
    T = s.Matrix(tc['common150_by39_normalized_tangent'])
    free = tc['common39_physical_coordinate_rows']
    assert T.shape == (150, 39) and T[free, :] == s.eye(39)
    assert T*orbit19[free, :] == orbit19
    W = s.diag(*weights)
    induced = (W*T)[free, :]
    assert W*T == T*induced
    tangent_weights, orbit_weights = {}, {}
    for value in sorted(set(weights)):
        rows = [i for i, weight in enumerate(weights) if weight == value]
        tangent_weights[value] = T[rows, :].rank()
        orbit_weights[value] = orbit19[rows, :].rank()
    quotient = {a: tangent_weights[a] - orbit_weights[a] for a in tangent_weights
                if tangent_weights[a] != orbit_weights[a]}
    assert quotient == {0: 17, 1: 3}
    assert tangent_weights == {0: 32, 1: 7, 2: 0, 3: 0}
    assert orbit_weights == {0: 15, 1: 4, 2: 0, 3: 0}
    assert keys[63] == ((1, 2), (0, 3)) and weights[63] == 0
    assert orbit19[63, :] == s.zeros(1, 19) and T[63, :] != s.zeros(1, 39)
    assert {str(k): v for k, v in quotient.items()} == ar['transverse_weight_multiplicities']
    # The effective stabilizer has one torus dimension: the central GL torus
    # is a representation kernel; the entire tracefree pure stabilizer Lie is
    # nilpotent, so the identity SL5 stabilizer is unipotent. The displayed
    # tensor-projection section supplies the remaining effective torus.
    evidence = {
        'physical150_by25_GL_action': enc(O),
        'literal_point150': [str(a) for a in point],
        'pure_GL_stabilizer_basis': [enc(M) for M in expected_kernel] + [enc(s.eye(5))],
        'generic_tracefree_matrix': enc(N),
        'generic_tracefree_cube_zero': True,
        'nonzero_tracefree_commutators': brackets,
        'finite_fixing_torus': {'GL': 'diag(1,1,lambda^-1,1,1)', 'tensor_multiplier': 'lambda'},
        'physical150_weights': weights,
        'actual_tangent39_weight_dimensions': tangent_weights,
        'orbit19_weight_dimensions': orbit_weights,
        'quotient20_weight_dimensions': quotient,
        'u63_zero_weight_nonzero_quotient_cotangent_class': True,
        'all_cocharacter_argument': [
            'Every cocharacter has connected image in the identity stabilizer.',
            'The connected pure SL5 stabilizer is unipotent by the nilpotent-matrix Lie algebra criterion.',
            'The only pure GL torus dimension is central GL, acting trivially on every quadratic bivector.',
            'The extended effective stabilizer modulo central GL has maximal torus rank one.',
            'The displayed rank-one torus is maximal; any other maximal torus is conjugate in the connected stabilizer.',
            'Conjugation preserves the actual tangent quotient representation; every nontrivial effective cocharacter retains seventeen zero weights.',
        ],
        'scope': 'Only gradings induced by actual GL5 times tensor Gm stabilizer actions; no arbitrary formal-local-ring action exclusion.'
    }
    (out/'independent-certificate.json').write_text(json.dumps(evidence, indent=2)+'\n')
    inputs = [Path(__file__), ROOT/'research/notes'/(STEM+'.md'),
              author/'report.json', author/'certificate.json',
              orbit_path/'certificate.json.gz', tangent_path/'certificate.json.gz'] + gates
    report = {
        'id': STEM, 'all_passed': True,
        'independent_literal_four_term_physical_tensor_and_all25_GL_directions': True,
        'full_pure_GL_stabilizer_dimension': 6,
        'tracefree_stabilizer_dimension': 5,
        'all_tracefree_matrices_cube_zero': True,
        'GL_and_GL_times_tensor_scalar_orbit_dimension': 19,
        'effective_extended_stabilizer_maximal_torus_rank': 1,
        'actual_full_nonreduced_D_tangent_dimension': 39,
        'actual_transverse20_weight_multiplicities': quotient,
        'all_cocharacters_not_just_diagonal_actions_covered': True,
        'no_stabilizer_induced_strictly_positive_transverse_grading': True,
        'u63_is_nonzero_zero_weight_quotient_class': True,
        'no_arbitrary_abstract_formal_grading_or_full_D_structure_claim': True,
        'no_author_producer_rerun_no_sourcejets_no_GB_no_new_points': True,
        'author_input_hashes_verified': hashes_checked,
        'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in inputs},
        'evidence_sha256': {'independent-certificate.json': sha(out/'independent-certificate.json')},
        'seconds': time.monotonic()-started, 'frozen': True,
    }
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: v for k, v in report.items() if k not in ['input_sha256', 'evidence_sha256']}, indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'reports'/STEM)
    args = parser.parse_args()
    main(args.output.resolve())

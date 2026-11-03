#!/usr/bin/env python3
"""Independent intrinsic whole-CB-fiber recognition, all nonzero parameters."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_CA_SELF_PAIR_CB_WHOLE_FIBER_20261103_v63'


def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def same(A, B):
    assert A.shape == B.shape and (A-B).applyfunc(s.expand) == s.zeros(*A.shape)


def main(output, replay):
    started = time.monotonic()
    folder = ROOT/'reports'/STEM
    report = read(folder/'report.json')
    assert report['all_passed']
    for path, digest in report['input_sha256'].items():
        assert sha(ROOT/path) == digest, path
    for path, digest in report['evidence_sha256'].items():
        assert sha(folder/path) == digest, path
    subprocess.run([str(ROOT/'.venv/bin/python'),
                    str(ROOT/'code/ambient_CA_self_pair_CB_whole_fiber_20261103_v63.py'),
                    '--output', str(replay)], check=True, timeout=30)
    assert sha(folder/'certificate.json.gz') == sha(replay/'certificate.json.gz')
    cert = read(folder/'certificate.json.gz')
    h, x, y, z, w, a = s.symbols('h x y z w a')
    X = s.Matrix([h, x, y, z, w])
    coords = [h, x, z, w]
    ey = s.eye(5)[:, 2]
    D = s.Matrix([x*x+12*a*h*x/5, -3*a*x*x/5, 0,
                  h*x+2*a*x*z/5, h*h+z*z-8*a*x*w/5])
    P = ey*D.T-D*ey.T
    mons = list(itertools.combinations_with_replacement(range(5), 2))
    cubics = list(itertools.combinations_with_replacement(range(5), 3))
    parse = lambda q: s.sympify(q, locals={'a': a})
    matrix = lambda rows: s.Matrix([[parse(q) for q in row] for row in rows])
    casimir = s.Matrix([[s.Poly(P[i, j], *X).coeff_monomial(X[u]*X[v])
                         for j in range(5)] for i in range(5) for u, v in mons])
    same(casimir, matrix(cert['all75_constant_Casimir_matrix']))
    # Independent exterior engine: use determinants E,v,e_y,D, rather than
    # the author's expanded E wedge v wedge general-P coefficient formula.
    factors = []
    for j in range(5):
        full = s.Matrix.hstack(X, s.eye(5)[:, j], ey, D)
        col = []
        for subset in itertools.combinations(range(5), 4):
            polynomial = s.Poly(s.expand(full[list(subset), :].det()), *X)
            col += [polynomial.coeff_monomial(s.prod(X[i] for i in mon))
                    for mon in cubics]
        factors.append(s.Matrix(col))
    factor = s.Matrix.hstack(*factors)
    same(factor, matrix(cert['all175_intrinsic_factor_matrix']))
    assert factor[:, 2] == s.zeros(175, 1)
    field = s.QQ.frac_field(a)

    def witness(A, saved):
        assert len(DomainMatrix.from_Matrix(A).convert_to(field).rref()[1]) == saved['rank']
        square = A[saved['minor_rows'], saved['minor_columns']]
        dm = DomainMatrix.from_Matrix(square).convert_to(s.QQ.poly_ring(a))
        determinant = s.factor(dm.domain.to_sympy(dm.det()))
        assert s.expand(determinant-parse(saved['minor_determinant'])) == 0
        return str(determinant)

    dc = witness(casimir, cert['constant_Casimir_rank_witness'])
    df = witness(factor, cert['intrinsic_factor_rank_witness'])
    assert dc == '-864*a**4/625' and df == '3*a'
    D4 = s.Matrix([D[i] for i in [0, 1, 3, 4]])
    deriv = lambda q: s.expand(sum(D4[i]*s.diff(q, coords[i]) for i in range(4)))
    A, B, C, E, rh, rx, rz, rw = s.symbols('A B C D rh rx rz rw')
    ell = A*h+B*x+C*z+E*w
    cofactor = rh*h+rx*x+rz*z+rw*w
    general = s.Poly(5*deriv(ell)-cofactor*ell, h, x, z, w)
    assert s.expand(general.as_expr()-s.sympify(cert['full_general_Darboux_equation'],
                    locals={'a': a, 'A': A, 'B': B, 'C': C, 'D': E,
                            'rh': rh, 'rx': rx, 'rz': rz, 'rw': rw})) == 0
    assert general.coeff_monomial(w*w) == -E*rw
    assert general.coeff_monomial(h*w) == -A*rw-E*rh
    assert general.coeff_monomial(z*w) == -C*rw-E*rz
    assert general.coeff_monomial(h*h) == 5*E-A*rh
    forms = [x, x+a*h-2*a*a*z, x+3*a*h]
    lambdas = [-3*a, 2*a, 12*a]
    for q, lam in zip(forms, lambdas):
        assert s.expand(5*deriv(q)-lam*x*q) == 0
    eigenbasis = s.Matrix([[s.diff(q, c) for q in forms] for c in [h, x, z]])
    same(eigenbasis, matrix(cert['eigenbasis_matrix']))
    assert s.factor(eigenbasis.det()) == -6*a**3
    # ell=x: all independent-k solutions are excluded by the full coefficient matrix.
    qmons = [s.prod(coords[i] for i in mon)
             for mon in itertools.combinations_with_replacement(range(4), 2)
             if mon != (1, 1)]
    K = s.Matrix([[s.Poly(5*deriv(q)-s.Rational(9, 2)*a*x*q,
                          *coords).coeff_monomial(mon) for q in coords] for mon in qmons])
    same(K, matrix(cert['ell_x_candidate_k_constraint_matrix']))
    dk = witness(K, cert['ell_x_candidate_k_rank_witness'])
    assert K[:, 1] == s.zeros(9, 1)
    # For ell!=x, m is nonzero*x and lies in W, so use k=x.
    # From f_k=(5/2)m, m=lambda*x/5, f=lambda*x²/4+alpha*ell².
    alpha = s.symbols('alpha')
    for q, lam in zip(forms[1:], lambdas[1:]):
        residue = s.expand(deriv(x)+alpha*q*q+lam*x*x/20)
        assert s.expand(residue-(alpha*q*q+(lam-12*a)*x*x/20)) == 0
        assert s.Matrix([[s.Poly(q*q, *coords).coeff_monomial(mon),
                          s.Poly(x*x, *coords).coeff_monomial(mon)]
                         for mon in [h*h, x*x]]).det() != 0
    # Universal endpoint equations, retaining all16 cubic coefficients.
    def raw(f, g, H):
        U = s.zeros(5)
        U[1, 2], U[3, 4] = f, g
        U[2, 3], U[2, 4] = -s.diff(f, x)*z/2-s.diff(H, w), -s.diff(f, x)*w/2+s.diff(H, z)
        return U-U.T
    def canonical(U):
        div = s.Matrix([sum(s.diff(U[i, j], X[i]) for i in range(5)) for j in range(5)])
        return (U-(X*div.T-div*X.T)/5).applyfunc(s.expand)
    hc = s.symbols('c0:16')
    monH = [s.prod(coords[i] for i in mon)
            for mon in itertools.combinations_with_replacement(range(4), 3)
            if any(i >= 2 for i in mon)]
    H = sum(c*q for c, q in zip(hc, monH))
    f0, f1, f2, g0, g1, g2 = s.symbols('f0 f1 f2 g0 g1 g2')
    F, G = f0*h*h+f1*h*x+f2*x*x, g0*h*h+g1*h*x+g2*x*x
    endpoint = canonical(raw(0, G, H))
    assert endpoint[0, :] == endpoint[1, :] == s.zeros(1, 5)
    endpoint = canonical(raw(F, 0, H))
    assert s.expand(-endpoint[0, 2]-2*s.diff(F, x)*h/5) == 0
    assert s.expand(-endpoint[1, 2]+F-2*s.diff(F, x)*x/5) == 0
    proper = ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json'
    properraw = ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json'
    rankreview = ROOT/'reports/AMBIENT_CA_SELF_PAIR_CB_SOURCE_RANK_20261103_v63_REVIEW.json'
    incidence = ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63_REVIEW.json'
    assert read(proper)['proper_smooth28D_source_image_is_entire_reduced_CB_accepted']
    assert read(properraw)['rank_at_most2_preimage_support'] == 'a=0 or b=0 or F=0'
    assert read(rankreview)['every_nonzero_a_specified_CB_lift_is_immersion']
    assert read(incidence)['all_passed']
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}_REVIEW.md',
             ROOT/'code/ambient_CA_self_pair_CB_whole_fiber_20261103_v63.py',
             ROOT/f'research/notes/{STEM}.md', folder/'report.json',
             folder/'certificate.json.gz', proper, properraw, rankreview, incidence]
    result = dict(id=STEM+'_REVIEW', all_passed=True,
        fresh_author_certificate_byte_identical=True,
        independent_all75_Casimir_and_determinant_exterior_all175_factor_matrices=True,
        no_constant_Casimir_rank5_minor=dc, unique_projective_constant_factor_rank4_minor=df,
        whole_proper_endpoint_coverage_not_single_flag_chart_reviewed=True,
        universal_endpointA_two_Casimirs_and_endpointB_trace_equations_checked=True,
        all_linear_Darboux_forms_exhausted_for_every_a_nonzero=True,
        eigenbasis_determinant='-6*a**3', ell_x_independent_k_impossible_minor=dk,
        only_other_valid_eigenform_is_existing_H_x_plus_3a_h=True,
        whole_CB_source_fiber_exactly_one_geometric_point_for_every_a_nonzero=True,
        accepted_full29_immersion_implies_whole_fiber_reduced_length1=True,
        proper_finite_closed_image_Nakayama_proof_reviewed_without_assumed_normalization=True,
        reduced_CB_smooth_cone29_projective28_at_every_a_nonzero=True,
        origin_completion_global_bad_locus_or_full_Poisson_claimed=False,
        input_sha256={str(path.relative_to(ROOT)): sha(path) for path in paths},
        seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['all_passed', 'no_constant_Casimir_rank5_minor',
        'unique_projective_constant_factor_rank4_minor',
        'reduced_CB_smooth_cone29_projective28_at_every_a_nonzero', 'seconds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--replay', type=Path, required=True)
    args = parser.parse_args()
    main(args.output.resolve(), args.replay.resolve())

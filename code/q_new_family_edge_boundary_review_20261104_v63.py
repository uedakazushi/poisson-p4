#!/usr/bin/env python3
"""Independent normal-form edge algebra and source B-valuation, no source32 products."""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
import math
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_NEW_FAMILY_EDGE_BOUNDARY_VALUATION_20261104_v63_REVIEW'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for a in iter(lambda: f.read(1 << 20), b''):
            h.update(a)
    return h.hexdigest()


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2000000000, 2000000000))
    out.mkdir(parents=True, exist_ok=False)
    edge = ROOT/'reports/Q_FIRSTJET_CROSS_BLOCK_EDGE_BOUNDARY_20261104_v63'
    source = ROOT/'reports/Q_NEW_FAMILY_EDGE_BOUNDARY_VALUATION_20261104_v63'
    er, ec = [read(edge/p) for p in ['report.json', 'certificate.json']]
    sr, sc = [read(source/p) for p in ['report.json', 'certificate.json.gz']]
    gates = [ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63_REVIEW/report.json',
        ROOT/'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63_REVIEW/report.json',
        ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_NORMAL_FORM_20261103_v63_REVIEW/report.json',
        ROOT/'reports/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63_REVIEW/report.json']
    assert er['all_passed'] and sr['all_passed'] and all(read(p)['all_passed'] for p in gates)
    x, y, a, g, h, i, j = s.symbols('x y a g h i j')
    S, delta = x*x+y*y, x*x-y*y
    D0, p4 = delta**2, x*x*y*y
    m = g*j-h*i
    ell = (a+x)**2-y*y
    d = -(a*x+S)/y
    cc, ff = m-ell, m-x*x*ell/y**2
    C = s.Matrix([[a, 1, g, h], [cc, -a, i, j], [-j, h, d, 1], [i, -g, ff, -d]])
    D = s.diag(x, -x, y, -y)
    qC, qD, qCD = s.trace(C*C), s.trace(D*D), s.trace(C*D)
    assert s.cancel(qC-2*qD) == 0 and s.cancel(qCD+qD) == 0
    # Newton coefficient formulas are an independent entrance into the
    # actual normal-form pencil, rather than reading edge_u/v substitutions.
    actual_p3 = s.cancel(qD*qCD/2-s.trace(C*D**3))
    p3 = -2*a*x*delta-2*x*x*S
    assert s.cancel(actual_p3-p3) == 0
    actual_p2 = s.cancel(qC*qD/4+qCD**2/2-s.trace(C*C*D*D)-s.trace(C*D*C*D)/2)
    B2 = -4*a*a*x*x-2*a*x**3-6*a*x*y*y+3*x**4+y**4
    kap = 2*x*y*(g*j+h*i)+S*m
    assert s.cancel(actual_p2-B2-kap) == 0
    actual_N = s.trace(C*D**2*C*D**2)-S*s.trace(C*C*D*D)+p4*s.trace(C*C)
    assert s.cancel(actual_N-2*D0*m) == 0
    A = p3+S*S
    assert s.expand(A+delta*(2*a*x+S)) == 0
    assert s.cancel(B2-(S*S-p3+D0-A*A/D0)) == 0
    K16 = D0*(actual_p2-S*S+p3-D0)+A*A
    assert s.cancel(K16-D0*kap) == 0
    uv = (kap*kap-2*S*kap*m+D0*m*m)/(16*p4)
    assert s.cancel(uv-g*j*h*i) == 0
    G32 = D0*(2*D0*m)**2-4*S*K16*(2*D0*m)+4*K16*K16
    assert s.cancel(G32-64*p4*D0**2*g*h*i*j) == 0
    F4, p2, p3v, p4v, N12 = s.symbols('F4 p2 p3 p4 N12')
    ss = F4/20
    dd = ss**2-4*p4v
    kk = dd*(p2-ss**2+p3v-dd)+(p3v+ss**2)**2
    abstract = s.Poly(dd*N12**2-4*ss*kk*N12+4*kk**2, F4, p2, p3v, p4v, N12)
    assert len(abstract.terms()) == 36
    assert all(sum(e*w for e, w in zip(exps, [4, 8, 8, 8, 12])) == 32 for exps, c in abstract.terms())
    assert s.expand(abstract.as_expr().coeff(N12, 2)-dd) == 0
    denominator = s.ilcm(1, *(s.denom(c) for c in abstract.coeffs()))
    content = math.gcd(*(int(c*denominator) for c in abstract.coeffs()))
    multiplier = s.Rational(denominator, content)
    primitive = s.Poly(abstract.as_expr()*multiplier, *abstract.gens)
    assert primitive.terms() == [(tuple(e), s.Rational(c)) for e, c in ec['primitive_integer_terms']]
    # Legacy metadata wrote this INTEGER with a decimal .0; parse it exactly.
    assert multiplier == s.Rational(ec['primitive_scale_over_displayed_G32'])
    expansion = ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63'
    names = ['F4', 'p0', 'p2', 'p3', 'N12']
    polynomials = {n: read(expansion/(n+'.json.gz')) for n in names}
    info = {}
    for n, p in polynomials.items():
        assert int(p['denominator_integer']) > 0
        info[n] = (max(t[4] for t in p['terms']), p['denominator_B_power'], p['denominator_D_power'])
    assert info == dict(F4=(0, 2, 4), p0=(0, 3, 8), p2=(1, 4, 8), p3=(1, 4, 8), N12=(1, 5, 12))
    def bound(factors, targetB, targetH):
        numerator = sum(info[n][0]*e for n, e in factors.items())
        db = sum(info[n][1]*e for n, e in factors.items())
        dh = sum(info[n][2]*e for n, e in factors.items())
        assert db <= targetB and dh <= targetH
        return numerator+targetB-db
    kparts = [{'p0': 1, 'p2': 1}, {'p0': 1, 'F4': 2}, {'p0': 1, 'p3': 1}, {'p0': 2},
        {'p3': 2}, {'p3': 1, 'F4': 2}, {'F4': 4}]
    k_bounds = [bound(f, 8, 16) for f in kparts]
    assert k_bounds == [2, 1, 2, 2, 2, 1, 0]
    # Treat the already-cleared K numerator as an unknown polynomial degree<=2;
    # no actual source polynomial products or expansion of G32 are performed.
    g_bounds = [info['p0'][0]+2*info['N12'][0]+16-(info['p0'][1]+2*info['N12'][1]),
        info['F4'][0]+2+info['N12'][0]+16-(info['F4'][1]+8+info['N12'][1]), 4]
    assert g_bounds == [5, 4, 4]
    P0 = polynomials['p0']['terms']
    leadN = [t for t in polynomials['N12']['terms'] if t[4] == 1]
    assert all(t[4] == 0 and int(t[0]) for t in P0) and P0 and leadN
    assert all(int(t[0]) for t in leadN)
    assert P0 == sc['leading_coefficient_factors']['p0_numerator_terms']
    assert leadN == sc['leading_coefficient_factors']['N12_leading_B1_coefficient_terms']
    # Their product is nonzero in QQ(a,Y,W); this follows from the domain
    # property without computing that product. Its B^5 cannot cancel lower B.
    evidence = dict(direct_generic_matrix_p2=str(actual_p2), direct_generic_matrix_p3=str(actual_p3),
        direct_generic_matrix_N12=str(s.factor(actual_N)), cofactor=64,
        primitive_G32_terms=len(primitive.terms()), primitive_exact_multiplier=str(multiplier),
        source_B_metadata={k: list(v) for k, v in info.items()}, K_cleared_B_bounds=k_bounds,
        G_cleared_B_bounds=g_bounds, p0_nonzero_leading_factor_terms=len(P0), N12_nonzero_leading_factor_terms=len(leadN),
        no_source_polynomial_product_or_source32_expansion=True,
        generic_NEW_edge_open='g*h*i*j !=0 on the accepted xy*delta*b*e chart',
        L_unit_or_NEW_multisection_degree_or_equation_not_claimed=True)
    (out/'independent-certificate.json').write_text(json.dumps(evidence, indent=2)+'\n')
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md',
        *[edge/p for p in ['report.json', 'certificate.json']],
        *[source/p for p in ['report.json', 'certificate.json.gz']],
        ROOT/'code/q_firstjet_cross_block_edge_boundary_20261104_v63.py',
        ROOT/'code/q_new_family_edge_boundary_valuation_20261104_v63.py',
        ROOT/'research/notes/Q_FIRSTJET_CROSS_BLOCK_EDGE_BOUNDARY_20261104_v63.md',
        ROOT/'research/notes/Q_NEW_FAMILY_EDGE_BOUNDARY_VALUATION_20261104_v63.md',
        *gates, *[expansion/(n+'.json.gz') for n in names]]
    result = dict(id=STEM, all_passed=True,
        independent_generic_matrix_Newton_p2_p3_N12_and_edge_cofactor64=True,
        all36_abstract_coefficients_and_exact_primitive_scale_verified=True,
        full_low_source_B_degrees_and_unique_noncancelling_B5_leading_factor_verified=True,
        G32_not_identically_zero_on_actual_NEW25_source=True,
        spectral_dominance_and_known_nonzero_gauge_imply_generic_Laurent_chart_nonempty=True,
        no_source32_product_expansion_new_point_or_high_degree_probe=True,
        L_unit_NEW_multisection_degree_genus_or_new_equation_claimed=False,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-certificate.json': sha(out/'independent-certificate.json')},
        seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'G32_not_identically_zero_on_actual_NEW25_source', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())

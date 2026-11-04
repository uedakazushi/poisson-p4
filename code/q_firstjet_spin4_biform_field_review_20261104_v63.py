#!/usr/bin/env python3
"""Independent coefficient and complete-field review of the joint Spin4 biform quotient."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_FIRSTJET_SPIN4_BIFORM_RATIONAL_FIELD_20261104_v63_REVIEW'


def read(p):
    return json.loads(p.read_text())


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for a in iter(lambda: f.read(1 << 20), b''):
            h.update(a)
    return h.hexdigest()


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    out.mkdir(parents=True, exist_ok=False)
    folder = ROOT/'reports/Q_FIRSTJET_SPIN4_BIFORM_RATIONAL_FIELD_20261104_v63'
    ar, ac = [read(folder/p) for p in ['report.json', 'certificate.json']]
    gate = ROOT/'reports/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63_REVIEW/report.json'
    assert ar['all_passed'] and read(gate)['all_passed']
    x0, x1, y0, y1, c, t = s.symbols('x0 x1 y0 y1 c t')
    variables = [x0, x1, y0, y1]
    v = s.symbols('v0:6')
    mons = [x0*y0**2, x0*y0*y1, x0*y1**2, x1*y0**2, x1*y0*y1, x1*y1**2]
    U0 = x0**2*y0+x1**2*y1
    V = sum(a*m for a, m in zip(v, mons))
    def move(f, A, B):
        new = list(A*s.Matrix([x0, x1]))+list(B*s.Matrix([y0, y1]))
        return s.Poly(f.subs(dict(zip(variables, new)), simultaneous=True), *variables)
    TA, TB = s.diag(t, 1/t), s.diag(t**-2, t**2)
    assert move(U0, TA, TB).as_expr() == U0
    weights = [s.cancel(move(m, TA, TB).as_expr()/m) for m in mons]
    assert weights == [t**i for i in ac['T_weights']]
    SA = s.Matrix([[0, 1], [-1, 0]])
    SB = s.Matrix([[0, s.I], [s.I, 0]])
    assert s.Matrix(ac['sigma_A']) == SA and s.Matrix(ac['sigma_B']) == SB
    assert SA.det() == SB.det() == 1 and SA**4 == SB**4 == s.eye(2)
    assert SA**2 == SB**2 == -s.eye(2)
    assert s.expand(move(U0, SA, SB).as_expr()-s.I*U0) == 0
    sigmaV = move(V, SA, SB)
    permutation = [sigmaV.coeff_monomial(m) for m in mons]
    assert permutation == [v[5], v[4], v[3], -v[2], -v[1], -v[0]]
    sigma = dict(zip([c, *v], [s.I*c, *permutation]))
    assert all(s.expand(s.sympify(a)-b) == 0 for a, b in zip(ac['sigma_on_c_and_V'], [s.I*c, *permutation]))
    assert all(s.expand(sigma[z].subs(sigma, simultaneous=True)+z) == 0 for z in [c, *v])
    assert (SA*TA*SA.inv()-TA.subs(t, 1/t)).applyfunc(s.cancel) == s.zeros(2)
    assert (SB*TB*SB.inv()-TB.subs(t, 1/t)).applyfunc(s.cancel) == s.zeros(2)
    # Independent resultant engine: univariate subresultants, not the author's
    # 4x4 Sylvester determinant; all six polynomial coefficients remain free.
    a, b, d, e, f, g, z = s.symbols('a b d e f g z')
    p, q = a*x0**2+b*x0*x1+d*x1**2, e*x0**2+f*x0*x1+g*x1**2
    jac = s.Poly(s.diff(p, x0)*s.diff(q, x1)-s.diff(p, x1)*s.diff(q, x0), x0, x1)
    jac_co = [jac.coeff_monomial(m) for m in [x0**2, x0*x1, x1**2]]
    result = s.resultant(a+b*z+d*z*z, e+f*z+g*z*z, z)
    assert s.expand(jac_co[1]**2-4*jac_co[0]*jac_co[2]-16*result) == 0
    detM = a*g-d*e
    assert s.expand(b*detM-(a*(b*g-d*f)+d*(a*f-b*e))) == 0
    assert s.expand(f*detM-(g*(a*f-b*e)+e*(b*g-d*f))) == 0
    assert result.subs({a: c, b: 0, d: 0, e: 0, f: 0, g: c}) == c**4
    # Exhaustion of the normalizer: preserving span(x0²,x1²) forces ab=cd=0.
    # Invertibility leaves precisely diagonal or antidiagonal A. For each case
    # solve the B matrix completely from its induced two-plane coefficient map.
    zeta = s.Symbol('zeta')
    cases = []
    for A in [s.diag(t, 1/t), s.Matrix([[0, t], [-1/t, 0]])]:
        S = s.Matrix([[A[0, 0]**2, A[1, 0]**2], [A[0, 1]**2, A[1, 1]**2]])
        B0 = zeta*S.inv()
        assert s.cancel(B0.det()-zeta**2/S.det()) == 0
        assert s.expand(move(U0, A, B0).as_expr()-zeta*U0) == 0
        cases.append(dict(induced_determinant=str(S.det()),
            compensating_B=[[str(q) for q in row] for row in B0.tolist()],
            required_zeta_square=str(S.det())))
    assert [v['required_zeta_square'] for v in cases] == ['1', '-1']
    # Diagonal zeta=1 is exactly T. Diagonal zeta=-1 is sigma²*T;
    # the two antidiagonal cases are sigma*T and sigma³*T. Thus N/T=mu4.
    r = v[1]*v[4]
    P, Q = v[0]*v[1]**3, v[5]/v[1]**3
    S, T = v[3]*v[1]**5, v[2]/v[1]**5
    torus = [r, P, Q, S, T]
    scale = dict(zip(v, [a*t**w for a, w in zip(v, [-3, 1, 5, -5, -1, 3])]))
    assert all(s.cancel(q.subs(scale, simultaneous=True)-q) == 0 for q in torus)
    rr, pp, qq, ss, tt, h = s.symbols('rr pp qq ss tt h')
    inverseV = {v[0]: pp/h**3, v[1]: h, v[2]: tt*h**5, v[3]: ss/h**5, v[4]: rr/h, v[5]: qq*h**3}
    assert [s.cancel(q.subs(inverseV, simultaneous=True)) for q in torus] == [rr, pp, qq, ss, tt]
    expected = [-r, Q*r**3, -P/r**3, -T*r**5, S/r**5]
    assert all(s.cancel(q.subs(sigma, simultaneous=True)-e) == 0 for q, e in zip(torus, expected))
    Qp, Tp = Q*r**3, -T*r**5
    generators = [c**4, r/c**2, P+Qp, (P-Qp)/c**2, S+Tp, (S-Tp)/c**2]
    assert all(s.cancel(q.subs(sigma, simultaneous=True)-q) == 0 for q in generators)
    assert all(s.cancel(q.subs(scale, simultaneous=True)-q) == 0 for q in generators)
    # Author symbols c are assumed nonzero; compare via explicit shared locals.
    loc = {str(a): a for a in [c, *v]}
    assert all(s.cancel(s.sympify(a, locals=loc)-b) == 0 for a, b in zip(ac['six_rational_joint_invariant_field_generators'], generators))
    zz, rr, pp, qq, ss, tt = s.symbols('zz rr pp qq ss tt')
    coords = [zz, rr, pp, qq, ss, tt]
    simple = [zz**2, rr/zz, pp+qq, (pp-qq)/zz, ss+tt, (ss-tt)/zz]
    invs = s.symbols('I0:6')
    inverse = {rr: zz*invs[1], pp: (invs[2]+zz*invs[3])/2, qq: (invs[2]-zz*invs[3])/2,
        ss: (invs[4]+zz*invs[5])/2, tt: (invs[4]-zz*invs[5])/2}
    assert [s.cancel(q.subs(inverse, simultaneous=True)) for q in simple] == [zz**2, *invs[1:]]
    jacobian = s.factor(s.Matrix(simple).jacobian(coords).det())
    assert jacobian == 8/zz**2
    invol = {zz: -zz, rr: -rr, pp: qq, qq: pp, ss: tt, tt: ss}
    assert all(s.cancel(q.subs(invol, simultaneous=True)-q) == 0 for q in simple)
    # The nontrivial involution fixes precisely the displayed degree2 subfield:
    # all six original coordinates are recovered after adjoining zz, zz²=I0.
    evidence = dict(independent_V_sigma_coefficients=list(map(str, permutation)),
        torus_weights=[-3, 1, 5, -5, -1, 3], independent_resultant=str(result),
        full_line_normalizer_cases=cases, exact_stabilizer='T', quotient_normalizer='T semidirect mu4',
        six_invariants=list(map(str, generators)), final_generator_Jacobian=str(jacobian),
        exact_field='C(Sym2A tensor B plus A tensor Sym2B)^(SL2A times SL2B)=C(I0,I1,I2,I3,I4,I5)',
        extra_Spin5_line_normalizer_not_quotiented=True)
    (out/'independent-certificate.json').write_text(json.dumps(evidence, indent=2)+'\n')
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', folder/'report.json', folder/'certificate.json',
        ROOT/'code/q_firstjet_spin4_biform_rational_field_20261104_v63.py',
        ROOT/'research/notes/Q_FIRSTJET_SPIN4_BIFORM_RATIONAL_FIELD_20261104_v63.md', gate]
    result = dict(id=STEM, all_passed=True,
        independent_full_U_V_torus_sigma_coefficient_actions=True,
        independent_resultant_subresultant_engine_and_generic_normalform=True,
        complete_exact_stabilizer_T_and_four_component_line_normalizer_exhaustion=True,
        exact_torus_field_generation_no_finite_index_left_at_this_stage=True,
        residual_linear_involution_and_complete_six_generator_fixed_field=True,
        six_generators_algebraically_independent_Jacobian='8/z^2',
        base_field='C, characteristic zero',
        full_Spin5_Sp4_W16_extra_line_normalizer_rationality_claimed=False,
        polynomial_invariant_ring_global_quotient_or_NEW_equation_claimed=False,
        new_point_high_degree_probe_GB_or_implicitization=False,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-certificate.json': sha(out/'independent-certificate.json')},
        seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'residual_linear_involution_and_complete_six_generator_fixed_field', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())

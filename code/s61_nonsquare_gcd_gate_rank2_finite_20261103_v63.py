#!/usr/bin/env python3
"""Exact finite certificate for the rank-two branch of the B,W gcd gate.

No Groebner basis, field reduction, or generic-point inference is used.
The source is the frozen cleared four-quadric coefficient report.
"""
from pathlib import Path
import hashlib
import json
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'work/nonsquare-v62/gcd-gate-exception-v54'

def main():
    started = time.monotonic()
    OUT.mkdir(parents=True, exist_ok=True)
    source = ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json'
    raw = json.loads(source.read_text())
    assert raw['all_passed']
    assert hashlib.sha256(source.read_bytes()).hexdigest() == '80f7c46b1bd4ffef4a78f492ad4f7a2bb40b2ffa27cddcd23d5b5d660db4a7e0'
    a, e, P = s.symbols('a e P')
    x = s.symbols('x0:6')
    z = s.symbols('z0:4')
    quads = []
    for i, row in enumerate(raw['cleared_qfour_sparse']):
        q = 0
        for mon, terms in row:
            c = sum(s.Rational(co)*a**(ex[0]-ex[2]+[0,2,2,1][i])
                    *e**ex[1]*P**ex[2] for ex, co in terms)
            q += c*s.prod(v**k for v, k in zip(x, mon[:6]))
        quads.append(q)
    Rc = 15*P**2-6*P*e+6*P+80*a
    Tc = -5*P+6*e-2
    Sc = 21*P*e+16*P+140*a
    Vc = 9*P*e+4*P+60*a
    Kc = 85*P**2-12*P*e+58*P+600*a
    J = 5*P-6*e-4
    B = 30*a*Rc*x[0]-9*P*Sc*x[1]+90*P*a*Tc*x[2]-18*P*Vc*x[3]-2*a*Kc*x[4]
    A = (B+10*a*P*J*x[4])/3
    Q2 = s.cancel((quads[2]-e*quads[1])/(a*a*J))
    assert s.expand(quads[1].subs(x[5],0)-A*B/6075) == 0
    assert s.expand(Q2.subs(x[5],0)-40*x[0]*B/6075) == 0
    assert s.expand(3*Sc-7*Vc-20*P) == 0
    t = (30*a*Rc*z[0]+90*P*a*Tc*z[1]-18*P*Vc*z[3]-2*a*Kc*z[2])/(180*P**2)
    cut = {x[0]:z[0], x[1]:3*t, x[2]:z[1], x[3]:-s.Rational(7,2)*t+z[3], x[4]:z[2], x[5]:0}
    assert s.cancel(B.subs(cut, simultaneous=True)) == 0
    G = s.Poly(s.cancel(P**4*quads[3].subs(cut, simultaneous=True)), *z)
    for c in G.coeffs():
        s.Poly(c, a, e, P, domain=s.QQ)  # no hidden parameter denominator
    r, U, V, T, R = s.symbols('r U V T R')
    inverse = {x[0]:r, x[1]:V,
               x[2]:R/a-U-6*(e+1)*V/(5*a)-(9*e+4)*T/(10*a),
               x[3]:(T-V)/2,
               x[4]:5*U+5*r-4*R/a+24*(e+1)*V/(5*a)+2*(9*e+4)*T/(5*a), x[5]:0}
    assert s.cancel(4*quads[0].subs(inverse, simultaneous=True)-U*R+5*V*T) == 0
    Bnormal = s.Poly(s.cancel(B.subs(inverse, simultaneous=True)), r,U,V,T,R)
    coeff = [Bnormal.coeff_monomial(v) for v in (r,U,V,T,R)]
    DB = 20*P**2+3*P*e+20*P+180*a
    assert s.expand(coeff[0]+20*a*DB) == 0
    eliminate_a = -P*(20*P+3*e+20)/180
    dual_norm = coeff[1]*coeff[4]-coeff[2]*coeff[3]/5
    N = s.Poly(s.cancel(-s.Rational(10125,4)*dual_norm.subs(a,eliminate_a)/P**2), P,e, domain=s.QQ)
    Ne = s.Poly(N.as_expr(),e)
    assert Ne.degree() == 4 and Ne.LC() == 72442188
    Gdb = s.Poly(s.cancel(G.as_expr().subs(a,eliminate_a)),*z)
    for c in Gdb.coeffs():
        s.Poly(c,P,e,domain=s.QQ)
    H = s.hessian(Gdb.as_expr(),z)/2
    for entry in H:
        s.Poly(entry,P,e,domain=s.QQ)
    data = dict(H3_on_DB=[[s.Poly(H[i,j],P,e).terms() for j in range(4)] for i in range(4)])
    # Only a univariate specialization is expanded; monicity supplies the finite-family argument.
    indices = [0,1,2]
    H1 = H.subs(P,1).extract(indices,indices)
    M1 = s.Poly(s.cancel(H1.det()),e,domain=s.QQ)
    N1 = s.Poly(N.as_expr().subs(P,1),e,domain=s.QQ)
    left,right,gcd = s.gcdex(N1,M1)
    assert gcd == s.Poly(1,e,domain=s.QQ)
    assert s.expand(left.as_expr()*N1.as_expr()+right.as_expr()*M1.as_expr()-1) == 0
    encode = lambda q: [[list(ex),str(c)] for ex,c in q.terms()]
    certificate = dict(N=encode(N),N1=encode(N1),M1=encode(M1),Bezout_left=encode(left),Bezout_right=encode(right),
                       global_H3_entries=[[encode(s.Poly(H[i,j],P,e)) for j in range(4)] for i in range(4)])
    (OUT/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    result = dict(all_exact_checks_passed=True, source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  coefficient_variables=['P','e'], DB=str(DB),a_elimination=str(eliminate_a),N=str(N.as_expr()),
                  N_total_degree=N.total_degree(),N_e_degree=4,N_e_leading_QQ_constant='72442188',
                  H3_entries_verified_polynomial_over_QQ_P_e=True, principal_minor_indices=indices,
                  exact_specialization_P=1, specialized_N_degree=N1.degree(),specialized_minor_degree=M1.degree(),
                  exact_QQ_univariate_Bezout_identity=True,finite_rank2_gcd_failure_branch_proved=True,
                  high_rank_proportional_branch_not_handled=True,full_gcd_gate_not_proved=True,
                  five_quadric_rank_gate_not_proved=True,all_core_rank24_or_product41_not_proved=True,
                  seconds=time.monotonic()-started)
    (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact smooth-quadric isotropic Log boundary extensions; Python 3 + SymPy 1.14.

Produces a small deterministic certificate. Does not use previous producers or
numerical sampling, and does not certify the 147-component classification.
"""
import argparse
import itertools
import json
from pathlib import Path
import sympy as s

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, default=Path('work/q3-isotropic-log.json'))
args = parser.parse_args()
assert s.__version__ == '1.14.0', s.__version__
X = s.symbols('u0 u1 z w0 w1')
u0, u1, z, w0, w1 = X
E = s.Matrix(X)
PAIRS = list(itertools.combinations(range(5), 2))
TRIPLES = list(itertools.combinations(range(5), 3))
xi = (u0, u1, -2*z)

def sign(word):
    return (-1)**sum(word[i] > word[j] for i in range(len(word)) for j in range(i+1, len(word)))

def log_lift(q):
    # ι_P Ω = dq ∧ ι_ξ(du0∧du1∧dz)
    gamma = {(1,2):xi[0], (0,2):-xi[1], (0,1):xi[2]}
    P = s.zeros(5)
    for i,j in PAIRS:
        rest = tuple(k for k in range(5) if k not in (i,j))
        value = 0
        for r,k in enumerate(rest):
            others = tuple(rest[t] for t in range(3) if t != r)
            value += (-1)**r*s.diff(q, X[k])*gamma.get(others, 0)
        P[i,j] = s.expand(sign((i,j,*rest))*value)
        P[j,i] = -P[i,j]
    return P

def divergence(P):
    return s.Matrix([s.expand(sum(s.diff(P[i,j], X[i]) for i in range(5))) for j in range(5)])

def jacobi(P):
    return {ijk:s.expand(sum(P[i,k]*s.diff(P[j,l],X[k]) +
                       P[j,k]*s.diff(P[l,i],X[k]) +
                       P[l,k]*s.diff(P[i,j],X[k]) for k in range(5)))
            for ijk in TRIPLES for i,j,l in [ijk]}

def pf_beta(P):
    output = []
    for omit in range(5):
        i,j,k,l = [n for n in range(5) if n != omit]
        output.append(s.factor((-1)**omit*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])))
    return output

def check(q, name, expected_rank, expected_groebner, expected_f, expected_wrank):
    G = s.hessian(q,E)/2
    assert G.det() != 0
    assert G.extract([3,4],[3,4]).rank() == expected_wrank
    P0 = log_lift(q)
    assert divergence(P0) == s.zeros(5,1)
    assert (P0*s.Matrix([s.diff(q,zeta) for zeta in X])).applyfunc(s.expand) == s.zeros(5,1)
    assert all(v == 0 for v in jacobi(P0).values())
    assert all(v == 0 for v in pf_beta(P0))
    zz = s.symbols('Z0:10')
    Z = s.zeros(5)
    for a,(i,j) in zip(zz,PAIRS): Z[i,j] = a;Z[j,i] = -a
    A = G.inv()*Z
    vec = A*E
    equations = []
    for i,j in PAIRS:
        expr = s.expand(sum(vec[k]*s.diff(P0[i,j],X[k]) for k in range(5)) -
                        sum(A[i,k]*P0[k,j]+P0[i,k]*A[j,k] for k in range(5)))
        equations.extend(s.Poly(expr,*X).coeffs())
    M = s.linear_eq_to_matrix(equations,zz)[0]
    assert M.rank() == expected_rank
    rr=list(M.T.rref()[1]); cc=list(M.rref()[1])
    exact_minor=s.factor(M.extract(rr,cc).det())
    assert exact_minor!=0 and len(rr)==expected_rank
    basis = M.nullspace()
    parameters = s.symbols('a b c')[:len(basis)]
    P = P0.copy()
    Kforms = []
    fields = []
    for p,b in zip(parameters,basis):
        sub = {zz[i]:b[i] for i in range(10)}
        Ai = A.subs(sub)
        Zi = Z.subs(sub)
        wi = (Ai*E).applyfunc(s.expand)
        # With v_T=T(dq), K_T=qT+1/2(E∧v_T) is strictly q-Casimir.
        Ti = G.inv()*Zi*G.inv()/2
        Ki = (q*Ti + (E*wi.T-wi*E.T)/2).applyfunc(s.expand)
        assert (Ki*s.Matrix([s.diff(q,zeta) for zeta in X])).applyfunc(s.expand) == s.zeros(5,1)
        assert divergence(Ki) == (s.Rational(3,2)*wi).applyfunc(s.expand)
        fields.append([str(a) for a in wi])
        Kforms.append(Ki)
        P += p*Ki
    coeffs = set()
    J = jacobi(P)
    for expression in J.values():
        for a in s.Poly(expression,*X).coeffs():
            if a: coeffs.add(s.factor(a))
    GB = s.groebner(list(coeffs),*parameters,order='lex')
    assert set(s.expand(poly.as_expr()) for poly in GB.polys) == set(s.expand(e) for e in expected_groebner)
    beta = pf_beta(P)
    f = s.factor(beta[0]/s.diff(q,u0))
    assert s.expand(f-expected_f) == 0, (f,expected_f)
    assert all(s.expand(beta[j]-f*s.diff(q,X[j])) == 0 for j in range(5))
    assert s.expand(sum(X[j]*beta[j] for j in range(5))-2*q*f) == 0
    print(name, 'qdet',G.det(),'Wrank',expected_wrank, 'stabilizer', len(basis), 'Jacobi ideal', GB)
    return {
        'name':name,
        'q':str(q),'determinant_of_quadric_matrix':str(G.det()),
        'rank_q_restriction_to_W':expected_wrank,
        'strict_stabilizer_matrix_shape':list(M.shape),
        'strict_stabilizer_matrix_rank':expected_rank,
        'rank_witness':{'rows':rr,'columns':cc,'minor':str(exact_minor)},
        'strict_stabilizer_fields':fields,
        'extension_parameters':list(map(str,parameters)),
        'all_Jacobi_polynomial_coefficients':len(coeffs),
        'exact_Jacobi_Groebner_basis':[str(z.as_expr()) for z in GB.polys],
        'pfaffian_cofactor_f':str(f),
        'rank4_reduced_equation':'f != 0',
        'scheme_reducedness':'nonreduced' if name == 'totally_isotropic_W' else 'reduced',
        'nonzero_P0_on_Q':any(s.expand(P0[i,j]) != 0 for i,j in PAIRS),
    }

p,q,r = s.symbols('a b c')
rows = [
  check(u0*u0+u1*u1+z*z+2*u0*w0+w1*w1,
        'rank_one_W_compatible_plane',8,[],s.Rational(3,2)*p*u0*u0*z,1),
  check(u0*u0+2*u1*w0+2*z*w1+w0*w0,
        'rank_one_W_transverse_lines',8,[p*(q-4)],
        p*z*(q*u0**2+2*q*u1*w0+q*w0**2+2*q*w1*z-4*u0**2+12*u1**2+16*u1*w0+8*w0**2+16*w1*z)/8,1),
  check(u0*u0+2*u1*w0+2*z*w1,
        'totally_isotropic_W',7,[p*r,r*r],
        -(-p*u1+r*u0)*(-r*(u0*u0+2*u1*w0+2*z*w1)+12*u1*z)/8,0),
]
assert rows[0]['strict_stabilizer_matrix_rank']==8
assert rows[1]['exact_Jacobi_Groebner_basis']==['a*b - 4*a']
assert rows[2]['exact_Jacobi_Groebner_basis']==['a*c','c**2']
certificate={
   'date':'2026-11-10','version':'v7', 'status':'PASS', 'sympy':s.__version__,
   'setup':'smooth Q in P4; resonant logarithmic alpha with xi=(1,1,-2) on (u0,u1,z); W=span(w0,w1)',
   'models':rows,
   'interpretation':'The totally-isotropic case has the nonreduced scheme Spec C[a,b,c]/(a*c,c*c), with reduced rank-four affine plane c=0, a != 0, and an embedded square-zero direction along a=c=0. Two rank-one W examples have different extension schemes; rank(q|W) alone does not determine the fibre.',
   'verification_scope':'Exact finite coefficient calculations with SymPy; general restriction-space theory is proved in the current manuscript, not by this certificate; classification of all q with degenerate restriction is NOT claimed.',
}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
print('SUCCESS',args.output)

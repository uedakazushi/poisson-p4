#!/usr/bin/env python3
"""Exact symbolic checks for the v18 multiplicity-three restrictions.

Uses sympy, no numerical sampling. This checks identities used in the proof;
normality and the geometric conclusions are proved in m3_restrictions.tex.
"""
import itertools
import json
from pathlib import Path

import sympy as S


def check_zero(value):
    assert S.expand(value) == 0, S.expand(value)


# Wedge with a nondegenerate bivector on a four dimensional vector space.
k01, k02, k03, k12, k13, k23 = ks = S.symbols("k01 k02 k03 k12 k13 k23")
vec = S.Matrix(S.symbols("v0:4"))
K = S.zeros(4)
for pair, val in zip(itertools.combinations(range(4), 2), ks):
    i, j = pair
    K[i, j] = val
    K[j, i] = -val
triples = list(itertools.combinations(range(4), 3))
wedge = S.Matrix([
    K[i, j]*vec[k] - K[i, k]*vec[j] + K[j, k]*vec[i]
    for i, j, k in triples
])
L = wedge.jacobian(vec)
pf = k01*k23-k02*k13+k03*k12
check_zero(L.det()-pf**2)

# Ambient homogeneous quadratic Jacobi equation for P1=w K(s,t).
# Keep all sixty pure-normal quadratic coefficients, rather than setting any
# unneeded coefficient to zero.
s, t, w, u, v = x = S.symbols("s t w u v")
norm_mons = [w*w, w*u, w*v, u*u, u*v, v*v]
P1, Q = S.zeros(5), S.zeros(5)
outer = [0, 1, 3, 4]
for num, (i, j) in enumerate(itertools.combinations(outer, 2)):
    aa, bb = S.symbols(f"ks{num} kt{num}")
    P1[i, j] = w*(aa*s+bb*t)
    P1[j, i] = -P1[i, j]
for num, (i, j) in enumerate(itertools.combinations(range(5), 2)):
    coeff = S.symbols(f"q{num}_0:6")
    Q[i, j] = sum(a*m for a, m in zip(coeff, norm_mons))
    Q[j, i] = -Q[i, j]
P = P1+Q
checked_jacobi = 0
for i, j, k in triples:
    i, j, k = outer[i], outer[j], outer[k]
    jac = sum(P[i, a]*S.diff(P[j, k], x[a])
              +P[j, a]*S.diff(P[k, i], x[a])
              +P[k, a]*S.diff(P[i, j], x[a]) for a in range(5))
    # First normal degree two is exactly the binary-degree-one part.
    restricted = S.Poly(S.expand(jac.subs(w, 0)), s, t)
    linear_binary = sum(coef*s**e[0]*t**e[1]
                        for e, coef in restricted.terms() if sum(e) == 1)
    expected = (Q[i, 2]*S.diff(P1[j, k], w)
                +Q[j, 2]*S.diff(P1[k, i], w)
                +Q[k, 2]*S.diff(P1[i, j], w)).subs(w, 0)
    check_zero(linear_binary-expected)
    checked_jacobi += 1

# B=0 mixed Pfaffian equations, evaluated with arbitrary abstract leading
# coefficients. The commutator relation ell*k=m*q is tracked explicitly.
p, q, r, k, ell, m = S.symbols("p q r k ell m")
cs, ct, Xu, Xv = S.symbols("cs ct Xu Xv")
c3 = -p*k+q*r
bs3 = -(r*w+m*v)*Xv+k*w*Xu
bt3 = (p*w+ell*v)*Xv-q*w*Xu
check_zero(q*bs3+k*bt3 + c3*w*Xv-(ell*k-m*q)*v*Xv)
det = S.Matrix([[k*w, -r*w-m*v], [-q*w, p*w+ell*v]]).det()
check_zero(det + c3*w*w-(ell*k-m*q)*w*v)
numerator_u = (p*w+ell*v)*cs*w**3+(r*w+m*v)*ct*w**3
numerator_v = q*w*cs*w**3+k*w*ct*w**3
check_zero(numerator_u-w**3*((p*cs+r*ct)*w+(ell*cs+m*ct)*v))
check_zero(numerator_v-w**4*(q*cs+k*ct))

# General B=0 first jet, including noncentral {s,t}, and arbitrary Q except
# for the two normal Hamiltonian components already proved above.
def linear(name):
    a,b = S.symbols(name+"0 "+name+"1")
    return a*s+b*t

pp,qq,rr,kk,ll,mm = [linear(name) for name in ("p","q","r","k","ell","m")]
aw,au,av = [linear(name) for name in ("aw","au","av")]
mu,nu,lambda_ = S.symbols("mu nu lambda")
B0 = S.zeros(5)
for (i,j),value in {
    (0,1):aw*w+au*u+av*v,
    (0,3):pp*w+ll*v,
    (0,4):qq*w,
    (1,3):rr*w+mm*v,
    (1,4):kk*w,
}.items():
    B0[i,j] = value
    B0[j,i] = -value
Q[2,3],Q[3,2] = mu*w*w+nu*w*v,-mu*w*w-nu*w*v
Q[2,4],Q[4,2] = lambda_*w*w,-lambda_*w*w

def pfaff4(matrix, ids):
    a,b,c,d = ids
    return matrix[a,b]*matrix[c,d]-matrix[a,c]*matrix[b,d]+matrix[a,d]*matrix[b,c]

def binary_degree(expr, degree):
    return sum(coef*s**e[0]*t**e[1] for e,coef in S.Poly(S.expand(expr),s,t).terms()
               if sum(e)==degree)

allP = B0+Q
bp = [(-1)**i*pfaff4(allP,[j for j in range(5) if j!=i]) for i in range(5)]
X,Y,Z = Q[0,2].subs(w,0),Q[1,2].subs(w,0),Q[3,4].subs(w,0)
check_zero(binary_degree(bp[0],0).subs(w,0)-Y*Z)
check_zero(binary_degree(bp[1],0).subs(w,0)+X*Z)
check_zero(binary_degree(bp[3],1).subs(w,0))
check_zero(binary_degree(bp[4],1).subs(w,0)-v*(ll*Y-mm*X))
lower_zero_sub = {a:0 for a in list(ll.free_symbols|mm.free_symbols) if a not in {s,t}}
for last,expected in [(3,rr*X-pp*Y+av*Z),(4,kk*X-qq*Y-au*Z)]:
    jac = sum(allP[0,a]*S.diff(allP[1,last],x[a])
              +allP[1,a]*S.diff(allP[last,0],x[a])
              +allP[last,a]*S.diff(allP[0,1],x[a]) for a in range(5))
    check_zero(binary_degree(jac,1).subs(w,0).subs(lower_zero_sub)-expected)

# Logarithmic linear derivations of a binary quadratic: both normal forms.
p0,p1,r0,r1,lam = S.symbols("p0 p1 r0 r1 lam")
pp,rr = p0*s+p1*t, r0*s+r1*t
node = S.Poly(S.expand(pp*t+rr*s-lam*s*t),s,t)
double = S.Poly(S.expand(2*s*pp-lam*s*s),s,t)
assert node.coeff_monomial(s*s)==r0
assert node.coeff_monomial(t*t)==p1
assert node.coeff_monomial(s*t)==p0+r1-lam
assert double.coeff_monomial(s*t)==2*p1
assert double.coeff_monomial(s*s)==2*p0-lam

# The squarefree B=0 resonance is an exact consequence of the quartic
# homogeneity and normal divergence equations.
alpha,beta,gamma,delta,rho,xx,yy,sigma = S.symbols("alpha beta gamma delta rho xx yy sigma")
mu0,lambda0,nu0 = -(alpha+beta)/3,-(gamma+delta)/3,-rho*(gamma+delta)/3
sigma0 = S.expand((rho*gamma+rho*delta+nu0)/2)
check_zero(sigma0-rho*(gamma+delta)/3)
F4bar = (s*yy-t*xx)*sigma0*v**4
required = rho*v*(gamma*s*yy*v*v-delta*t*xx*v*v)
diff = S.expand(S.diff(F4bar,v)-required)
check_zero(S.Poly(diff,s,t,v).coeff_monomial(s*v**3)-rho*yy*(gamma+4*delta)/3)
check_zero(S.Poly(diff,s,t,v).coeff_monomial(t*v**3)+rho*xx*(4*gamma+delta)/3)

# Double-root c=s^2: the same identities give two nonzero-delta resonances
# and one separate delta=0 branch. Check the rational identities before
# taking cases, then all substitutions used in the proof.
eta,theta,iota,epsilon,zeta = S.symbols("eta theta iota epsilon zeta")
sigma_double = (eta+3*iota)/6
F4double = (s*yy-t*xx)*sigma_double*v**4
required_double = v*(eta*s*yy*v*v-(theta*s+iota*t)*xx*v*v)
ddiff = S.Poly(S.expand(S.diff(F4double,v)-required_double),s,t,v)
check_zero(ddiff.coeff_monomial(s*v**3)-((4*sigma_double-eta)*yy+theta*xx))
check_zero(ddiff.coeff_monomial(t*v**3)+(4*sigma_double-iota)*xx)
double_sub = {theta:eta*epsilon/delta,iota:eta*zeta/delta}
check_zero(S.expand(((4*sigma_double-eta)*yy+theta*xx).subs(double_sub))
           -eta*((-delta+6*zeta)*yy+3*epsilon*xx)/(3*delta))
check_zero(S.expand(((4*sigma_double-iota)*xx).subs(double_sub))
           -eta*(2*delta+3*zeta)*xx/(3*delta))
# First resonance: X=0, zeta=delta/6, Y nonzero.
check_zero(((-delta+6*zeta)*yy+3*epsilon*xx).subs({zeta:delta/6,xx:0}))
# Second resonance and the allowed upper-triangular change t'=t+a*s.
resonance_sub = {zeta:-2*delta/3,yy:3*epsilon*xx/(5*delta)}
check_zero(((-delta+6*zeta)*yy+3*epsilon*xx).subs(resonance_sub))
check_zero(((2*delta+3*zeta)*xx).subs(resonance_sub))
check_zero((yy-epsilon*xx/(delta-zeta)).subs(resonance_sub))
# In the delta=0 branch eta=0; if iota!=0, X=Y=0 follows from these
# triangular equations. If iota=0 and theta!=0, the residual equation is
# theta*X=0, while Y remains unrestricted by these necessary equations.
check_zero((4*sigma_double-iota).subs(eta,0)-iota)
check_zero(((4*sigma_double-eta)*yy+theta*xx).subs(eta,0)-(2*iota*yy+theta*xx))
check_zero(((4*sigma_double-eta)*yy+theta*xx).subs({eta:0,iota:0})-theta*xx)

result = {
    "all_passed": True,
    "wedge_determinant": "Pf(K)^2",
    "generic_pure_normal_quadratic_coefficients": 60,
    "heisenberg_jacobi_restrictions_checked": checked_jacobi,
    "b0_pfaffian_linear_system_and_cramer_numerators": True,
    "b0_mod_w_pfaffian_and_jacobi_identities": True,
    "binary_quadratic_logarithmic_derivations": ["s*t", "s^2"],
    "b0_squarefree_one_to_minus_four_resonance": True,
    "b0_double_root_two_resonances_and_exceptional_branch": True,
    "scope": "Exact polynomial identities only; geometric proof in TeX."
}
Path(__file__).with_name("verification.json").write_text(
    json.dumps(result,indent=2)+"\n", encoding="utf-8")
print(json.dumps(result,indent=2))

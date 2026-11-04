# Actual proper core: five lines and generic Artin algebras

This new packet uses the independently accepted actual full seven-function core
and its three-generator presentation. It performs no Groebner basis or primary
decomposition search. All statements concern the formal germ at the origin.
The labels below name coordinate lines; they do not identify additional ambient
Poisson components or determine the original D0 ring.

Write x=v2, y=v6, r=x-a-d, s=a+x and P=d²+5. The actual core is

    T=QQ[[a,d,x,y]]/(Q0,Q1,s*r²).

The source-bound Q0,Q1 have 11 and 27 terms. Their relation to the original
seven functions in the proper-source chart is checked again using only the
accepted units D=1-d² and P. The earlier positive source identities identify
the full seven-function ideal with this ideal; the earlier initial regular
sequence proves that T is one-dimensional CM and Gorenstein, of multiplicity12.

## Exhaustive reduced support

Every prime over the core contains r or s. At r=0, set

    h=−a²−3ad+dy,  j=a²d+2ad²−5a+5y.

Then Q0=8h, Q1=−8aj and d*j−5h=a*P*(a+2d). If d belongs to the prime,
h forces a=0, giving N. Otherwise the elimination forces a=0 or a+2d=0,
giving C or R. At s=0, Q0=4d*(−a+2d+2y). If d=0, Q1=10a*(a−2y),
giving N or P. Otherwise the second factor vanishes and Q1 reduces to
2d*(d+2y)*P, giving L. Here P is a formal unit. Literal substitution into
all three source equations verifies the converse for every line:

| Label | Parametrization (a,d,x,y) |
|---|---|
| N | (0,0,0,t) |
| C | (0,t,t,0) |
| R | (−2t,t,−t,−2t) |
| P | (2t,0,−2t,t) |
| L | (t,t,−t,−t/2) |

Thus the radical is exactly the intersection of these five linear primes.
The elementary support argument also works in the polynomial open D*P≠0;
the CM and nilpotence conclusions here use the accepted formal-germ CI gate.

## Actual generic local algebras

All calculations are symbolic over QQ(t), extending to the generic completed
residue fields of the lines. No numerical parameter is chosen. Two-variable
implicit eliminations use explicit2×2 determinants; the simple cases use fixed
3×3 Jacobians. Every claimed second coefficient is checked in the original Q0,Q1.

At N, fix y=t and use a=u. The implicit unit minor for (d,x) is −80t², and

    d=4u²/t+O(u³),  x=−3u+2u²/t+O(u³),
    s*r²=−32u³+O(u⁴).

Hence the local algebra is QQ(t)[[u]]/(u³), with Hilbert function(1,1,1).
In particular d is nonzero and has order2 in this length3 algebra.

At C, fix d=t and set a=u, r=v. The coefficient of y in Q0 is the unit8t,
so y is eliminated by an exact polynomial formula. The resulting Q1 has
order2 and u² coefficient −16(t²+5). Since s=t+2u+v is a unit, the other
relation is v². The two quadratic initial forms are a regular sequence.
The local algebra has Hilbert function(1,2,1), length4, and nilpotence index3.
This statement does not guess a diagonal quadratic normal form.

At R the Jacobian of (Q0,Q1,r) in the three normal variables is a unit, and
s=−3t is a unit; the algebra is QQ(t)[[r]]/(r²). At P, fix y=t and use d=u.
The implicit minor for (a,x) is −240t². The first derivatives are a'=1,x'=−1,
and the checked quadratic coefficient of s is −1/(3t). Since r=−4t is a unit,
the algebra is QQ(t)[[d]]/(d²). Finally the Jacobian of (Q0,Q1,s) is a unit
at L and r=−3t, so the algebra there is reduced of length1.

The generic lengths in order N,C,R,P,L are therefore3,4,2,2,1. Their sum12
agrees with the accepted actual multiplicity. Equivalently, the r divisor has
lengths1,2,1 at N,C,R; the s divisor has lengths1,2,1 at N,P,L. The CM core
has no embedded primes. It injects into its five generic localizations, so the
computed generic nilpotence bounds prove Nil(T)³=0 globally. Its square is
nonzero at N and C, proving that the exact index is3.

## Whole normalized E4

Retain the five multiplied pivot relations. The whole ring is

    E4=S/(dJ),  S=QQ[[u1,...,u5,a,d,x,y]],
    J=(u1,...,u5,Q0,Q1,s*r²).

The exact sequence0→S/J --times d→ E4→S/(d)→0 shows that its associated
primes are exactly (d) and the five core line primes (with the u variables).
The line primes N and P contain d and are embedded. The other three line
primes and (d) are minimal.

At N and P, the localized radical of dJ is (d). Under multiplication by d,
the nil ideal is the corresponding core module; its square corresponds to
d*T and its cube to d²*T. The computed generic models give d²=0 and d≠0
at both lines. At C,R,L the parameter d is a unit, so the nilpotence bound
is the core bound. At the generic prime(d) the ring is reduced. A submodule
vanishing at every associated prime is zero. Consequently Nil(E4)³=0 and
Nil(E4)²≠0: its exact nilpotence index is3. This is normalized E4 only.
Its previously proved dimension8/depth1 do not determine the original D0
nilradical, original intersection depth, or any full Poisson scheme.

The finite certificate stores the source restrictions, exact elimination
identities, all five literal substitutions, generic unit minors, second
coefficients and local Hilbert functions. The single bounded run uses only
small fixed polynomial arithmetic under60s/2GB/32MiB/150000 terms.

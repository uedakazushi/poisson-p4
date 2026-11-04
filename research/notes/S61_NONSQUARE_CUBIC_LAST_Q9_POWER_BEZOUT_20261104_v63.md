# A fixed a^29 Bezout certificate for the last quadratic row

The nine queried coefficients of the actual last quadratic source row
admit the explicit finite-DAG identity

\[
 \sum_{i=0}^8 H_i h_i=a^{29}.
\]

Every `H_i` is a polynomial over QQ in `a,e,d`, represented by the saved
circuit. Its verified degree upper bound is at most 207. Thus
`(H_i/a^29)` is an actual unit functional on these nine entries of the
full remaining row on FULLU; the other thirteen entries have functional
coefficient zero. No free target complement is inferred from pointwise
rank. A target completion has not yet been executed.

The code and report stem is
`S61_NONSQUARE_CUBIC_LAST_Q9_POWER_BEZOUT_20261104_v63`.
The successful fixed construction took 0.749 seconds and saved 25,572
scalar nodes, 145 polynomial atoms and a 591 KB certificate. It uses
no GB, point, coefficient-degree guess, minor search or unknown inverse.
Independent review is pending. The physical chart application is on
FULLU only; polynomial continuation of the coefficient formulas is not
an assertion about physical points with C=0.

## 1. The source binding, including both target changes

Use the old eight-step checkpoint rows 4, 5 and 9. The last source
combination before the target changes is

`F = row9 - c*row5 + (c*d5-d9)*row4`,

where `c` is the explicitly saved nine-term polynomial, and `d5,d9`
are the earlier GL6 source-clearing coefficients. Row 4 is a combination
of `q0*x0` and `1600*q5*x0`; its entire polynomial vanishes at x0=0.
The code tracks the ten quadratic source rows through every fixed saved
operation and verifies this fact. None of those operations introduces
a `q5*L` term into another source row. The coefficient of `1600*q5*L`
in the surviving relation remains exactly one.

The first GL6/tail change acts only within the x0-divisible target
subspace. The second GL3/tail change need not preserve the ideal
`(x0)`, and that preservation is not used. Instead the exact source
combination above is reconstructed before both target changes.
Let the two original L coefficient triples be `b,s`, and put
`v=(6(18e+7)/(5a),-18,-9)`, so `b*v=1`, `c=s*v`. Then its three
remaining L coefficients are `s-c*b`. Their dot with `v` is zero.
The final target labels 22 and 25 give its first two original L
coefficients; the third is recovered using the unit `v[2]=-9`.

Write the queried coefficients in their frozen order
`h=(h11,h13,h14,h15,h17,h18,h19,h22,h25)`. At x0=0 the source has the
form `A(y)+L*B(y)` modulo L², where

\[
\begin{aligned}
A={}&h_{11}y_1^2y_2+h_{13}y_1y_2^2+h_{14}y_1y_2y_3
 +h_{15}y_1y_3^2+h_{17}y_2^2y_3+h_{18}y_2y_3^2+h_{19}y_3^3,\\
B={}&h_{22}y_1y_2+h_{25}y_1y_3
 +(v_0h_{22}/9-2h_{25})y_3^2.
\end{aligned}
\]

Let `P=q0|x0=L=0`, `R=1600*q5|x0=L=0`, and let `ell(y)` be the
coefficient of L in `1600*q5|x0=0 mod L²`. These are read directly
from the initial source-bound quadratic rows. The saved source
cofactors give a linear form `m`, another linear form `n` and scalar
`lambda`. All coefficients of the following identities are checked:

\[
 A=Rm+Pn,\qquad B=R+\ell m+\lambda P.
\]

This binds the monic-L source argument to the actual nine coefficients,
including all earlier source and tail operations.

## 2. Geometric unit-ideal proof without a chosen dual

The pair P,R has the accepted literal coefficient dual. The symmetric
matrix of P on `(y1,y3)` has determinant `-25/16`; therefore P has
rank at least two over every characteristic-zero residue field and
is geometrically squarefree. If all nine h vanish in a residue field,
the source identities imply

`R*m=0 mod P`, `R=-ell*m mod P`, and hence `R²=0 mod P`.

Squarefreeness gives `P|R`. Since P and R have the same degree, this
makes them proportional, contradicting their coefficient dual.
Thus the nine-entry ideal is the unit ideal on FULLU. This argument
alone would not compute a particular dual, its denominator, or a
chosen free complement. The following construction supplies them
explicitly at the level of the row functional.

## 3. A constant-determinant conic frame and its functional

Make the polynomial coordinate change

\[
 S=-5y_1/2+(9e+4)y_2/5,\quad T=y_2,
 \quad W=y_3-S/5-(3e+8)T/25.
\]

Its only inverse denominators are nonzero rational constants. Direct
coefficient comparison gives

\[
 P=SW+qT^2,\qquad q=(625a+108e^2+156e+48)/125.
\]

Use the monic reduction `SW=-qT²`. Write the normal form of R as

`Rbar=r0*S²+r1*S*T+r2*T²+r3*T*W+r4*W²`.

Transport the full accepted pair functional phi, with `phi(P)=0`
and `phi(R)=1`, to these five monomials. Its coefficients `t_i`
have denominators dividing a². The code verifies that
`T_i=a²*t_i` are polynomials and checks the literal identity

`sum(T_i*r_i)=a²`.

No unknown q, conic discriminant or parameter factor is inverted.

## 4. Five explicit universal power witnesses

The nine normal coefficients of `Rbar² mod P`, in the ordered basis
`S4,S3T,S2T2,ST3,T4,T3W,T2W2,TW3,W4`, are

\[
\begin{split}
J=(&r_0^2,\ 2r_0r_1,\ r_1^2+2r_0r_2,
2r_1r_2-2qr_0r_3,\\
&r_2^2-2qr_1r_3+2q^2r_0r_4,
2r_2r_3-2qr_1r_4,
r_3^2+2r_2r_4,\ 2r_3r_4,\ r_4^2).
\end{split}
\]

The certificate contains explicit nine-component polynomial witnesses
for membership of

`r0², r1⁴, r2⁶, r3⁴, r4²`

in `(J)`. All five identities are checked in the universal integer
polynomial ring `ZZ[q,r0,r1,r2,r3,r4]`. They are obtained by the fixed
monic substitutions for `r0²,r4²,r1²,r3²,r2²`, not by computing a GB.
The recorded rewrite counts are `1,3,20,3,1`.

The bounds `(2,4,6,4,2)` have `sum(bound_i-1)=13`. Consequently every
degree-14 monomial in r contains one of those five powers. For each
of the 3,060 exponent tuples summing to 14, the code chooses the
first applicable fixed witness and attaches its exact multinomial
coefficient. This gives an explicit coefficient circuit `C_j` with

`sum(C_j*J_j)=(sum(T_i*r_i))^14=a^28`.

The finite proof is verified term by term through the five saved
universal identities and all multinomial coefficients. Large core
parameter polynomials are not expanded to verify this equality.

## 5. Return to the actual nine coefficients

The checked source relation gives the identity

`R*B-ell*A = R² + P*(lambda*R-ell*n)`.

Normal reduction by P therefore expresses every J as a linear
combination of the nine actual h. The recovered third L coefficient
has only an a denominator. The code constructs and verifies the
polynomial matrix identity

`a*J = M9x9*h`

in every coefficient, with no denominators other than rational
constants. Substituting it into the fixed power certificate yields

`sum(H_i*h_i)=a*sum(C_j*J_j)=a^29`.

The complete 25,572-node circuit, all 145 polynomial atoms, all five
universal witnesses, every actual source-membership coefficient and
all 3,060 term choices are saved. The reported degree bound 207 is
obtained from the actual atom degrees by addition/max propagation
through the circuit. The cofactors themselves remain unexpanded.

## 6. Next completion is explicit but not executed

The functional `(H_i/a^29)` extends by zero to the other thirteen
coordinates of the full actual row. If the remaining row is written
as `(h9,b13)` and the unit functional as `c9`, an elementary target
change can clear all thirteen tails using `h9*c9=1`, leaving
`(h9,0^13)`. One zero target then supplies the explicit involution

`[[I9-c9*h9,c9],[h9,0]]`,

whose square is the identity and determinant is -1. It sends
`(h9,0)` to a unit vector. This avoids an unspecified projective
complement. That transport is a separate next stage and was not
performed here. Full cubic57, full quartic239 and kernel38/gap42
are unchanged by the present dual certificate.

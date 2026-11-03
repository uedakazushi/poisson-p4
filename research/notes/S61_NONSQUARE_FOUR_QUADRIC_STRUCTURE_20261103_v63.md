# Small structural model for the raw four quadrics

The source is the four cleared quadrics of
`S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62`, after d=P/a and row
multipliers a^[0,2,2,1]. Denote them Q0,Q1,Q2,Q3 in light variables x0,...,x5;
x4 is the earlier Z coordinate and x5 the earlier W. Work on the original
known-unit core a*P*F*B0*C*J!=0, J=5P−6e−4. All identities below are exact QQ
polynomial identities. There is no new field or large GB calculation.

## Verified normalization and generic geometry

Replace Q2 by Q2'=(Q2−eQ1)/(a²J). This is a polynomial quadric of core degree3,
with light support in (x0,x5), replacing a core-degree7 row by a core-degree3 row.
Only the known unit a²J is divided. The ideal and linear-syzygy dimension are
preserved at each point of this open core.

Set U=x4−5x0+4x2, V=x1, T=x1+2x3 and
R=a(U+x2)+(21e+16)x1/10+(9e+4)x3/5. Then

4Q0=U R−5 V T.

Together with r=x0 and z=x5 these are invertible coordinates over a!=0; the
forward determinant is ±2a and both inverse compositions are checked exactly.
Thus Q0 has rank4 and a two-dimensional radical uniformly, with no unknown
factor inverted. This gives two P1 families of maximal Q0-isotropic P3s, useful
for a later common-linear-space test; no all-parameter absence is claimed here.

Generic complete-intersection behavior is possible and is certified. At
(a,e,P)=(1,1,1), where every known unit is nonzero, restrict by
x4=y0+y1+2y2+3y3, x5=2y0+3y1+5y2+7y3, xi=yi for i<4.
The complete degree5 Macaulay map of the four quadrics has exact QQ rank56
(shape56-by80). A selected56-by56 nonzero rational determinant and the entire
matrix are saved. Hence the four quadrics form a regular sequence over
QQ(a,e,P). This does NOT prove regularity at all core points.

## Exact factorization and a conditional syzygy reduction

Write Rc=15P²−6Pe+6P+80a, Tc=−5P+6e−2,
Sc=21Pe+16P+140a, Vc=9Pe+4P+60a,
Kc=85P²−12Pe+58P+600a. Define

B=30a Rc*x0−9P Sc*x1+90Pa Tc*x2−18P Vc*x3−2a Kc*x4,
A=(B+10aPJ*x4)/3.

On x5=0 the verified identities are

Q1=A B/6075,  Q2'=40x0 B/6075,  B−3A=−10aPJ*x4.

Both B and A are independent of x5; A is independent of x0 modulo x5 on P!=0,
since 3Sc−7Vc=20P. Define the actual polynomial quadric

D=(40x0 Q1−A Q2')/x5.

At a fixed core point, suppose:

1. Q0 and Q3 are coprime after imposing B=x5=0.
2. Their restrictions to x0=A=x5=0 are linearly independent.
3. The five quadrics Q0,Q1,Q2',Q3,D are linearly independent.

Then the original four quadrics have no linear syzygy. Indeed reduce a relation
sum Li Qi=0 modulo (B,x5). Coprimality and degree force
L0=c0 B+d0 x5, L3=c3 B+d3 x5. Reduction modulo x5 and division by B, then
restriction to (x0,A), force c0=c3=0 by gate2. Coprimality of the independent
linear forms x0,A gives L1=40k x0+d1 x5 and L2=−k A+d2 x5. Division of the
whole relation by x5 leaves the constant relation

d0 Q0+d1 Q1+d2 Q2'+d3 Q3+kD=0.

Gate3 forces every coefficient to vanish. Conversely any nonzero constant
relation of the five quadrics constructs a nonzero linear syzygy by these same
formulas. Thus, under gates1–2, gate3 is equivalent to the required raw rank24.
No gate has been assumed merely from generic behavior.

The five-quadric coefficient matrix is only5-by21, with846 core coefficient
terms. Exact known-unit/constant-span operations using ONLY a,P,J give two
safe pivots and a3-by19 residual. Full pointwise rank5 is not certified in this
packet. The source, native residual, exact helper and logs are saved for the
next small-minor computation.

## Explicit limits of easy cuts and of the restricted-rank shortcut

The fixed easy cut family x4=lambda*x0, x5=mu*x0 has an x0=0 branch. A known-P
coordinate chart reduces it to two binary quadrics. Their fixed-degree
homogeneous resultant is

9*P^4*a^5*J²*H/1600000,

with a nonunit residual polynomial H of degree12, saved exactly. H cannot be
inverted as a known unit. Its vanishing concerns this chosen cut and does not
prove a common component of the full four-quadric ideal.

On B=x5=0, generic Q0,Q3 are coprime, but the uniform gate remains open. In the
split Q0 coordinates the B coefficient of the radical coordinate r is
−20a*DB, DB=20P²+3Pe+20P+180a. If restricted Q0 has rank at most2, then DB=0
and the dual norm of its remaining B covector is zero. Eliminating a by DB
reduces the latter equation to an explicit quartic in(P,e), up to known factors.
That exceptional curve is retained as a case equation, never inverted. Neither
the generic gcd nor the generic complete-intersection receipt removes it.

The exact row reduction, normal form, factorizations and generic CI receipt are
proved here. The uniform three gates, all-core raw rank24/original product41,
quartic coverage and ambient removal remain unproved. Independent review of
this structural packet is pending.

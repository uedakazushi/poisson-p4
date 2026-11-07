# Every finite-wall lift belongs to the main T37 Lagrangian-curve closure

2026-11-07. Pure hand geometry using the frozen exact source equations and the finite-field proof. No new local, kernel, Hessian, node or anchor computation is needed. This is a new completion of the historical finite-wall frontier; the old unresolved checkpoints and their statements remain immutable.

## 1. The actual lifted source and its finite projection

Let `u=a4,v=a5`, and let B be the reduced rationally irreducible factor

\[
B=373248u^5+82944u^3v-139968u^2v^4+3888uv^2-31104v^5+1.
\]

Define X_B by the **full original six Lagrangian equations** in the exact cell012 plane chart, together with B=0. The saved seven-polynomial Gröbner basis is an equivalent presentation of those equations. Let `π:X_B→B` forget a0,a1,a2,a3.

The first three Gröbner rows solve a0,a1,a2 polynomially in a3,u,v (their own coefficients are one). The fourth row, index3, is quadratic in a3 with leading coefficient v. Thus over `v≠0`, a3 is integral of degree at most two over the base. Over the open where the coefficient c4 of a3 in row index4 is nonzero, that row solves a3 linearly and the projection is finite as well.

These two base opens cover B. Indeed, up to a nonzero rational constant,

\[
c_4=72u-136048896v^8-6480v^3.
\]

At v=0 this is 72u, whereas `B(u,0)=373248u5+1`; B has no point with v=c4=0. Finiteness is local on the base, so **π is finite on the whole B curve**, not only a generic chart.

The exact generic rational source verified in `T37_DOUBLE_CONIC` supplies a point over Q(B). Let Y be its reduced closure in X_B. Then Y is finite over B and dominant. A finite morphism is closed, hence its image is the whole B. In particular **every arithmetic closed point of B has a nonempty Y fiber**. No assertion that the projected curve is smooth or nodal enters this argument.

## 2. Where can a fiber contain more than the generic lift?

Outside the simultaneous zeros of c4 and the a3 coefficient c5 of row index5, one of those linear rows determines a3 uniquely. The first three rows then determine the remaining variables. Every such fiber has at most one support point. Since Y surjects onto B, that support point belongs to Y.

The saved exact univariate patch calculation gives the entire simultaneous-zero locus:

\[
u=\phi(v)=v^3(1889568v^5+90),
\qquad r(v)=408146688v^{10}+42768v^5+1=0.
\]

The possibility v=0 adds no point of B, as above. Substituting this relation into all source equations yields exactly the saved algebra

\[
R_{20}=\mathbb Q[v,a_3]/\left(
 v^{10}+\frac{11v^5}{104976}+\frac1{408146688},\;
 a_3^2-(349920v^7+15v^2)a_3+288684v^9+\frac{99v^4}{2}
\right).
\]

The first three rows give the recorded actual a0,a1,a2 map. Every original Lagrangian equation and every Gröbner row reduces to zero in this algebra; the extra-relation ideal has rank zero. Thus this is the **whole actual fiber source**, not just an ambient algebra containing it and not a selection of twenty anchors.

## 3. The arithmetic fiber has no proper nonempty closed subset

`T37_FINITE_QUARTIC_FIELD_DESCENT.md` proves R20 is a field of degree20, without assuming reducedness. Its base field at the wall is

\[
F_{10}=\mathbb Q[v]/(r),\qquad [F_{10}:\mathbb Q]=10.
\]

For completeness, put `w=v5`, `U=34992w`. Then `Q(U)=Q(sqrt(13))` and

\[
N_{\mathbb Q(U)/\mathbb Q}(w)=2^{-8}3^{-13}.
\]

Some place above 2 has valuation of w not divisible by five. The valuation equation for a root of v5=w forces extension ramification divisible by five, so `[Q(v):Q(U)]=5`. This proves the degree10 statement. The established quartic-field argument then gives `[R20:F10]=2` and irreducibility of its quadratic a3 equation.

The wall is therefore **one arithmetic closed base point** `Spec F10`, and its entire actual fiber is `Spec R20`. The Y fiber over this point is nonempty by finite surjectivity. It is a closed subset of `Spec R20`. A field spectrum has no proper nonempty closed subset, so this fiber is the whole fiber.

After extension to C, all twenty geometric embeddings are consequently contained in Y. This uses arithmetic field irreducibility to control all conjugate points at once. It does not select a complex root and does not infer PGL equivalence from Galois conjugacy. Equivalently, the closed-fiber ideal in the field R20 is zero; this conclusion persists after base change.

Combining this with the singleton fibers outside the wall proves

\[
(X_B)_{\mathrm{red}}=Y
\]

on support. In particular there are no extra isolated support components above the simultaneous-zero locus. We do not need to identify the analytic branches or the scheme's embedded/nilpotent structure.

## 4. The exact generic LQ² identity extends to every actual B-plane

The frozen generic-field computation proves

\[
D|_{\mathrm{plane}}=LQ^2
\]

by all21 ternary coefficient identities over Q(B). This is an actual divisor identity, not a cofactor comparison. The projective multiplication image of `(L,Q)↦LQ²` is closed because its source is proper. The actual restricted quintic is never zero in this cell: its C5 coefficient is the constant1/3. Hence the divisor map from X_B to projective quintic space is regular everywhere, including the wall.

The generic image of Y satisfies this closed condition. Therefore all of Y satisfies it; by §3 this means **every actual lifted B-plane**, including the whole finite wall source and every rational-function denominator complement.

This completes the particular uniform normal-cone entrance obstruction for the reduced H plus irreducible quartic targets E6 smooth-flex, A5 smooth-hyperflex and Ga9 linear-w. Those reduced1+4 divisors cannot equal LQ². It also excludes the two stable divisors `H²LQ_irreducible` and `H²L1L2L3`: neither admits a squared degree-two factor. All actual plane parameters, the full original six equations and every source-chart complement are retained.

The conclusion is about this whole plane-source branch and its uniform normal-cone construction. It does **not** exclude the entire T37 component as a parent via another degeneration, other plane branches, moving genuine periods, or nonuniform coordinate contractions. No new component dimension, membership or count is asserted.

## 5. Consequences for the stopped calculations

The strict-linear-symmetry QQ kernel calculations were unnecessary for this support statement. Their two real120-second failures remain accurately recorded, with no kernel conclusion inferred. The new finite/field/closed-image argument uses their completed finite source presentation, not the uncompleted matrix calculation.

A native factor multiplication replay remains useful as independent arithmetic evidence. The root's exact zero plane-Jacobian determinant is consistent with these wall points belonging to the Lagrangian curve closure, but the closure proof does not depend on that determinant or on the uncompleted native factorization. No new matrix retry is requested.

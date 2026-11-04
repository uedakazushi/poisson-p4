# Independent additional norm, formal pencil and necessary residual-cover review

This review covers the frozen additional-minor norm and norm-pencil
packets, including their projection-curve/residual corollary. It reuses
the accepted B35 source-minor and earlier pair/monic/norm receipts.
No author worker, GB, nested resultant, point or prime is run.

The new source k is exactly normalized generator20, corresponding to
original minor33 in the lexicographic 35-minor order. Its 88 source
coefficients are independently matched to that primitive minor and
the accepted source-unit normalization. The monic degree-four H
is the exact preceding source-bound H, with only a and rational
constants inverted. The review independently performs the four
monic divisions for both g and k and compares all 32 multiplication
matrix entries and all eight quotient witnesses.

For the new k norm, a separate native arithmetic checker uses a
different monomial order and last-row Laplace expansion. It verifies
the fixed cleared four-by-four determinant, all sixteen adjugate
cofactors, all sixteen matrix-adjugate products, the saved first
adjugate column, and the complete Laurent identity
`Normk = cofactorH*H + cofactork*k`.
This checks all coefficients, including d cancellations. It does
not infer a determinant from a single adjugate product alone.

The native arithmetic source is newly written for this review;
no author arithmetic helper is imported. Polynomial operations use
exact QQ coefficients, with a Laurent offset only in a. Its fixed
arithmetic task file is saved for replay. All removed factors are
checked by literal multiplication, including their quotients in
`U=adFBCJ`; hence every removed factor is a unit on FULLU. Final
coprimality with U and all common-curve gcds are independently
checked. Unknown curve or determinant factors are never inverted.

For the formal pencil, the checker forms its five coefficients
by sixteen row-multilinearity determinants, rather than the
author's permutation/linear-lambda expansion. It verifies all
five complete QQ coefficient arrays, not samples in lambda.
The endpoints are exactly `(5000*a)^4*oldNormg` and the new Normk.
Every content factor, all 130 saved unit-removal steps across the
two packets, every common-curve quotient and the three final gcd
steps are checked. The final common G is unchanged and has 49
terms, e-degree12 and a-degree6.

Over an algebraically closed characteristic-zero field at a!=0,
the determinant of multiplication by `g+lambda*k` in the monic-H
algebra is the product of the linear factors
`g(r)+lambda*k(r)` over H roots with multiplicities. All five
coefficients vanish exactly when at least one H root annihilates
both g and k. Repeated roots cause no problem; etaleness is not
assumed. Using just the two endpoint norms would allow different
H roots, which is why the three mixed coefficients are retained.
Neither that common-root condition nor the surviving G curve is
sufficient for all the other B minors to vanish.

The necessary set remains the union of two pieces. If h is the
old stripped norm and K the new stripped norm, the exact gcd gives

`V(h,K)=V(G) union V(h/G,K/G)`.

Both pieces are retained. The residual quotients have total degrees
16 and25 and exact gcd1. Plane Bezout therefore bounds the number
of geometric residual projection points by400. The monic degree-four
H has at most four distinct d roots over each projection point,
giving at most1600 geometric residual points before FULLU filtering.
This is an upper bound for a necessary exceptional set, not a proof
that any such bad point exists.

The term "curve48" has the precise module meaning

`QQ[a,a^-1,e,d]/(G_monic_e12,H_monic_d4)`

is finite free of rank48 over `QQ[a,a^-1]`, with basis `e^i*d^j`
for `i<12,j<4`. It is not an assertion that the ordinary projective
degree of that affine curve is48. The review checks the constant
e-leading coefficient and the exact saved monic cover. Adding
the remaining FULLU inverses is not automatically finite over the
base; no flatness/Nakayama argument is inferred from these norms.

The accepted outcome is only the source-bound necessary curve-cover
plus finite-residual containment. There is no B-bad existence,
emptiness, actual B-curve component, global H7, full cubic57,
quartic239, or source-kernel conclusion.

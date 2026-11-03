# Uniform product rank at least35 on the nonsquare J≠0 chart

2026-11-03. This proves an all-point lower bound, not product rank41,
full cubic rank57, a free ambient cokernel, quartic coverage or ambient removal.
Work over characteristic0 on

    U = {adFBCJ ≠ 0},   J=5ad−6e−4.

In serialized matrices the third core variable `d` means **P=ad**. Thus
`B0=(50a+18e²−9e−2)P+(60e−40)a=aB` and `J=5P−6e−4`.
The source is the already checked exact six-light quotient by ell; its
four retained quadrics are q0,q1,q3,q5, with six coordinates
`(v0,v2,v4,v7,Z,W)`. Each saved cleared quadric differs from the actual
one by a unit on U. Substituting original d=P/a and multiplying the four
rows by a powers `(0,2,2,1)` gives a polynomial coefficient matrix M of
size4×21 over QQ[a,e,P].

## Polynomial certificate for independence

The saved21×4 polynomial matrix L satisfies, in all16 entries,

    M L = −a²P²(5P−6e−4) identity4.

The verifier reconstructs M from the public coordinate report, checks the
original-d/P conversion term by term, and multiplies M L by independent
Python Fraction sparse arithmetic. L has2461 terms. The multiplier is a
unit on U, so the four restricted quadrics are independent at **every**
geometric point of U. The discovery used M2's module Gröbner basis and
annihilator; the final certificate does not rely on trusting its output
as a Gröbner basis. The certificate itself is an ordinary QQ polynomial
identity.

## Elementary shadow bound

For four independent quadrics in six variables over any field, row reduce
their coefficients to obtain four distinct leading degree2 monomials.
Each leading monomial has six degree3 multiples. For two distinct degree2
monomials, their sets of degree3 multiples intersect in at most one monomial:
a common cubic must equal their least common multiple. Inclusion-exclusion
through pairs therefore gives at least `4·6−binom(4,2)=18` distinct leading
monomials among the24 linear products. Choosing one product for each of
these monomials proves linear-product rank at least18.

This applies pointwise to the four quadrics modulo ell.

## Restoring the seven light variables

The accepted pair cover gives ell,Z,W independent at every point of U,
and q0,q1,q3,q5 independent modulo (Z,W). The original six-quadric ideal
is exactly `(q0,q1,q3,q5,ell Z,ell W)` via a unit row transformation.

Restriction of its degree3 span modulo ell has rank at least18. Its kernel
contains `ell·(Z,W)·S1`, of dimension13 in seven variables. It also contains
the four cubics `ell qi`. Their classes modulo the preceding13-dimensional
space are independent: cancelling ell reduces any relation to a constant
combination of the qi lying in `(Z,W)`, excluded by the pair cover. Thus
the restriction kernel has dimension at least17. Consequently

    rank(original six quadrics × seven light variables) ≥ 18+17 = 35

at every geometric point of U. This concerns the original42 product rows,
not merely the free abstract row module. The known universal syzygy still
gives the separate upper bound41.

## Consequence for finding the next constructive pivot

The primary full-cubic reduction has31 product-row unit pivots. Its ten
remaining product rows therefore have pointwise rank at least4 on U.
Let N be this10-row residual product matrix. For a fixed point p, the
projective set of coefficient vectors `[lambda]` with `lambda N(p)=0`
has dimension at most5. Since dim U=3, the full incidence in `U×P9` has
dimension at most8. Its image closure is a proper closed subset of P9.
A general constant QQ row combination therefore has entries generating
the unit ideal of QQ[U]. A coefficient of lambda is a nonzero QQ unit,
so this combination can be introduced by an invertible row replacement.
This is an existence theorem for a further scalar unit row, not an
explicit Bézout lift or a theorem that all remaining rows are unimodular.
No claim about the J=0 chart follows from this argument.

## Separate constructive scalar-row progress

On the sufficient five-light20×35 criterion, the saved11×26 residual's
row6 has a QQ polynomial Bézout lift of `−a²P²J`:163 terms, cofactor degree3.
Its full original-row identity was checked by M2 and the native QQ helper.
No cofactor is a declared unit, so the image column is appended; one
unit pivot yields10×26, not10×25. The result has124869 terms. This is the
tenth known-unit pivot of that sufficient criterion. It does not complete
its rank20 proof. The previous degree1 searches remain valid refusals
only in their bounded coefficient class.

Replay of the main polynomial certificate:

    .venv/bin/python code/s61_binary_nonsquare_quad_shadow_20261103_v62.py --output work/quad-shadow-replay.json

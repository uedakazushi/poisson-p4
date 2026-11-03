# Nonsquare: bounded exact unit lift on the thin product presentation

This is an exact algebraic reduction on `adFBCJ != 0`, not an all-point rank41
or cubic rank57 theorem. Here `J=5ad−6e−4`. The mod-ell chart and the pinned
24×56 product matrix are those of the accepted coordinates/product packets.
Their ten known-unit pivots give a 14×46 matrix M. Its zero-based row1 has
nineteen nonzero entries p_j.

The new exact identity is

    sum_j c_j p_j = a J,    deg(c_j) <= 3.

There are154 nonzero coefficient terms. The older ordinary-GB lift had199
terms and coefficient degree7. The new certificate is independently checked
by native QQ polynomial arithmetic and Python Fraction multiplication. At
source column22 its witness coefficient is a declared unit. Replacing that
column by M*c is therefore invertible after localization. One exact unit
pivot leaves a 13×45 presentation with225656 terms, giving eleven constructive
unit pivots in the mod-ell 24×56 map. The saved residual is reconstructible;
its size does not establish full rank13 at every core.

The reusable producer constructs the finite Macaulay coefficient matrix for
bounded multipliers of the source polynomial vectors, with the negative
target as the final row and its coefficient fixed to1. Sparse elimination
modulo p=2^61−1 chooses independent equations and pivot unknowns. If small
rational reconstruction fails, an exact QQ Dixon solve of that square system
produces a candidate. Success requires verification of every original QQ
polynomial identity. Thus the modular calculation only selects a system;
no semicontinuity, parameter specialization, or flatness inference is used
for the positive certificate. For this lift the matrix is423×381 and the
selected square system is258×258.

A refusal is not QQ nonmembership. In particular an unlucky-prime fixture
refuses the true identity p*a*(1/p)=a, while a large-coefficient fixture forces
and passes the exact QQ Dixon fallback. The finite bounds on the recorded
positive certificate are literal bounds, not a claim about all syzygies.

The stronger isolated-column target M*c=(aJ)*e_1 was also checked at multiplier
degree<=5. Its finite coefficient matrix has2577 columns, and reduction modulo
the same prime has full column rank2577. Since all input coefficients are
p-integral, QQ rank is at least2577 and at most2577. Consequently that particular
isolated unit-column identity has no solution within degree<=5. This does not
exclude a larger degree or a different unit target. The scalar row identity
above deliberately permits nonzero images in the other thirteen rows.

No free-cokernel, all-core product-rank, all-point cubic lower bound, quartic
length65 coverage, or ambient exhaustion follows from this packet. The generic
QQ(e) module attempts remain separate capped work and are not certificates.

# Nonsquare: two short exact unit pivots (V9)

On the existing open adFBCJ≠0, J=5ad−6e−4, two new constant-coefficient identities give
actual rank-preserving reductions. They do not prove all-point rank41/57, a free cokernel,
the quartic bound65 or ambient exhaustion. Existing V8 is unchanged.

For the pinned mod-ell product14×46 matrix M, zero-based row10 satisfies

    (−14985/8) M_(10,11) + (2997/8) M_(10,13) = a³ J².

Replacing column11 by that same combination in EVERY row is invertible: its coefficient
at the replaced column is the nonzero QQ constant−14985/8. The target entry is a unit
on the declared open. One exact native unit pivot leaves13×45 with44143 terms, versus
the earlier154-term scalar witness whose pivot left225656 terms. This gives eleven
constructive unit pivots in the mod-ell24×56 presentation. Full rank24 is still open.

For the MAIN full-cubic27×61 matrix N in V6's P=ad chart, zero-based row4 satisfies

    −(N_(4,5)+N_(4,30))/50 = C J,
    C=100a+18e²−9e−2,  J=5P−6e−4.

The six declared units there are a,F,B0,C,J,P, with B0=aB. Column5 replacement is
invertible over QQ, then the known-unit pivot gives26×60 with1290701 terms. Combined
with V6's thirty earlier pivots, full cubic rank is31 plus the rank of this26×60
matrix. Its all-point lower rank remains unproved. The global free ROW module57
result does not make this cokernel free or preserve pointwise rank automatically.

Both displayed polynomial identities were checked independently by Python Fraction
arithmetic on the original saved entries. The main source SHA is the pinned V6
P-clean27 SHA; the product source SHA is in the public PRODUCT_ALT evidence. Native
replay uses the same combinations in all rows, checks known-unit factors before any
row division, and saves the full residual SHA. The constant replacements have
nonzero rational determinants; they are not claimed to have determinant1.

## Discovery and reproducibility

The new search removes ONLY the declared known-unit factors from each scalar entry,
solves a QQ constant span in those primitive entries, then restores the common
product of the active denominators to obtain an actual polynomial column witness.
It verifies the entire original target identity before success. An optimized version
selects a square coefficient system from the source span and tests known-unit targets
without row-reducing all target columns at once. A search refusal is not a statement
about the full localized row ideal or arbitrary-degree multipliers. In the final two
identities the restored witness coefficients are simply the two constants above.

`research/handoff/20261103_v62/v9_short_unit_pivots/` contains the frozen sources,
witnesses, original target rows, small product input, native logs, QQ identity/factor
receipts and omitted-output hashes. The89MB main input and94MB residual are regenerated
by V6 and the one-pivot replay. The active known-unit lists contain no extra factor
from the separate regular-sequence experiment.

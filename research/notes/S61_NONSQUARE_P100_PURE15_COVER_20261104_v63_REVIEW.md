# Whole P100 rank149: independent fixed-minor cover review

The whole-chart lower bound is accepted from the actual source-bound minor
cover. The ring is exactly

    QQ[z,1/(z(z+20)(z+60)(203z+2700)(z+12))].

The independent review reuses the accepted original q5H/G18 eight-source frame
and the accepted full corrected K2 seven-source frame, pinning both peer
reports and all their source hashes. K2 is the actual corrected C23 source,
not an abstract pure vector. The old source row order and the eight target
columns are unchanged. The target98 columns are(0,1,2,3,4,5,6,9), with original
122 labels(0,1,2,4,5,6,7,10).

A separate FLINT checker verifies both old8×8 adjugate products in every entry,
128 polynomial identities. For all seven new source rows it builds the LCD
from the five known chart-unit factors and the frozen exponent lists. Each
of the56 selected rational coefficients is cleared by an exact polynomial
division and remultiplication. Its denominator is therefore a divisor of a
known whole-chart unit; no unknown polynomial inverse is used.

Only the fixed fifteen replacements are computed: new physical light v0 at
old selected row positions0,...,7, then v1 at positions0,...,6. For each one
the exact polynomial is the dot of the cleared source row with the indicated
column of the old adjugate. Forward and reverse summation agree. The rows
retain their old order, so no additional permutation sign is introduced.
No other replacement, source coefficient or rank is searched.

Every computed full minor is compared to the frozen polynomial. Its saved
QQ scalar, five known-unit multiplicities and monic primitive factor are
reconstructed literally. The primitive factors are checked not to have any
remaining displayed linear unit factor. This is named-factor arithmetic,
not a GCD computation.

The five saved primitive contributors, indexed0,7,9,10,15, satisfy the literal
polynomial identity sum Ai*primitive_i=1. Independently the five actual-minor
contributors satisfy

    sum Bi*fullMinor_i = E,
    E=z^21*(z+20)^93*(z+60)^6*(z+12).

All coefficients are checked exactly by FLINT, including the long rational
strings. E is a unit throughout the fixed P100 chart. At every geometric
point at least one of the five cleared8-minors is therefore nonzero. All row
LCD factors are chart units, so these minors witness eight actual independent
pure directions. Adding the accepted117 quadratic,23 H1 and one H2 direction
proves original pointwise quartic rank at least149 on the whole P100 chart.

No single original149-minor, column completion, full239 rank, global FULLU
bound, excluded physical point or generic-exhaustion result is asserted.
There is no producer replay, GCD, Gaussian/rank search, unknown-factor
inversion, new source/NF query, point, prime or GB.

The independent executable and its reproducible input are in ignored work.
The committed review contains only small source/hash/arithmetic receipts,
without duplicating the18 MB native JSON. Bounds are60s/2GB/32MiB and150000
live polynomial slots including the conservative44870-slot literal reserve.

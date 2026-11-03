# Exact obstruction to raw h-target cofactors through degree four

Use the literal QQ24-by56 raw four-quadric matrix M and h=a^5*(5P−6e−4)^3
from the previous degree-at-most3 bounded-lift packet. The single new bounded
coefficient-degree4 test has8736 coefficient rows,1960 unknown cofactor
coefficients and24 added coordinate targets h e_i. It includes every monomial
of total degree at most4 in each of the56 unknown entries, and every output
component/monomial of degree at most11. Source denominators are32003-units.

Exact FLINT RREF overGF(32003) gives rank1960 for the source coefficient matrix
and rank1984 for its augmentation by all24 targets: the ENTIRE augmented matrix
has full column rank. Its columns are consequently independent overQQ. Therefore
no nonzero constant QQ combination of h e_i has a polynomial cofactor vector of
total degree at most4. All individual coordinate targets are excluded. This
extends the earlier finite-degree obstruction; it does not imply that higher
cofactors fail, any pointwise rank defect, or failure of QQ localized surjectivity.

The test takes5.70seconds with an RSS upper bound222432KiB, under a60-second /
2200000-KiB cap. No degree5 or larger test is included. The coefficient matrix is
rebuilt directly from the complete literal source. It uses the same coefficient
indexing and exact modular backend as the previous positive/negative diagonal
fixtures, which remain saved there. Source, runner and rank receipts are retained.
This packet awaits independent review.

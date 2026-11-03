# Direct raw four-quadric product criterion over GF(32003)

The source is the four actual cleared quadratic rows in the frozen coordinate
report S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62. Substitute d=P/a
and clear only row powers a^[0,2,2,1]. Form all six light-variable products of
each quadric, using all 56 cubic monomials in the six light variables. The raw
matrix is 24-by-56 with 4824 coefficient terms (804 before taking six products).
Its core degree is at most seven and its QQ coefficient heights are at most
20 numerator / 17 denominator bits. Every rational coefficient denominator is
a unit at 32003. There is no preceding Schur normalization or prefix-integrality
assumption in this raw finite-field computation.

A single Singular process computes its field standard basis with 160 columns,
checks that every original generator reduces to zero, computes the annihilator
of the entire 24-dimensional cokernel, and saturates that annihilator by exactly
U=a*P*F*B0*C*J, with F=500a+90e²+3e−24,
B0=(50a+18e²−9e−2)P+(60e−40)a, C=100a+18e²−9e−2, J=5P−6e−4.
The annihilator has seven generators and dimension two; its U-saturation is
(1). Therefore this raw finite-field matrix is surjective over GF(32003)[U^-1].
The run takes 13.71 seconds, maximum RSS 67516 KiB, under a 60-second /
2200000-KiB cap. The full native basis, annihilator and saturation are saved with
explicit entry-by-entry exports, together with the complete literal QQ input.

This is the STRONGER four-quadric criterion, not the earlier actual twenty-row
residual. Its QQ analogue would suffice for original product rank 41 by the
separately proved seventeen kernel directions, subject to checking the actual
source coordinate and unit bridge. That QQ statement has NOT been inferred
from this finite-field saturation. No QQ rank, full cubic, quartic or ambient
claim is made here. A separately bounded direct QQ calculation uses the literal
Fraction input in a different directory and will be frozen separately.

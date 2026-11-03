# Earlier safe21 actual product module: bounded homogeneous pilot

This distinct pilot replaces the later, coefficient-heavy 31-pivot residual by
an earlier actual residual of the SAME QQ product-rank-41 problem. The saved
36-by-70 polynomial matrix has twenty product rows and sixteen genuine cubic
rows. Its first twenty rows have exactly seven zero columns; deleting only those
columns gives the actual 20-by-63 product matrix. The original row indices and
row powers of a used to clear the P=ad substitution are recorded. The source is
plain native format: its header `36 70 36` is followed immediately by entries,
with no metadata lines. The third header integer is not a metadata flag.

Over QQ on the known unit chart, the preceding twenty-one accepted unit pivots
make pointwise full rank twenty of this residual sufficient for original product
rank forty-one. This QQ implication is not transferred to positive characteristic
without checking all preceding pivot constants. Neither this implication nor a
localized finite-field unit assertion is used in the pilot's conclusion.

The twenty source rows contain 30731 terms, total degree at most fourteen, and
QQ numerator/denominator maxima of 36/47 bits. The later ten-row source contains
205165 terms, degree at most twenty-six and maxima 217/222 bits. Each entire
source column is homogenized to its fixed QQ maximum degree BEFORE reduction.
All row shifts are zero, the variable order is (a,e,P,w), and Singular uses
`(dp,C)`. The field is GF(2147483647); all source denominators are units.

The single authorized run completes in 18.69 seconds, maximum RSS 72720 KiB,
with a monic homogeneous reduced basis of shape 20-by-1050, 532491 terms and
maximum degree twenty-one. All original homogeneous generators reduce to zero.
A separate native checker verifies monic leading coefficients, column homogeneity
and full reduction with 28453920 monomial-divisibility comparisons. The full
finite-field basis is preserved, not a truncation or a serialization of one row.
The cap was 120 seconds / 2200000 KiB, with one Singular process and zero workers.

This is a finite-field computational result and an input-conditioning comparison.
One prime gives no QQ reconstruction, no QQ product-rank-41 theorem, no full
cubic-rank-57 or quartic/ambient conclusion. The complete QQ GB, original-module
containment, p-integral monic coefficients and same-good-prime graded equality
checks described in the homogeneous modular certificate remain necessary before
any QQ membership test is accepted. No additional prime batch is part of this
immutable pilot packet.

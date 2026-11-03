# A33rd product unit pivot by transporting the original two-row certificate

2026-11-03. This continues the frozen scalar-product-row packet, on
`aPFB0CJ≠0`, P=ad. It proves a concrete rank-preserving reduction to24×58,
not all-point product rank41, full cubic rank57, quartic65 or ambient removal.

The preceding ORIGINAL rows0,7 had a60×2 polynomial right inverse up to a.
After its preferred row7/column43 replacement, transport column0 of this
right inverse through the inverse column change. Its polynomial59-entry
vector has138 terms. Its target is row6 of the25×59 residual, which is the
old row0 after row interchange and elimination. The verifier reconstructs
this transport from the preceding frozen certificate.

It also reconstructs every target-row entry: take
`a*oldrow0 − (oldrow0·cofactor7)*oldrow7`, omit column43, and divide by the
actual saved row-cleaning unit. That unit is a nonzero QQ constant times a.
All59 complete polynomial identities are checked with Fraction arithmetic.
The transported vector therefore has target a times a nonzero QQ constant,
a declared unit. It is not inferred merely from a numerical rank or a jet.

None of these cofactors is individually a declared unit. However an exact
QQ constant combination of components50,53,54,55,57 equals1. If alpha is
that constant vector, replace source columns by

    A'_50 = A_50/alpha_50,
    A'_j = A_j − (alpha_j/alpha_50) A_50  (j≠50).

The inverse change acts on cofactors by `c'_50=sum alpha_j c_j=1`, with
other cofactors unchanged. This is an invertible QQ column operation.
The native implementation calculates the FULL polynomial image before
and after the change, asserting equality in all25 rows. Thus the exposed
coefficient permits column replacement and one unit pivot, without a new
redundant column. The result is24×58 with2,755,382 terms,33 total constructive
product-row pivots. The remaining eight product rows and sixteen genuine
cubic rows still require work.

The generic utility's optional eighth argument is a constant-combination
file: column count followed by that many rational constants. The seven
previous arguments match the frozen one-witness helper. A fresh2×2 fixture
with cofactors `(e,e+1)` and constants `(-1,1)` verifies the hidden-unit branch;
its residual is `(e)`. Only named core units and nonzero QQ constants are
inverted. The large generated24×58 cache is omitted from repository storage;
its SHA256, byte count and replay data are retained.

The exact acceptance verifier reconstructs the target row and transported
vector from the preceding original two-row certificate, checks the constant
combination1 and the scalar unit image, and pins the native output.

    .venv/bin/python code/s61_binary_nonsquare_second_scalar_pivot_20261103_v62.py --output work/second-scalar-pivot.json

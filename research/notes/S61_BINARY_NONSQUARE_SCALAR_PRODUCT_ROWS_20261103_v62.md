# Two scalar unit lifts in the primary full cubic residual

2026-11-03. Work on adFBCJ≠0, J=5ad−6e−4. The serialized core is (a,e,P=ad),
with B0=aB. This packet continues the frozenV9 primary26×60 matrix after31
constructive product-row unit pivots. It proves exact further algebraic reductions;
it does not prove product rank41, full cubic rank57, quartic65 or ambient removal.

## Original polynomial identities

The coefficient matrix is pinned by SHA256
`6d21ebfbeade3c5a0a7bd2dc091bde5ec274622ddf899add905242164d616ea3`.
Its first ten rows are remaining genuine quadratic-linear products, followed by
sixteen genuine cubic rows. The following refer to zero-based rows in this matrix.

The entries of row0 generate an ideal whose reduced QQ Gröbner basis is

    (a, eP², 3e²P+2eP, 9e³−4e).

Row7 similarly has reduced basis

    (a, P², 3eP+2P, 9e²−4).

No per-entry unit stripping occurs in either of these two rows: all of their
original entries have zero valuations for the named units. The scalar Gröbner
bases discover the unit a. The final certificates are ordinary polynomial vectors:

    row0 L0 = a,      row7 L7 = a.

L0 has64 terms, coefficient degree4 and maximum256 coefficient bits.
L7 has34 terms, coefficient degree3 and maximum233 bits. Both complete identities
are independently replayed with Python Fraction arithmetic, and the certificate
rows are compared to the entire pinned original matrix. No numerical extrapolation
or assumption that a computed list is a Gröbner basis is needed for acceptance.

## Preferred explicit32nd pivot

L7 has a nonzero QQ constant in component43. Replacing source column43 by the
full image of L7 is therefore an invertible column operation: its determinant
is this nonzero constant. Its row7 entry is a, a known unit. One exact native
Schur pivot gives25×59, with1,801,170 terms and32 total constructive product-row
unit pivots. All26 original image entries are computed over QQ and logged.
The output is about158MiB and is excluded from the compact saved packet; its
SHA256, byte count and exact replay recipe are preserved.

For comparison, L0 itself has no individually known-unit component. Appending
its image and pivoting gives25×60. Two P-linear components have a nonzero QQ
constant in their constant span; an exact invertible two-column QQ change
exposes a unit component and permits25×59 instead. This is valid but produces
2,101,063 terms. The L7 route has better conditioning and uses the already
accepted native helper without this extension.

## A stronger two-row certificate

For the ORIGINAL2×60 submatrix M of rows0 and7, a separate60×2 polynomial
matrix T with206 terms satisfies

    M T = a identity2.

All four identities are checked independently with Fraction arithmetic. Thus
these two coordinate rows define a split surjection on the named open. This
packet has concretely executed only the preferred first scalar pivot, so it
records32 constructive native pivots. The two-row certificate can be transported
through that pivot to find the next one without another unit search.
No freeness of the full ambient cokernel follows from the split projection.

Per-entry unit stripping is safe for a scalar entry ideal after localization.
It is not valid independently on entries of a multirow matrix. An exploratory
three-row primitive matrix was discarded for precisely this reason; it supplies
no claim here. The two accepted original rows had no stripped factors.

## Replay and storage

`certificate.json.gz` contains both original source rows, both scalar vectors,
the full two-row right inverse and native output hashes. Reconstruct the pinned
source with the immutable V6/V9 recipes, then run

    .venv/bin/python code/s61_binary_nonsquare_scalar_product_rows_20261103_v62.py --output work/scalar-product-rows.json

The native preferred step uses the frozenV9 `witness-image-stats.cpp` interface:

    witness-stats-bin FULL26 OUTPUT P-units 1 ROW7_WITNESS 7

`ROW7_WITNESS` is the60-polynomial scalar vector in the certificate, in the existing
sparse native format. Success must log `EXACT_WITNESS_REPLACE 43 TARGET 7`,
`TERMS 1801170` and `PASS 1`. Existing setup/runtime is reused.

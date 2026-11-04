# Source-only finite-H15 border with the new C68W source

This preflight changes only the bottom source of the previously fixed
case9/C-plus28 border. The failed K2W minor is not recomputed. The top eight
sources and their known-unit LCDs are unchanged; the new bottom row is the
saved literal normalized C68*v3, original319 row226.

Physical v3 is exactly joint W in the accepted frame. The new row's complete
98-component denominator profile uses only the five known curve units.
Its full-row LCD is saved, and every selected entry denominator divides
the corresponding old/new LCD. The nine old local98 columns remain
`[0,1,2,3,4,5,6,9,28]` with their original122 labels saved explicitly.

The support pattern has a perfect matching. Inputs total1057 RF terms,
with1508 terms as the cleared-entry bound. The determinant degree is at most
218, giving219 polynomial terms and162936 bytes as a conservative raw JSON
bound. All input QQ coefficient denominators pass the old-prime32003 gate.
These are cost/source facts only; no determinant or finite remainder is
evaluated.

The proposed next job is one global QQ[z] Bareiss determinant, a literal
remainder modulo the exact monic H15, and gcd/Bezout of the **whole** finite
special fiber `(Hbar,Rbar)` at32003. It never inverts an old8-minor or
finite-ring matrix pivot. Monic H anchors finiteness, so no leading-coefficient
condition on the new determinant is required. A unit result would remove the
old finite chosen-cover exception and prove whole150; a nonunit result is
inconclusive and authorizes no automatic change of source, column or prime.

Code/report stem:
`S61_NONSQUARE_P100_K7W_FINITE15_PREFLIGHT_20261104_v63`.
The source-only pass took0.038 seconds. All bytes are immutable; no new9-minor,
QQ/Fp gcd, source query, rank or factor inventory is performed.

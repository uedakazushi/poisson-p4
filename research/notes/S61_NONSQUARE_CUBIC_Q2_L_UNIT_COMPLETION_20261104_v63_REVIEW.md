# Independent review: actual Q2 L-unit completion

The preceding GL6 completion has an independent accepted review. This audit
checks the next frozen GL3 packet and its three-minor frontier without
repeating any Gaussian elimination. In the actual 16-by-23 presentation,
rows labelled 5 and 9 at old columns 21,22,25 are literal checkpoint leaves.
The old GL6 operations only changed no-L coordinates and cleared with a row
having zero L coefficients. All six leaf identities are checked directly.

Over QQ(a,e,d), the actual row b=(b0,b1,b2) satisfies
T=900a-126e-121, b0=aT/150 and
2b1+b2=(18e+7)T/1125-1/9. Its dot with
(6(18e+7)/(5a),-18,-9) is exactly one. The only new nonconstant denominator
is the already inverted a. A separate free-polynomial engine checks both
inverse products for U0, E and U0E, the full row formula
(t,(1-t)b0,(1-t)b1), and determinant -9. No author helper is invoked.

The old 2,138 scalar nodes are unchanged. All 3,914 resulting nodes are
independently replayed as commutative polynomial expressions over abstract
checkpoint and pair-dual leaves. This binds every GL3 entry, the embedded
23-target operators and inverses, all fifteen tail clearings, the elementary
16-source operator, and all 330 entries of the 15-by-22 answer. None of the
five heavy/psi targets 27..31 is dropped or changed by the target operator.

Substituting the proved actual dot makes the entire pivot row the unit vector
at old column 21. Clearing this column from the other sources has determinant
one. Deleting this unit source and target preserves the cokernel, so the
actual 16-by-23 presentation becomes the displayed transformed 15-by-22
presentation. It is not asserted to be a literal old submatrix. There are ten
proved global unit pivots in this transformed projected presentation.

For completeness the audit also reconstructs exactly the three previously
selected L-block determinants and their allowed-factor reconstructions.
They are nonzero and their residual polynomials are not individual named
units. This is only a fixed-minor frontier: it proves neither rank drop nor
properness of the ideal they generate. The positive row Bezout identity is
compatible with that frontier.

No further pivot, unknown inversion, point, prime, GB, full-FULLU cubic57,
quartic239, or kernel-count change is claimed. Source and proof hashes are
recorded in the dedicated independent review report.

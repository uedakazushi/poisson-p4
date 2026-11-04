# Actual residual Q3: support and bounded candidate costs

This read-only preflight uses the accepted eight-step checkpoint of
`S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63`. It does not evaluate a
minor, dot product, new point or Gaussian step. The actual source rows
are 4, 5 and 9, labelled `scaled_q5*x0`, `scaled_q5*y1` and `scaled_q5*L`.
These labels identify the saved transformed relations, not untouched
original rows.

The exact three-by-24 matrix is saved in the matching certificate. It
has 9, 19 and 17 nonzero entries respectively. All coefficient
denominators are rational constants; numerator degrees are at most
4, 6 and 9. The union of supports has 19 pure target coordinates.
The five remaining heavy/psi coordinates are zero on all three rows.
There are 16 active no-L columns and three L columns, with the latter
numbered 21, 22 and 25. The first row has no L terms. The other two
have respectively 15 and 43 numerator terms in those three columns.

The computation only counts potentially nonzero determinant terms.
Four lowest-cost fixed candidates are columns `(2,21,22)`, `(2,21,25)`,
`(7,21,22)` and `(7,21,25)`. Each has only two potentially nonzero
permutation products, raw multiplication support bound 315 and degree
bound 11. No value, factorization or nonvanishing conclusion is
attached to those candidates.

A more constructive observation supersedes the need to start with
those minors. Row 4 was modified only by the old `q0*x0` pivot. Its
six columns 4 through 9 are `x0*Sym2(y1,y2,y3)` and equal the accepted
scaled-q5 pair minus a multiple of the q0 pair. The full second column
of the accepted pair dual is therefore a candidate unit functional
on these six coefficients. It has a literal zero in local slot 1.
The preflight saves this candidate but does not evaluate its dot.
The separately authorized completion packet verifies that identity
and constructs the finite unit change of coordinates.

The global rank-three result for these rows comes from the accepted
rank-19 Z=W quadratic-product theorem and the seven saved quadratic
pivots. It does not give a computed target basis or an explicit
right inverse. No whole-FULLU cubic-rank-57 or full-quartic claim
is made here.

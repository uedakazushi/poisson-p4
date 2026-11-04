# Nine small coefficients of the last quadratic relation

This packet queries exactly nine coefficients of the actual remaining
quadratic source row labelled 9 in the global `15x22` presentation.
It does not expand the other coefficients of that row or another
source row. All nine coefficients are independently evaluated from
their 114 scalar-DAG ancestors and matched to the source formulas
below. The calculation took 0.193 seconds.

In the old eight-step checkpoint, let `b` and `s` be rows 5 and 9
on original L-target columns `(21,22,25)`, and put

`v=(6(18e+7)/(5a), -18, -9)`, `c=s*v`.

The exact dot `c` simplifies to a polynomial with nine terms. The
two preceding target changes involve only x0-containing targets
and the L-target triple. The first pivot row has zero coefficients
on all pure-y-cubic targets. Consequently the seven final
coefficients on original columns

`[11,13,14,15,17,18,19]`

are exactly `old_row9_j-c*old_row5_j`. The last two queried
coefficients, on final target labels 22 and 25, are respectively
`s[0]-c*b[0]` and `s[1]-c*b[1]`. The shifted L-target indices
are intentional: they are the first two old coordinates after the
GL3 completion, with its pivot coordinate 21 removed.

The raw multiplication support bound is 5,890, below the authorized
6,000-term gate. The nine exact outputs have a total of only 373
numerator terms, degrees at most nine, and rational constant
denominators. The seven pure coefficients have d-degree at most
two. The two L coefficients have d-degree one and have 14 and 23
numerator terms. Thus a proposed simultaneous affine-linear
elimination in d for all nine coefficients would be incorrect.

The certificate saves all nine polynomials, the input triples and
the nine-term scalar `c`, the exact DAG roots and their allowed
checkpoint leaves. It also saves the full source-cofactor recipe
relative to the old eight-step source frame:

`row9 - c*row5 + (c*d5-d9)*row4`.

Here `d5,d9` are the scalar column coefficients in the preceding
GL6 source clearing. Their accepted finite-DAG roots remain lazy;
they are not discarded or expanded. The full invertible target
changes are still those in the two preceding completion packets,
so pure and heavy tails remain part of the presentation.

Only powers of the named FULLU units `a,d,F,B,C,J` are stripped.
The first eight numerators have one factor `a`; the final coefficient
has no such named factor. None of the nine residual polynomials
is a nonzero constant. All factorizations are checked by exact
reconstruction. This says that these entries are not individually
recognized units. It does not say their ideal is proper or that
the actual row vanishes anywhere.

No Gaussian, determinant, new point, GB, unknown factor inversion,
extra coefficient query or automatic next operation was performed.
An explicit Bezout identity and target completion for the remaining
row are still open. Whole-FULLU cubic57, quartic239 and the quartic
kernel count38 with gap42 are unchanged.

# Fixed pure-functional tests and the affine annihilator obstruction

The authorized constant and affine QQ functional tests are complete.
For each of the correctly source-bound rows G18, G39 and G68, none
of the eight fixed known-unit targets is attained. The calculation
took 2.542 seconds and its exact certificate is 1.8 MB. No higher
degree, point, GB, new minor or source/target transport was tried.

The input is the preceding reviewed source-frame preflight:
`Q10*w=0`, `Gi_clear*w=C*u`, with targets

`1,a,d,F,B,C,J,a*d`.

Constant functions have 26 rational unknowns. The fixed affine
cofactor basis `(1,a,e,d)` has 104. The three constant combined
systems have respectively 169,169,152 coefficient rows and the
affine systems have 338,338,311. Their exact column ranks are
26 and 104 in all three cases. Each saved span basis is reconstructed
from its input-column witness and every source column reduces to
zero. All 48 target reductions have nonzero exact remainders.

## A stronger conclusion already contained in the same certificate

Every normalized pivot in all six saved span bases belongs to one
of the ten quadratic constraints, whose equation labels are 0
through 9. None belongs to the genuine equation, labelled 10.
The basis vectors have triangular identity coefficients at their
successive pivot rows. Their projections to the Q10 constraint
rows are therefore already independent. It follows that the Q10
coefficient matrices alone have full column ranks 26 and 104.

Consequently there is no nonzero polynomial tuple w of core
degree zero or one satisfying `Q10*w=0`. This is stronger than
failure of the eight target values on these three sources: in
these two coefficient classes, there is no nonzero functional
annihilating the quadratic image to begin with. This conclusion
uses the recorded pivots and witnesses, not an additional rank
computation. A common unit denominator cannot change it, because
that denominator clears from the homogeneous equation `Q10*w=0`.

This is a finite polynomial-degree statement. It does not say that
the quadratic image lacks annihilators over the parameter field,
or rule out higher-degree polynomial or regular functional tuples.
No degree-two or larger test is authorized by, or automatically
launched from, this negative result.

## The genuine H*L lift remains part of the correct candidate recipe

The input preflight proves that all six projected-heavy coefficients
of these three sources vanish before the G11 elimination. Their
actual G11 source-cleanup coefficients are zero. A future candidate
pure w nevertheless must annihilate the deleted mixed G11 row
when lifted to the earlier presentation. Its only added component
is the deleted `H*L` coordinate

`- dot(G11_pure,w) / G11_HL`.

The pivot `G11_HL=2*F*B²*C²/15625` is checked directly against the
original source. The denominator is a known FULLU unit, and this
one-component lift leaves both Q10 annihilation and the selected
Gi value unchanged. The code contains the all-QQ coefficient
checks for such a candidate, including the lambda/C source identities.
There were no candidates, so no lift, unit normalization or pivot
was executed; the candidate-verification report flags are vacuous
over an explicitly empty candidate set.

All pure and heavy tails remain in the actual presentation. The
global `14x21` same-cokernel theorem and all accepted kernel/row
counts are unchanged. Neither global cubic57 nor quartic239 is
proved or disproved by this test.

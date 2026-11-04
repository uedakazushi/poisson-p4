# Compact input for a later P100 pure-rank test

This source-only receipt does not execute a field rank or minor calculation.
The already saved eight-by98 block consists of the actual q5H-cleared pure
row followed by the seven literal G18-times-light rows. It has644 nonzero
entries,18193 numerator/denominator terms and at most214-bit coefficients.

The explicit GL2 unit split proves, over QQ(z), that its rank is one plus
the rank of the seven-by97 residual block. Thus this compact original
block can replace the22.6MiB Bezout-derived rows in a future generic-rank
test, without changing that question.

Only existing whole-curve unit LCD profiles are used. Clearing each row's
five-factor LCD has23618 polynomial terms as a conservative total bound;
no clearing polynomial is evaluated in this receipt. This gives a small
input for one separately authorized exact QQ(z) row-rank/selected-minor
transcript. Unknown polynomial pivots would certify generic rank only.
Their possible finite zeros require a separate whole-curve invariant-factor
or determinantal-ideal argument. No rank or new lower bound is inferred here.

Code: `code/s61_nonsquare_P100_pure8_rank_source_preflight_20261104_v63.py`.
Evidence: `reports/S61_NONSQUARE_P100_PURE8_RANK_SOURCE_PREFLIGHT_20261104_v63`.
The source-only pass took0.091 seconds; all files are immutable.

# Unexecuted optional-row index clarification

2026-11-04. The frozen `P100_H1_DUAL` computation correctly stopped at the
NF117 descent failure. Its sixteen raw-H annihilations, six original
quadratic values and q1 unit value are unaffected.

Its unexecuted optional q5H branch had an incorrect cache-row label. The
controlling producer explicitly defines

```
complete134 = [q5H] + cubic112 + quartic21.
```

Therefore q5H is row **0**, cubic products are rows **1 through 112**, and
quartics are rows **113 through 133**. Row112 is the final cubic product.
The old optional q5H row112 input was parsed but never evaluated because
the descent gate failed; no optional scalar result was claimed. Quartic
row118 remains correctly bound to genuine source label55.

The delta7 source preconditioner rows43 through49 and rows8 through14 are
correct: they correspond to cubic genuine indices6 and1 respectively,
with this initial q5H offset of one. New right-inverse code must use
`cache_row=1+7*genuine_index+light_position` and preserve that actual source
preconditioner. Frozen producer and result bytes remain unchanged. The
old optional branch must not be generalized or executed using its incorrect
q5H label; a new producer with the corrected index is required for that use.

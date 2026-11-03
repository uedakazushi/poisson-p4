# J=0 cubic cokernel: exact preparation, unfinished PID continuation

The actual 57-row cubic basis uses the alternative ad basis, followed
by 21 a,F,B,C unit pivots. The saved 36x70 source never inverts J.
Its hash matches the frozenv6 handoff manifest exactly.

ExactJ0 substitution d=(6e+4)/(5a) clears each row onlybya^M.  An independent
binomial implementation compared ALL2520QQ polynomial entries (108530terms).
An additional13known-unit pivots gives23x57 over
QQ[a,e,1/(a*F*Bs*C*(3e+2))].  That regenerated37.9MB output remains inwork;
its hash/reproduction recipe is retained rather than duplicating it here.

Both traced coordinate-field Smith runs timed out at180CPU/190wallseconds
without a Smith form, annihilator or unit result.  A built-in rational-field
extgcd failure emitted an INVALID nonunit candidate; it is explicitly excluded
from all mathematical claims.  Explicit verified Euclidean arithmetic passes
a nontrivial rational-coefficient fixture and its saved proof replay.

The valid QQ(e)[a] localizedrow-gcd run saved12cofactor steps, but conservative
scope is11completed plus12thwitness saved.  The restored60CPU/65wall prefix
replay also timed out before final confirmation/export.  No finiteness,
all-point rank57, actual defect curve or ambient nilpotence is proved here.

`continuation.tar.gz` restores curated files at their originalwork paths;
compiled binaries are omitted.  `report.json` fixes every member hash and
clearly records the unsuccessful computations.  Restore at repository root:

```sh
tar -xzf reports/S61_NONSQUARE_J0_PID_20261103_v63/continuation.tar.gz
```

Proof/frontier: research/notes/S61_NONSQUARE_J0_PID_20261103_v63.md.
Next: replay the SAVED cofactors with per-step progress, export the24x70
state, then continue under a new justified bounded budget.  Do not repeat
the first12gcd searches or unchanged tracedSmith computations.

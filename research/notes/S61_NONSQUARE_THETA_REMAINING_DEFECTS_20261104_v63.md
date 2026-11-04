# The last Theta defect supplies a generic extra direction

Date: 2026-11-04. Exact QQ primary certificate; independent review pending.
This packet evaluates only the two remaining authorized fixed defects.

The source B rows are ordered `[7,9,10,11,12,13,14]`; its recorded pivot
rows are `[7,13,10]`. Before using the rational Theta chart, all six full
QQ[a,e,d] coefficient identities for raw rows 9 and 14 are verified:

```
Delta*B_i=sum_j lambda_ij*B_pivot_j,  i=raw9,raw14.
```

The actual seven-by-four Theta matrix remains `[B|S]`, with S the B^4-cleared
Schur column after the raw8 unit pivot. Native exact substitution evaluates
only S_raw9, S_raw14, S_raw13, and the one necessary old clearing coefficient
V. The already frozen S_raw7 and S_raw10 are reused. New term counts are
157,157,127,8 respectively. In particular V has a nonzero image in the
actual rational Theta domain; it is not declared a unit.

The fixed defects are

```
D_i=Delta*S_i-sum_j lambda_ij*S_pivot_j.
```

All four product polynomials for each tested row are saved. The exact result
is

```
D_raw9=0,
D_raw14 != 0, with 269 Laurent polynomial terms.
```

The latter is retained as a completely unfactored unknown polynomial. No
global inverse of Delta or D_raw14 is introduced. There is no automatic test
of another row, rank/minor family, GB, point or prime.

The earlier raw12 defect is zero. For raw11 the previous small certificate
gave `V*S11-U*S7=0` and `Delta*U-lambda*V=0`. Therefore

```
V*(Delta*S11-lambda*S7)
 = Delta*(V*S11-U*S7)+(Delta*U-lambda*V)*S7=0.
```

Since the actual Theta ring is a domain and the saved image of V is nonzero,
its Delta-coordinate defect is also zero. Thus exactly three of the four
tested nonpivot defects vanish and the fourth is nonzero. The proposed
alternative conclusion that every four-by-four minor vanishes is
inapplicable.

The nonzero Delta image and nonzero D_raw14 have nonzero product in the
rational Theta domain. Over its fraction field, B has rank three and the
raw14 row gives the extra S direction. Hence the seven-by-four matrix has
generic rank four. The accepted raw8 unit pivot adds one, and the two
earlier FULLU raw-H unit pivots add two. Consequently, on this branch,

```
generic rank(heavy14x5)=5,
generic rank(original raw16x7)=7.
```

There is also a pointwise criterion requiring no Delta inverse: the globally
defined source combination
`kappa=Delta*e_raw14-sum lambda_14j*e_pivot_j` has `kappa*B=0` and
`kappa*S=D_raw14`. Wherever B has rank three and D_raw14 is nonzero, S adds
an independent direction and the heavy block has rank five. The possible
heavy defect is therefore contained in the union of the actual B-defect
locus and `V(D_raw14)`. The unknown polynomial is not yet simplified, and no
claim that this union is empty or finite is made.

This concerns only heavy coefficients. Pure components can compensate for
a heavy defect, and the full cubic57 and quartic239 conclusions remain
open. No negative full-rank conclusion follows.

Producers: `code/s61_nonsquare_Theta_remaining_defects_20261104_v63.py` and
`.cpp`. Evidence: `reports/S61_NONSQUARE_THETA_REMAINING_DEFECTS_20261104_v63`.
The fixed batch completed in 2.271 seconds including compilation, with 3965
stored polynomial terms under 57-second/2-GB/32-MiB/150,000-term guards.
All source, native input/output and prerequisite hashes are pinned.
Certificate SHA-256:
`a70ff73a2d9cffaca13949856a181cdbcce61e54c8ca62bd8328ca01effc0e5f`.

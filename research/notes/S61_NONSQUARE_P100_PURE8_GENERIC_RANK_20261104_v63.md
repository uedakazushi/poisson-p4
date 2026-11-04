# An exact generic8-minor in the P100 pure source block

Date: 2026-11-04. Exact characteristic-zero primary certificate; independent
review pending. This uses only the eight already saved pure98 rows. No new
source/oracle value or further genuine cubic is evaluated.

The first row is the actual q5H-cleared pure source; the other seven are
G18 times physical lights `[v0,v1,v2,v3,v4,v6,v7]`. Clear each row's known
five-factor P100 unit LCD. All784 coefficient identities between the saved
RF rows and these polynomial rows are checked. The complete input remains
the compact eight-by98 block; it is never fraction-free Schur-eliminated.

## Incremental minor and certificate

Keep only a selected square S, its determinant D and adjugate H. For the
next candidate border `[S b; c d]`, compute

```
h=D*d-c*H*b.
```

If h is nonzero as an exact QQ[z] polynomial, append that row and column.
The new adjugate is

```
[[ (h*H+(H*b)*(c*H))/D, -H*b ],
 [ -c*H,                    D ]].
```

Every division displayed here is checked as an exact polynomial identity.
Both left and right adjugate products are checked at every stage. Starting
with the empty determinant1 proves inductively that each saved D is the
literal selected determinant with the saved row/column order and sign.
No selected unknown polynomial is declared a whole-curve unit.

The algorithm reaches eight rows. The selected source positions are
`[0,2,3,5,7,6,1,4]`; selected local98 target positions are
`[0,1,2,3,4,5,6,9]`, corresponding to original122 labels
`[0,1,2,4,5,6,7,10]`. The final nonzero polynomial determinant, selected
eight-by8 matrix and final adjugate are saved. Since the original matrix
has eight rows, this proves its exact rank8 over QQ(z). The accepted GL2
unit-split transport then gives exact rank7 for its seven-by97 residual.

## Scope in the actual quartic map

The selected polynomial determinant differs from the selected RF minor
only by the product of eight known whole-curve row LCD units. Therefore
the actual eight pure sources are independent at every P100 point in
D(det). Together with the accepted117 quadratic,23 H1 and one H2
directions, the original weighted quartic map has pointwise rank at least
**149 on D(det)**, a nonempty open of P100.

The unknown zeros of det have not been excluded or tested. The independently
accepted whole-P100 bound remains142. There is no uniform149, pure98,
full239, other-curve, allFULLU or automatic next-job claim.

## Resources and files

The one job took5.879 seconds including compilation; native arithmetic
took0.787 seconds. It checked19 bordered candidates,140 exact divisions,
all eight two-sided adjugate identities and all784 row clearings. The
instrumented peak was58141 polynomial terms including literal inputs;
10957 terms were serialized. The largest temporary polynomial length365
includes pre-division adjugate-update numerators.

Producers: `code/s61_nonsquare_P100_pure8_generic_rank_20261104_v63.py` and
`.cpp`. Evidence: `reports/S61_NONSQUARE_P100_PURE8_GENERIC_RANK_20261104_v63`,
including `native.json`, `source-labels.json`, selected source/target bindings,
the native input and report. All primary bytes are immutable. No full
eight-by98 elimination, new point, prime, GB or automatic follow-up ran.

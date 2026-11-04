# Three pure source rows: correctly bound functional input

This preflight prepares a finite coefficient-linear test for a pure
target functional. It performs 78 known-unit divisions and records
source cleanup and input sizes. It does not solve a QQ span problem,
choose a point, run a GB or expand any POWER coefficient.

The sources have original augmented labels `10,13,16`, relative
genuine labels `0,3,7`, and original cubic-equation labels `18,39,68`.
Their initial 26 pure coefficients have respectively
`1059,1055,794` terms. The frozen ten quadratic-product rows have
317 terms.

## The genuine H*L pivot is included in the audit

The code checks all six projected-heavy coefficients of each source
**in the initial augmented matrix**, before the old genuine pivot.
They are all literally zero there, and remain zero in the old
checkpoint. This is a statement about the six projected-heavy
coordinates of the augmented matrix, not an assertion that all
seven original raw H coefficients vanish.

The old genuine pivot is source label 11, target label 26=`H*L`.
The three sources have no entry in that pivot column, and none
appears in that step's elimination list. Their exact G11 source
subtraction coefficients are therefore zero. No mixed genuine
source correction is silently discarded.

Every earlier change to these three rows is a multiple of one of
the quadratic-product rows. The code replays only that fixed source
cofactor recipe as a finite DAG: normalized quadratic pivots and
saved pre-elimination coefficients. It saves the ten quadratic
cofactor roots for each source, with its own original source
coefficient equal to one and both other genuine coefficients zero.
Thus if `Q10*w=0`, the old checkpoint value on the source equals
its initial pure value. No later degree207 target completion is
needed to formulate this test.

## Correct known-unit normalization and clearing

All three initial source-clearing factors are exactly
`lambda=-64000*C`. Each of the 78 divisions by lambda is checked
as a rational-function identity. The resulting coefficients are
not all polynomials: their denominators can contain one factor C.
The preflight retains these denominators explicitly rather than
assuming exact polynomial division.

Write `Gi_norm=Gi_initial/lambda`. The proposed equation
`Gi_norm*w=u` is cleared correctly as

`Gi_clear*w=C*u`, where `Gi_clear=C*Gi_norm`.

The cleared source term counts are again `1059,1055,794`; the
normalization does not justify claiming a smaller polynomial input.
A success would give the original source value `lambda*u`, a
known FULLU unit, provided that the simultaneous quadratic
annihilation constraints are satisfied.

## Fixed finite linear-algebra classes, not yet solved

The target vector w has 26 entries. Constant QQ entries give 26
unknowns; the fixed affine basis `(1,a,e,d)` gives 104. In each
case impose every coefficient of all ten equations `Q10*w=0`
and one correctly cleared source equation. The fixed target
units are

`1,a,d,F,B,C,J,a*d`.

The exact combined row and serialized input counts, using a common
coefficient basis for all eight targets, are:

| Source | Constant rows / unknowns | Constant terms | Affine rows / unknowns | Affine terms |
|---|---:|---:|---:|---:|
| G18 | 169 / 26 | 1435 | 338 / 104 | 5563 |
| G39 | 169 / 26 | 1431 | 338 / 104 | 5547 |
| G68 | 152 / 26 | 1170 | 311 / 104 | 4503 |

These counts include the known C factor on each target. No matrix
rank, nullspace or membership solve has been performed. A future
candidate must be checked against all original QQ coefficients,
including all ten quadratic annihilation equations. Testing a raw
genuine coefficient alone would not establish independence modulo
the quadratic image.

The actual source and pure/heavy tails remain intact. A successful
functional can be transported through the accepted invertible
operations without expanding their large core coefficients. The
present packet claims neither a new unit pivot nor cubic57,
quartic239, a kernel count change or failure of any higher-degree
functional class. The fixed preflight took 0.503 seconds.

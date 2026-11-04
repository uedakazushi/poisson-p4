# The fixed raw12 Schur-column defect is zero

Date: 2026-11-04. Exact QQ primary certificate; independent review pending.
Only the one authorized raw12 expression is evaluated here.

In the generic B7-by-3 source certificate, selected rows are block indices
`[0,5,2]`, or raw source rows `[7,13,10]`. For block row 4 (raw12), its
coordinates satisfy `lambda_1=0` and `lambda_2=Delta` literally. Before
substitution, this packet independently checks all three full polynomial
identities

```
Delta*B_raw12 = lambda_0*B_raw7 + Delta*B_raw10.
```

They hold in the original QQ[a,e,d] coefficient ring. Delta is the saved
407-term polynomial, not a newly computed minor and not a global unit.

The accepted rational Theta coordinate map is then used to substitute only
the saved cleared Schur entries S of raw rows `[12,7,10]`. Native FLINT
arithmetic uses the literal coordinate numerators and known powers of d.
The resulting Laurent polynomials have term counts `[138,89,138]`.
Using the previously saved transformed Delta and lambda_0 gives the exact
bivariate identity

```
Delta*S_raw12-lambda_0*S_raw7-Delta*S_raw10 = 0.
```

All three product polynomials are saved, so the zero is a full coefficient
identity, not an evaluation or interpolation. The job stores 1506 polynomial
terms and completes in 1.459 seconds including compilation under
57-second/2-GB/32-MiB/150,000-term guards.

Thus this second fixed B-source relation also supplies no extra heavy
direction on Theta. Together with the earlier raw11-vs-raw7 zero, this still
does not determine the rank of the entire 7-by-4 block: the raw9 and raw14
coordinate defects have not been evaluated. No inference concerning the
full cubic or quartic presentation is made. Delta has never been inverted
to make a global claim. The named-unit stripping branch was not entered
because the exact defect is zero; there was no factor search, rank/minor
solver, GB, point or prime.

Producers: `code/s61_nonsquare_Theta_raw12_defect_20261104_v63.py` and `.cpp`.
Evidence: `reports/S61_NONSQUARE_THETA_RAW12_DEFECT_20261104_v63`.
The report pins source matrices, coordinate map, Delta/lambda profiles and
both native input and output. Certificate SHA-256:
`ec78dd6efd20084864fa25c3be549e8c3d774a70b869442e39d366f1ed41592b`.

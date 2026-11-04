# A whole-P100 right inverse for the H1 projection

Date: 2026-11-04. Exact characteristic-zero primary certificate; independent
review pending. This packet gives explicit coefficients from the 112 cubic
products to all 23 H1 targets. It does not invert the full quartic matrix.

Work over the complete rational P100 curve ring

```
C=QQ[z,1/(z(z+20)(z+60)(203z+2700)(z+12))].
```

The actual raw H coefficient rows have rank six here, with the normalized
dual eta of the preceding packet satisfying `eta(v3)=1`. The exact source
quadratic value is

```
u=q1(eta)=137500/(243z^2),   u^-1=243z^2/137500.
```

This is a unit on the entire curve. The raw-H defect therefore need not
produce an H1 defect after imposing the quadratic relations.

## Explicit finite source construction

Put W=v3 and `li=vi-eta_i*W` for the seven original physical light
variables `[v0,v1,v2,v3,v4,v6,v7]`. The coefficient vectors expressing each
li in the sixteen raw H-linear source rows are saved explicitly. They use
the fixed six-row unit basis with original genuine indices
`[1,6,8,9,11,14]`. There is no new Gaussian or minor choice.

The source solve follows the actual unit operations: the eta-functional and
Cdelta row, the unnormalized second checkpoint pivot divided by its known
entry, the Theta row8 unit, and the B adjugate divided by the proven unit
Delta. All initial and later fraction-free source row divisors are included.
The actual transformed selected rows are checked against all thirty
checkpoint entries before the solve. Every one of the seven source lifts
is then checked in all seven physical coefficient slots.

For any two linear forms f,g the identity

```
f*g-eta(f)*eta(g)*W^2
 = (f-eta(f)W)*g + eta(f)*W*(g-eta(g)W)
```

gives explicit cubic-product coefficients. Apply the same identity term by
term to q1, writing its resulting source combination as F(q1). For each
retained target Q, use

```
Q = F(Q)-[eta^2(Q)/u]*F(q1) + [eta^2(Q)/u]*q1.
```

Every coefficient on the right is saved. The target Q is the corresponding
quadratic monomial in the actual joint light coordinates
`(L,ell,x0,x3,r,Z,W)`, transported to the original physical lights by the
pinned forward frame and constant inverse frame. All 23 identities are
checked in every one of the 28 physical quadratic coefficient slots.

These are identities involving the **H-linear parts** of the cubic sources.
Actual cubics also have pure tails. Multiplying by H and replacing H-linear
parts by full source rows adds the corresponding pure quartic tails. Those
tails remain in the lazy pure98 operator; they are neither set to zero nor
silently discarded. No assertion `H*li=0` or `H*W^2=0` in the full ideal is
made.

## Actual NF117 and source-label binding

The five eliminated H-quadratic directions are rebuilt from the full
corrected joint prime rows, including the Z/W corrections. The fixed four
monic pivots and the q0H pivot are normalized using whole-curve units only.
All five original q0H through q4H relations are reduced explicitly to zero
in their 23 retained coefficients. In particular `NF117(q1H)=0`.

The retained 23 monomials match the actual parent target labels. The
normal form has the same five-dimensional source kernel and is the identity
on this retained complement, so it is the actual H-degree-one block of the
parent NF117 operator. Applying it to the 23 source identities gives the
identity matrix on the H1 targets.

The complete 134-row cache is ordered

```
[q5H] + cubic112 + quartic21.
```

Thus the right-inverse rows are actual cache rows 1 through 112, not 0
through 111. The actual delta7 preconditioner replaces cubic genuine6 by
`delta*(C6-C1)`, and the cofactor transport is

```
c_new(6,j)=c_old(6,j)/delta,
c_new(1,j)=c_old(1,j)+c_old(6,j),
delta=-64000*C_core^2.
```

Both source frames are checked directly in all 23 times 28 coefficients,
giving 1288 full QQ(z) coefficient checks. Every nonzero saved cofactor
denominator is verified to be a unit of the complete P100 curve ring.

Consequently the native file's **112-by-23 coefficient matrix is an actual
right inverse for the H1 projection** from these cubic products on the
whole curve. No quartic row was queried; the 21 quartics and pure98 part
remain present in the full source model.

The known 117 quadratic pivots, the separate H2 unit and these 23 H1
directions imply pointwise rank at least 141 for the original full quartic
matrix at every P100 point. Row combinations may be appended without
altering the source image to exhibit this unit block. This does not assert
a single fixed original 141-minor or an additional free source complement.
Full rank239, full cubic57 on the entire curve and all-parameter exhaustion
remain outside this packet.

## Evidence

Producers: `code/s61_nonsquare_P100_H1_right_inverse_20261104_v63.py` and
`.cpp`. Evidence: `reports/S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63`.
`native.json` saves all seven source lifts, the five normalized H-quadratic
rows, all 112-by-23 cofactors and all q1 coefficients. The fixed run took
4.683 seconds, peak univariate polynomial length 31, and 18,019 stored terms
under the 57-second/2-GB/32-MiB/150,000-term guards. There was no full134
Gaussian, POWER expansion, new point, prime, GB or minor search.
Certificate SHA-256:
`609faa1eb59d6d6e3079b0d97b68de58ec45366a3cc4942cf40112dbc1a13549`.

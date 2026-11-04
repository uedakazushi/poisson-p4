# A unit row and the remaining heavy block on Theta=0

Date: 2026-11-04. Status: exact characteristic-zero primary certificate;
independent review pending. This follows the p-branch profiles and uses the
same FULLU ring `R=QQ[a,e,d,1/(adFBCJ)]` and raw checkpoint row labels.

Let `g1=(3e-2)d+20`, `g2=(6e+1)d+40`, and
`Theta=B^2-C*g1*g2`. The p-branch certificate has already shown that g1=0
and Theta=0 are disjoint closed pieces of p=0 on FULLU by the explicit B-unit
comaximality identity. This note concerns only the Theta piece.

The first heavy coefficient of raw row 8 has the exact remainder

```
f = g1*g2^2/37500 modulo Theta.
```

The saved source quotient proves this against the original raw coefficient.
Its inverse on `R/(Theta)` is explicit:

```
f^-1 = 37500*C^2*g1/B^4.
```

Indeed, writing `n=37500*C^2*g1`, the exact polynomial identity is

```
f*n-B^4 = Theta*(Theta-2B^2).
```

For the original, unreduced row entry the certificate adds its saved quotient
times n to the right-hand cofactor. Hence the inverse is valid on the entire
Theta scheme, without assuming irreducibility or reducedness. Only the known
unit B occurs in its denominator. The monic reduction of Theta uses its
literal leading coefficient `2500d^2` in a, also a known unit on FULLU.

All three middle heavy entries `[2,3,4]` of row 8 are literally zero in the
original checkpoint. Use this unit first coefficient to eliminate the first
column from the other rows. The middle three columns remain unchanged.
Multiply the last target column by the known unit B^4. For each remaining
row i its last coefficient is then

```
B^4*M[i,5] - M[i,0]*n*M[8,5] modulo Theta.
```

Rows `[0,2,3,4,5,6]` vanish after specialization: their middle entries were
already zero and both outer entries have zero saved Theta remainder. Deleting
the unit pivot row leaves exactly the seven raw rows
`[7,9,10,11,12,13,14]` and four target columns. The first three are precisely
the existing B7-by-3 block modulo Theta. The last is the displayed Schur
column with its known-unit target scaling.

All 28 source expressions are reduced using the fixed monic quadratic
`Theta/(2500d^2)`. Each is saved with its full quotient/remainder identity
against the original Theta, not merely evaluated at a point.

| Column | Nonzero entries | Total remainder terms |
| --- | ---: | ---: |
| middle 2 | 7 | 290 |
| middle 3 | 7 | 334 |
| middle 4 | 7 | 290 |
| cleared last Schur column | 7 | 2271 |

The largest remainder has 357 terms and a-degree at most one. The last column
is not identically zero; no dependence or independence relative to B has
been inferred. At every residue field on this branch the exact block
criterion is

```
rank(heavy14x5) = 1 + rank(saved7x4).
```

Thus rank five here requires rank four of the actual saved 7-by-4 block.
The next source-bound test can use the accepted cleared row-coordinate
relations of B to examine this extra column on a recorded minor chart.
That would be a new test; it has not been run here. No unknown minor has been
inverted, and no rank/minor/kernel solver, GB, new point or prime was used.
No full cubic- or quartic-rank conclusion follows from this profile alone.
All pure tails remain in the original source representation.

Producer: `code/s61_nonsquare_raw_H_Theta_unit_profile_20261104_v63.py`.
Evidence: `reports/S61_NONSQUARE_RAW_H_THETA_UNIT_PROFILE_20261104_v63`.
The exact run completed in 0.561 seconds under 57-second/2-GB/32-MiB guards.
The report pins the source checkpoint, prior profile and accepted bridge.
Certificate SHA-256:
`c31e3a1e665abb3be428e21f4acf412c800b86f01fa8cc6abb24780e1b302e02`.

# Heavy-column profiles on the two components of p=0

Date: 2026-11-04. Status: exact characteristic-zero primary certificate;
independent review pending. This is a fixed coefficient-profile calculation,
not a new rank or minor computation.

Let

```
C = 100a + 18e^2 - 9e - 2,
B = (50a + 18e^2 - 9e - 2)d + 60e - 40,
F = 500a + 90e^2 + 3e - 24,
J = 5ad - 6e - 4,
R = QQ[a,e,d,1/(adFBCJ)].
```

All source coefficients are those of the frozen raw H checkpoint
`S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63/checkpoint.json`.
The accepted `CUBIC_GENUINE14_H5_PREFLIGHT` and its independent review identify
the remaining fourteen rows and five columns, up to recorded known-unit row
and column scalings, with the actual final cubic presentation's heavy block.
Here a row label always means the raw checkpoint row, not the final row label.

## Exact scheme decomposition

Put

```
g1 = (3e-2)d + 20,
g2 = (6e+1)d + 40,
Theta = B^2 - C*g1*g2.
```

The selected raw checkpoint coefficient is
`p = matrix[2][0] = matrix[4][0]`. Direct full QQ coefficient comparison gives

```
750000*p = g1*Theta,
1 = Theta/B^2 + (C*g2/B^2)*g1.
```

Thus `(g1)+(Theta)=R`, and the Chinese remainder theorem gives the exact
scheme decomposition

```
R/(p) = R/(g1) x R/(Theta).
```

In particular no intersection of these two pieces has been discarded:
the explicit comaximality certificate uses the already known unit B.
Neither g1 nor Theta is being declared a unit on R.

## Saved full quotient identities

The fourteen raw rows are
`[0,2,3,4,5,6,7,8,9,10,11,12,13,14]`. For each row, the first and last heavy
columns `[0,5]` are reduced in both branch models. Each of all 56 entries is
saved with its literal identity

```
source_entry = quotient * original_branch_polynomial + remainder.
```

On g1=0 the monic divisor is `g1/(3d)`, monic of degree one in e over
`QQ[a,d,d^-1]`. It imposes `e=2/3-20/(3d)` and uses only the known unit d.
On Theta=0 the monic divisor is `Theta/(2500d^2)`, monic of degree two in a
over `QQ[e,d,d^-1]`. Its asserted leading coefficient `2500d^2` is checked
literally. These are quotient models before the remaining FULLU localizations;
all displayed identities persist after those localizations.

| Branch | Nonzero rows in either profiled column | Total remainder terms | Largest entry |
| --- | --- | ---: | ---: |
| g1 | 5,7,8,9,10,11,12,13,14 | 302 | 26 |
| Theta | 7,8,9,10,11,12,13,14 | 1079 | 117 |

The two old pair rows 2 and 4 have zero first-column remainder and equal
last-column remainder in both branches; the saved profiles in fact show both
profiled columns vanish in those rows. The rows whose other three heavy
coefficients are literally zero before either specialization are
`[0,2,3,4,5,6,8]`. On g1=0 the surviving such rows are 5 and 8. On Theta=0
only row 8 survives in the two profiled columns.

## Scope and next test

The other three heavy columns `[2,3,4]` retain their accepted seven-row B-block
provenance; they were not recomputed or reduced in this job. No alternate pair
determinant, rank solver, 14-by-5 minor family, GB, new prime or new point was
used. These profiles do not prove heavy rank five on either branch. A B-block
rank drop alone does not imply a rank drop of the full cubic presentation;
the pure components can compensate.

The next small test is to check the two surviving sparse rows 5 and 8 on the
g1 branch against short expressions in `a*d^2-2d+8`, and then check their one
fixed 2-by-2 determinant. That test has not been executed in this certificate.
The Theta branch requires an additional source-bound argument.

## Reproduction and evidence

Producer: `code/s61_nonsquare_raw_H_p_branches_20261104_v63.py`.
Run with `--output` pointing to a new directory; the frozen inputs need no
mutation. The completed primary took 0.124 seconds with a 57-second alarm,
2,000,000,000-byte address-space bound and 32-MiB output guard.

Evidence folder: `reports/S61_NONSQUARE_RAW_H_P_BRANCH_PROFILES_20261104_v63`.
`report.json` pins the producer and all source/bridge hashes.
The certificate SHA-256 is
`6d3bfed276a6d061c3f114f685e2756f269d614cce9125f8e09f9c912a5fa633`.
All arithmetic is QQ Laurent polynomial arithmetic; no finite-field-to-QQ
promotion is used.

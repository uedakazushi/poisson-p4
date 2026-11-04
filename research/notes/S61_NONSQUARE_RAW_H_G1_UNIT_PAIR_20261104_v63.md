# A known-unit heavy pair on the g1 branch

Date: 2026-11-04. Status: exact QQ primary certificate; independent review
pending. The calculation is one fixed pair check, not a minor search.

Use the FULLU ring and raw checkpoint conventions of
`S61_NONSQUARE_RAW_H_P_BRANCH_PROFILES_20261104_v63`. Write

```
g1=(3e-2)d+20,   T=ad^2-2d+8.
```

On this branch `e=2/3-20/(3d)`. The first and last heavy coefficients of the
two raw checkpoint rows 5 and 8 are exactly

```
M = [ T^3/(100d^3)   T^4/(10d^4) ]
    [ T^2/300        T^3/(6d)    ].
```

All four formulas are checked against the saved monic remainders and the
original coefficients. For every coefficient the certificate saves and checks
the full identity `original = g1*quotient + displayed_entry`. All six entries
of these two rows in the other three heavy columns `[2,3,4]` are literally zero
in the original checkpoint, before restriction to either p=0 branch.

There is a literal source-ring identity

```
B-50T/d = g1*(6e+1-20/d).
```

Thus T is a unit on FULLU/(g1). Direct coefficient checks give

```
det(M) = T^6/(750d^4) = B^6*d^2/(750*50^6) modulo g1.
```

The difference between the two determinant representatives is saved with its
explicit quotient by g1, using the six-term difference-of-powers identity.
The denominator in the following inverse is therefore a known unit, using
only the already inverted B and d:

```
M^-1 = [ 125d^3/T^3    -75/T^2      ]
       [ -5d^4/(2T^4) 15d/(2T^3)  ].
```

The saved inverse uses the equivalent common-denominator adjugate form. Both
products, in both orders and all four entries each, are verified by exact
Laurent polynomial comparison after clearing the determinant.

Because these two rows have zero middle three coefficients, they form an
invertible block in the first/last heavy columns. Eliminate those two columns
from the other rows using M^-1. This leaves their middle three coefficients
unchanged. The middle block is exactly the accepted seven-row B block, so at
every residue field of FULLU/(g1),

```
rank(heavy14x5) = 2 + rank(B7x3).
```

Consequently any independently certified finite-or-empty cover of the B-block
rank-defect locus restricts the possible heavy-rank defect on this branch.
This statement does not assert that the B-block defect is nonempty, or that
it is empty. It makes no claim about the Theta branch. No implication from
heavy-rank defect to full cubic-rank defect is used; pure components can
compensate. Pure tails are retained in the original source data and have not
been discarded or expanded into a new full cubic presentation here.

Producer: `code/s61_nonsquare_raw_H_g1_unit_pair_20261104_v63.py`.
Evidence: `reports/S61_NONSQUARE_RAW_H_G1_UNIT_PAIR_20261104_v63`.
The exact run took 0.041 seconds, with 57-second/2,000,000,000-byte/32-MiB
guards. No solver, new point, prime, GB, alternate minor family or Theta
operation was used. The report pins all source and bridge hashes.
Certificate SHA-256:
`dae597437fd49d852e80f585921812fe6e37342d40a47b52e77abdd7bba296a9`.

# The recorded B-minor chart survives on Theta

Date: 2026-11-04. Status: exact QQ primary certificate; independent review
pending. This packet substitutes only one previously saved B minor and its
three specified row-coordinate groups. It does not evaluate a new minor or
any extra-column defect.

The frozen rational Theta isomorphism is used with target coordinates
`(zeta,d)`:

```
a=(zeta+20)*((zeta+20)(zeta-10)+50d)/(500d^2),
e=2/3+(zeta^2-400)/(60d).
```

The input is the exact 407-term Delta and three row-coordinate groups of
`S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63/certificate.json`. Its selected
B rows are block indices `[0,5,2]`, equivalently raw source rows `[7,13,10]`.
The three groups concern block indices `[1,4,6]`, equivalently raw rows
`[9,12,14]`. Every group has three coefficients lambda in that selected-row
order. The matrix convention remains `Delta*row_i=sum_j lambda_ij*row_pivot_j`.

The native FLINT calculation operates in exact bivariate rational Laurent
polynomials. For each source polynomial it clears only a known power of d;
all rational constants are retained exactly. It reconstructs the two
coordinate numerators directly and checks them against the pinned chart.
The complete coefficient inputs and outputs are saved.

The image of Delta is **nonzero**, with 56 terms, zeta-degree 27 and d
exponents between -7 and 0. Thus the previously recorded B-rank-three chart
does meet the Theta surface generically. This nonzero polynomial is still
an unknown, non-global-unit chart factor; it is not inverted on all FULLU.

The transformed lambda groups have term counts

| Raw row | lambda 0 | lambda 1 | lambda 2 |
| --- | ---: | ---: | ---: |
| 9 | 86 | 62 | 69 |
| 12 | 70 | 0 | 56 |
| 14 | 123 | 69 | 78 |

All ten stored transformed polynomials together have 669 terms. The largest
native intermediate has 357 terms. There was no determinant, factor, GCD,
rank, kernel or GB computation, and no point or prime specialization.

The packet also computes coefficient-free support bounds for a *future*
substitution of the already saved seven Schur entries S. It does not
substitute their coefficients. From those support sets the three proposed
cleared expressions

```
Delta*S_i - sum_j lambda_ij*S_pivot_j
```

have the following conservative costs:

| Raw row | Uncollected product-pair bound | Result support bound |
| --- | ---: | ---: |
| 12 | 60,480 | 758 |
| 9 | 91,918 | 789 |
| 14 | 108,201 | 834 |

The raw12 expression is the cheapest next fixed test. None of these three
expressions was evaluated here. Even a nonzero result would first establish
an extra direction on the recorded minor chart; its zero locus and the
Delta=0 part would still require treatment. No uniform heavy-rank or full
cubic/quartic-rank conclusion is made.

Producers:
`code/s61_nonsquare_Theta_Delta_lambda_profiles_20261104_v63.py` and `.cpp`.
Evidence: `reports/S61_NONSQUARE_THETA_DELTA_LAMBDA_PROFILES_20261104_v63`.
The fixed job took 1.135 seconds including compilation, within the
57-second/2-GB/32-MiB and 150,000 stored-term guards. Native input/output and
all source hashes are pinned in `report.json`. Certificate SHA-256:
`0f9f4ef0f844d15c883748d432b8f667538737443c4e094e8aac588cc1cf3900`.

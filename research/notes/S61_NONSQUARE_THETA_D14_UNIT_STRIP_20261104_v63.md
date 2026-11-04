# Exact unit factors of D14 and Delta on Theta

Date: 2026-11-04. Exact QQ primary certificate; independent review pending.
Only gcd and exact division against already certified FULLU unit localizers
are used. No general factorization of an unknown polynomial is performed.

The inputs are the frozen nonzero 269-term D14 and 56-term Delta Laurent
polynomials in `(zeta,d)`. The permitted localizers are d, zeta, zeta+20,
L, K, Fnum and Jnum, together with nonzero QQ constants. All are units in
the actual rational Theta ring by the coordinate-isomorphism certificate.

After extracting the rational coefficient content and known d monomial,
native exact FLINT gcd/division is applied to each of the six polynomial
localizers. At every step both identities are checked and saved:

```
current = gcd * source_quotient,
localizer = gcd * localizer_quotient.
```

The complete product reconstruction is checked at the end. This also
certifies a partial divisor of a named localizer as a unit, without declaring
any unknown factor invertible.

For Delta the twenty-one steps remove eight zeta factors, seven zeta+20
factors and six L factors. The remaining polynomial is the nonzero rational
constant -1. Its d unit exponent is -7. Hence **Delta is a unit everywhere
on FULLU/(Theta)**. In particular the selected three-by-three B minor is
invertible on this entire branch, so B has rank exactly three at every
residue field there. This is an additional exact result; Delta was not
assumed to be a unit in the previous profile or defect packets.

For D14 the forty-two steps remove fifteen zeta factors, fourteen zeta+20
factors and thirteen L factors. Its d unit exponent is -16. The remaining
unknown polynomial is

```
Q6 = -9zeta^4-360zeta^3-3600zeta^2
     -540d*zeta^2-10800d*zeta-8000d^2.
```

This packet retains Q6 as an unknown nonconstant polynomial. The rational
contents are saved exactly in `native.json`; they are nonzero ground-field
units. No inverse of Q6 is introduced, and this calculation does not prove
D14 to be a global unit. The visible possible two-factor expression for Q6
requires its own subsequent fixed identity check.

## Binding D14 to an actual minor

The frozen Delta source determinant is replayed directly with its original
row order `[7,13,10]` and B-column order `[0,2,1]`, equivalently raw heavy
columns `[2,4,3]`. Its full QQ coefficient polynomial agrees with Delta.
The already verified original identities `Delta*B14=lambda*B_pivots` then
give the universal bordered determinant identity

```
D14 = det(rows [7,13,10,14], columns [B0,B2,B1,lastSchur]).
```

This follows by the adjugate identity over the original coefficient field
and therefore is a polynomial identity; it does not require a global Delta
inverse or an additional choice of minor. On Theta, undoing the row8 unit
pivot and the B^4 last-target scaling gives

```
det(raw heavy rows [8,7,13,10,14], columns [0,2,4,3,5])
 = (original_row8_first/B^4)*D14.
```

The prefactor is a certified known unit on this branch. Thus D14 is an
actual maximal-heavy-minor certificate with known unit scalings, not an
unbound source relation. This makes no assertion about the pure
components or full cubic57/quartic239.

Producers: `code/s61_nonsquare_Theta_D14_unit_strip_20261104_v63.py` and
`.cpp`. Evidence: `reports/S61_NONSQUARE_THETA_D14_UNIT_STRIP_20261104_v63`.
The exact run took 1.032 seconds including compilation; 7596 output terms
and peak polynomial size 269 satisfy the 57-second/2-GB/32-MiB/150,000-term
guards. No other minor search, GB, point or prime was used.
Certificate SHA-256:
`5cb20953165c5570b8aeec461be041c3792e765600e9804e395bc3cf1a563844`.

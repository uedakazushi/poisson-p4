# A whole-P100149 lower bound from a fixed minor cover

Date: 2026-11-04. Exact characteristic-zero primary certificate;
independent review pending. The ring is the entire previously fixed chart

```
T100=QQ[z,1/(z(z+20)(z+60)(203z+2700)(z+12))].
```

The actual weighted quartic map has **pointwise rank at least149 over every
geometric point of this chart**. This does not assert rank239, a single fixed
149-minor, another chart or physical excluded points.

## Fixed sources, targets and signs

Use only the already saved15 pure sources `[q5H-cleared;G18*7;K2*7]`.
K2 is the corrected full source `C23-sum_p alpha2,p*C_p`. Keep the eight
old selected local98 target columns `[0,1,2,3,4,5,6,9]`, original122 labels
`[0,1,2,4,5,6,7,10]`. The old8-by8 matrix S has source row order
`[0,2,3,5,7,6,1,4]`. Its determinant D and adjugate H were already proved
and independently reviewed. This job checks both old adjugate products again.

Clear the known whole-curve LCD of each of the seven new K2 rows. All56
selected-entry source identities are checked exactly. For a new row v,
the matrix obtained by replacing selected source row position i has
literal determinant

```
D_(v,i)=sum_(k=0..7) v_clear[k]*H[k,i].
```

All other selected rows and their order are retained, so this formula has
no extra permutation sign. The source LCDs are units. Thus these polynomial
determinants differ from the actual projected source minors only by known
whole-chart units. Every computed dot is checked in both summation orders.

## Exact gcd and polynomial certificate

First strip D using only the five displayed linear unit factors and a
nonzero QQ leading coefficient. Its unit multiplicities are
`[21,93,6,0,1]`; the remaining unknown monic polynomial has degree79.
Every factor division and full reconstruction is checked. Unknown factors
are retained throughout.

Compute the fixed replacement candidates in order: physical light v0,
replaced positions0..7, then physical v1, positions0..6. After15 candidates
the exact gcd is1, so the authorized computation stops. No other new minor
or source is tested. The gcd degrees progress

```
79,79,79,79,79,79,60,60,10,1,1,1,1,1,0.
```

Each replacement's known-unit multiplicities, scalar and primitive part
are saved. The literal extended-gcd identity among primitive polynomials
is verified, and every primitive polynomial divides by the final gcd.
Then all unit denominators are cleared, giving an **actual polynomial
identity among actual projected source minors**:

```
sum_i A_i(z)*D_i(z)
 = E(z)
 = z^21*(z+20)^93*(z+60)^6*(z+12).
```

E is a unit in T100. Only five minors have nonzero coefficients: old D,
and replacements `(physical light,selected row position)`
`(v0,6),(v1,0),(v1,1),(v1,6)`.
Their saved polynomial cofactor degrees are241,248,224,104,8 respectively.
All rational coefficient divisions are by nonzero QQ constants; no unknown
polynomial is inverted. The cleared identity is checked coefficientwise.

Consequently at every geometric point at least one of these eight-by8
minors is nonzero. This gives eight independent actual pure directions,
together with the accepted117 quadratic,23 H1 and one H2 directions.
Therefore the original map has rank at least117+23+1+8=149 throughout
P100. A possible failure of a particular old minor is covered by the other
four; no physical specialization or field irreducibility assumption enters.

## Evidence and resources

Producers: `code/s61_nonsquare_P100_pure15_cover_20261104_v63.py` and `.cpp`.
Evidence: `reports/S61_NONSQUARE_P100_PURE15_COVER_20261104_v63`:

- `source-binding.json` fixes old/new source coefficients, source LCDs,
  row order, local target columns and original122 labels;
- `native.json` saves every computed minor, full unit factorization,
  primitive extended-gcd coefficients and cleared actual-minor coefficients;
- `native-input.txt` pins the old square/adjugate and selected new RF entries;
- `certificate.json` and `report.json` state the resulting chart scope.

The single job took6.564 seconds including compilation; native arithmetic
took1.359 seconds. The actual peak was62684 polynomial terms including a
conservative44870-term literal reserve that counts the entire older18193
input, earlier K2 serialization and old witness. It serialized5745 terms;
the full packet is20203941 bytes, below32MiB. A27MiB native JSON cap is
enforced. No cap was raised or arithmetic job repeated.

Report SHA-256:
`3a7e735027b75f7e9c1ac6d7c77615710da0db387f2129c0f4ffe77c8e8cf6d2`.
Native output SHA-256:
`b2811cdb0a0651339313ac21b98e88bd7a64b1170455a986c7db724de18d126c`.
All primary bytes are immutable. No new NF coefficient, Gaussian/rank
search, unknown-factor inversion, two-row replacement, further source
block, Laurent series, point, prime or GB is included.

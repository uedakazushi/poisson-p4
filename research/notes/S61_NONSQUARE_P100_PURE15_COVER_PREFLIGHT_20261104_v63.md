# A fixed-minor cover preflight for the saved15 pure source rows

This source-only packet performs no polynomial dot product, minor,
known-factor strip, gcd or rank computation. It uses the accepted nonzero
eight-by8 selected minor S from `[q5cleared;G18*7]` and the seven newly
saved K2 products. A large generic15-by98 elimination would repeat a weaker
question than the existing geometric generic theorem. The following fixed
minor family instead targets a whole-curve lower bound.

Keep the old eight target columns `[0,1,2,3,4,5,6,9]` and its saved adjugate
H=adj(S). For each new row v, clear its known whole-P100 unit LCD. Replacing
selected source row i by v gives the literal determinant

```
D_replace(v,i)=sum_(k=0..7) v_clear[k]*H[k,i].
```

This provides56 fixed candidates with no arbitrary minor search or unknown
denominator. If the exact QQ[z] gcd of old D and these candidates is a unit
after stripping only the five known curve factors, this finite family
covers every point of P100. At each point one selected8-minor is nonzero,
so the actual quartic map has pointwise rank at least149 throughout P100.
The original whole142 bound and open149 proof would then have a stronger
uniform counterpart. No such gcd result is asserted here.

A remaining unknown gcd defines only an unresolved finite cover. At a point
where S has corank at least2, every one-row replacement can vanish while
new sources still restore rank through two or more replacements. Thus a
nonunit gcd would not prove actual quartic rank failure.

Only source denominator profiles, polynomial degrees, integer coefficient
heights and support bounds are computed. The56 output polynomials have a
total8743-term upper bound, with830024 naive monomial-product pairs and a
conservative3766797-byte uncompressed JSON bound. The latter uses common
QQ coefficient denominators and integer numerator-height bounds, including
the exact known-unit LCD quotient height. No polynomial multiplication is
used to produce a candidate.

Before any future dot/gcd job, stripping known units from old D alone may
reduce its200-degree univariate input. Any extended-gcd witness still needs
its own output guard. There is no automatic56-dot calculation, higher-row
replacement, full15-by98 rank test, point, prime or GB.

Code: `code/s61_nonsquare_P100_pure15_cover_preflight_20261104_v63.py`.
Reports: `S61_NONSQUARE_P100_PURE15_COVER_PREFLIGHT_20261104_v63`, whose
certificate gives all56 source/degree/height bounds separately. The
source-only pass took0.076 seconds; all primary bytes are immutable.

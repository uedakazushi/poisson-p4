# A small fixed candidate for the next whole-P100 unit cover

This is source/cost inspection only. No9-minor, gcd, new residual coefficient
or rank has been computed.

Keep the five accepted weighted8-source selectors of the whole149 cover.
Append the actual source K2*v3 (row11 of the saved15-source matrix) to each,
and append one target column to the fixed eight selected target columns.
Comparing only cached-entry support and sizes selects **old local98 column25**
as the cheapest common candidate with a perfect-matching support in all five
9-by9 matrices. No polynomial value or rank guides that choice.

The fixed9-source orders are the saved five selector orders followed by
row11. The extra source has its whole-row known-unit LCD, as do all the
selected old rows. Columns are C followed by25, with no permutation sign.
All five determinants are therefore actual9-minors multiplied only by known
whole-P100 row units.

Their determinant-degree bounds are `[231,226,214,227,233]`, giving1136
polynomial terms in total. Integer numerator/common-denominator height
bounds give891374 bytes for the five raw JSON determinant polynomials.
The source entries before LCD clearing total6164 RF terms when counted with
repeated uses; the bound for all cleared entries is8206 polynomial terms.
No coefficient expansion is performed. A future job can
stream one9-by9 matrix at a time and save exact determinant/source identities.
Any gcd/Bezout part still needs its own resource/output guard.

If the exact gcd of these five fixed9-minors becomes a known curve unit,
the original quartic map has rank at least150 on the whole P100 chart.
A nonunit result would mark only an unresolved finite chosen-minor cover;
no actual rank failure would follow, and no new column/source is automatic.

These determinants also avoid evaluating the large global-section
coefficients. If N_i is the defined weighted9-minor and L_w is the new
source LCD, then

```
E*L_w*R[11,17] = sum_i a_i*N_i.
```

Here R is the typed15-by90 residual, whose columns follow the complement T.
Its column17 corresponds to **old local98 label25**. The earlier certificate
string `R_w,25` uses that old target label, not a90-column position; the
supplemental `column-index-clarification.json` records this exact map.

Code: `code/s61_nonsquare_P100_next9_cover_preflight_20261104_v63.py`.
Reports: `S61_NONSQUARE_P100_NEXT9_COVER_PREFLIGHT_20261104_v63`.
The source-only pass took0.187 seconds. All primary bytes and the separate
index clarification are immutable. No determinant/gcd job is launched.

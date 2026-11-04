# Source and cost preflight for70 pure quartic products on P100

Date: 2026-11-04. No new coefficient value, quartic normal form, rank, point,
prime or oracle query is evaluated. The inputs are the ten independently
reviewed actual pure cubic source kernels and the accepted light/NF frames.

## Exact source and proposed coefficient computation

Each kernel K_i is saved with its original16 and actual delta7 source
cofactors. Multiplying it by each physical light v_j gives an actual
combination of the112 original cubic-product sources, using the original16
cofactors at rows `174+7*s+j`. Use the actual delta7 cofactors only for the
preconditioned112-row H1 cache. Its full H1 and H2 parts are zero. Transport the whole
pure cubic and the physical light through the accepted inverse light frame,
form joint210, and apply the fixed NF117 pure98 operator. The normal form
differs from that actual source only by certified quadratic source rows.
No polynomial image independence or70-row rank is assumed.

One future job should handle **one cubic and its seven rows**, writing and
releasing each output before the next. It can reuse the fixed244 scalar and
547 linear NF nodes, while retaining just the one joint cubic84. It need not
query a source DAG or expand core parameter polynomials.

The source-only forecast uses exact numerator degrees and known five-factor
denominator profiles of the saved curve values. Sym^3 of the light inverse
and NF operations are propagated as bounds; no RF coefficients are
multiplied. The seven-row retained-output bounds are:

| Kernel index | Actual source label | Seven output term bound | Max live bound |
|---|---:|---:|---:|
| 0 |18|43746|51755|
| 2 |23|45384|53925|
| 3 |39|43746|51755|
| 4 |40|44880|53669|
| 5 |43|45930|54505|
| 7 |68|42654|50536|
| 10 |92|47568|56407|
| 12 |115|46476|55167|
| 13 |120|47022|55745|
| 15 |134|47568|56407|

Live bounds include literal inputs, retained84 source and joint coefficients,
the inverse frame, scalar NF values, the current210 input, linear-node peak,
an extra98 output copy and10000 temporary terms. Every one-row live bound
passes150000. All70 retained outputs would total454974 terms, so retaining
all70 in one arithmetic job is outside this plan. Output bytes must still
be guarded in an eventual coefficient job; these term estimates alone are
not a32MiB byte guarantee.

## Sparse operator option

The symbolic NF98-by210 coefficient support has925 entries. The symbolic
Sym^3(B)84-by84 support has2698 entries. Composing NF with multiplication by
one physical light gives a98-by84 operator on joint cubic coefficients.
Its support bounds for lights `[0,1,2,3,4,6,7]` are
`[377,1162,1162,149,284,1162,1162]`; numeric RF-term bounds are
`[7555,31588,28346,1500,5109,31588,28346]`.

These are support/degree bounds, not expanded numeric operators. The total
numeric bound134032 leaves insufficient room for all other live input and
output caches under150000. Thus the safe default is the shared typed NF DAG,
or at most one separately bounded sparse operator at a time. The W=v3
operator is particularly small, but no separate W-row job is launched.

The first previously unevaluated candidate is K_2 (actual23) times all seven
lights, with45384 output and53925 peak bounds. This is a proposed next
coefficient input, not an authorized or executed follow-up in this packet.

## A future formal Laurent rank witness

For a later large generic-rank test, formal expansion at z=0 can avoid a
large QQ(z) adjugate. This is an excluded pole of the rational chart, not a
new physical source point. The old prime32003 may be used only after all
saved inputs and every needed denominator leading coefficient are proved
integral and unit there. Saved kernel QQ coefficients pass the integrality
check; the remaining complete-circuit leading-denominator gate is explicitly
pending.

Precision9 would permit a positive fixed-minor/LU witness only when the
total precision loss, or certified determinant valuation, is at most8.
This is a total bound, not8 per pivot. The selected minor must be bound to
the actual rational rows and its nonzero leading coefficient must be saved.
Such a positive residue certificate would imply a QQ(z) rank lower bound.
Zero coefficients or insufficient precision would supply no rank upper or
geometric negative. No Laurent series or rank computation is done here.

## The21 genuine quartic sources still needed later

Their complete H1-cache indices are113..133. Labels are

```
psi(77,131),17,22,38,48,55,58,63,67,78,91,96,98,101,
109,112,118,119,121,126,128.
```

Use the accepted Qhat frame, not raw T34 coefficients:
`Qhat_l=(50/B)^2*F_l(T34*diag(Tden)*v)`. The pair is
`chi=u*Qhat77+v*Qhat131`, with H2 coefficient U=(50aBC)^2,
and `psi=-b131*Qhat77+b77*Qhat131`. Every other retained source is first
`Qhat_l-(b_l/U)*chi`. To obtain a pure source it then needs its actual112-
cubic H1 lift subtracted, with all cubic pure tails retained. Those21 pure
quartic images remain lazy. Their source operators are pinned by the
accepted H_SQUARE_H1_PREFLIGHT certificate; no new quartic root queries
are prepared or evaluated in this preflight. The actual cache ordering
remains `[q5H row0]+112 cubics+21 quartics`.

Code: `code/s61_nonsquare_P100_pure70_source_preflight_20261104_v63.py`.
Reports: `S61_NONSQUARE_P100_PURE70_SOURCE_PREFLIGHT_20261104_v63`.
The source-only pass took0.725 seconds. Report SHA-256:
`2128875be6afa8be37bd65da8e29a23025638c6fccafefe350a79a8a231dc3e8`.
Certificate SHA-256:
`762039b9258ce9d4dd19b6529530cd19f01a2381bd5d68a6793b28dbafc9e40e`.
All primary bytes are immutable. No70-row rank or full239 claim is made.

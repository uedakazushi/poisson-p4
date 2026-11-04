# One actual pure98 source row on P100

Date: 2026-11-04. Exact characteristic-zero primary certificate; independent
review pending. Only one complete pure row is evaluated. No matrix-rank or
subsequent coefficient job is performed.

## Actual source and normalization

The source is the **actual producer q5H cache row0**, not a hypothetical
unscaled Hq5. All23 cached H1 coefficients are compared literally over
QQ[a,e,d] with row57 of the saved QPAIR matrix at its remaining target
columns. Its98 pure coefficients and H2 coefficient are literally zero.
Thus the separate H2-unit source deletion leaves this row unchanged.

Multiply this actual H1 row by the independently verified whole-P100
112-by-23 right inverse. Only42 cubic-product sources occur. The saved
delta7 source change is undone before taking original pure tails:

```
k_old6=delta*k_new6,
k_old1=k_new1-delta*k_new6.
```

The other coefficients are unchanged. All42 cofactor sums are checked in
both summation orders and the inverse source transport is checked.

Only the six original genuine cubics with indices
`[1,6,8,9,11,14]`, original135 labels `[20,61,85,86,99,129]`, are needed.
Their504 pure cubic coefficient roots were prepared in the frozen
source-only preflight. Their values are now evaluated on the full P100
rational chart. The exact frame adapter uses

```
Pi_i = F_actual(T34*D*v) / frozen_basis57_row_unit.
```

The row-unit operation is **division**. The normalization, monomial order
and source-root mapping are pinned by the preflight input hash. No full
source graph or old coefficient query is regenerated.

The actual pure quartic source is

```
Psi = -sum_(i,j) k_old(i,j)*Pi_i(v)*v_j.
```

The minus sign comes from subtracting the full cubic lift from q5H, whose
pure part is zero. This keeps all actual cubic pure tails. Cubics have zero
H2 part, and their H1 projection is exactly the cached right inverse, so
the full source combination has zero H1 and H2 parts. No H-linear equation
is set to zero as an ideal shortcut.

## Coordinate transport and normal form

The six pure cubics are transported through the accepted light inverse to
joint coordinates `(L,ell,x0,x3,r,Z,W)`. Their42 light cofactors use exactly
the same inverse. All210 joint pure quartic coefficients are assembled and
verified a second time coefficient by coefficient using the target-centric
sum over divisors of each monomial.

Only the accepted pure98 normal-form closure is then evaluated:244 scalar
and547 linear nodes. Its input is the actual joint210 pure vector, with
all other weighted target components zero. No quartic source, all134 H1
row family or POWER circuit is queried.

The result has **98 nonzero components**. None of its individual components
is certified a unit after testing the known P100 linear localizers. Their
total numerator and denominator term counts are3098 and1501; maximum
degrees are42 and38. The row is therefore a nonzero pure source on a
nonempty open of P100, giving at least one generic pure direction beyond
the141 unit-block lower bound. This gives a generic rank lower bound142
on this curve, and pointwise on the open where any chosen nonzero component
does not vanish. No uniform pure pivot has yet been exhibited.

Absence of an individual unit component does not prove the component ideal
proper: an exact univariate gcd/Bézout combination could still give a unit.
That test has **not** been executed. There is no rank98 or full239 claim.

## Reusable evidence

- `cubic504.json`: all six original normalized pure tails, with84 coefficients
  each in the preflight's original physical monomial order.
- `cofactor42.json`: actual q5H row0 values and both source frames of its42
  lift coefficients. Light order is `[v0,v1,v2,v3,v4,v6,v7]`.
- `joint-cubic504.json`: all six transformed cubic tails. Its monomial order
  is lexicographically increasing seven-tuples of total degree3, in joint
  variable order `(L,ell,x0,x3,r,Z,W)`.
- `joint210.json`: all210 checked coefficients, with explicit monomial labels.
- `pure98.json`: all98 exact normal-form values and their original122 target
  labels; all denominators are certified whole-curve units.
- `source-binding.json`: the actual row0/QPAIR57/H2-deletion comparison and
  the exact FrameAdapter/root provenance.

These caches can support a later structured pure-kernel construction
without repeating these504 coefficient values. No additional pure source
has been evaluated here.

## Resource and validation record

Producers: `code/s61_nonsquare_P100_q5_pure98_row_20261104_v63.py` and `.cpp`.
Evidence: `reports/S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63`.
The single job completed in8.699 seconds including compilation; native
arithmetic took3.720 seconds. It evaluated19,960 cubic DAG nodes, checked
all42 lift coefficients and all210 combination identities, and produced
all98 normal-form outputs. The optional transformed cubic cache also fit.
It serialized41,312 terms in about2.8MiB.

The live RF polynomial counter includes retained atoms, values, outputs,
temporaries and an8746-term literal-input reserve. Its recorded peak is
77,397 terms. One unit-test gcd temporary can hold at most6 extra terms;
the supplemental budget receipt gives the conservative actual bound77,403,
well below150,000. No limit was raised or job repeated.

There was no new point, prime, GB, Gaussian, all33-quartic or134-row
evaluation, and no automatic follow-up query. The source and output files
are immutable. Report SHA-256:
`5d3ce231f5fab85b2b26543f2fe80a6b37b44850d36f9580dc2960d69be0f68d`.
Pure98 output SHA-256:
`a8873b6c5d88786328f69d7fdbd8bbf4d3424c45d2abffbbb82cc819ac95ffe4`.

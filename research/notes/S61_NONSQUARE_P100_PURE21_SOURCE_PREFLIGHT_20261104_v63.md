# Full-source preflight for21 pure quartic rows on P100

Date: 2026-11-04. This is an input and cost receipt, with no new curve
coefficient values or normal forms evaluated.

## Actual full cubics

Stream the immutable normalized `remaining.txt.gz` in
`S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62`. Its full header is
`57 91 57`; the target-rank57 field is metadata, not an all-point-rank
assertion. The uncompressed full matrix hash agrees with the accepted
producer. Basis positions41,44,48 are normalized original63 rows42,45,50,
and correspond to actual135 equations18,39,68. Their genuine16 indices
are0,3,7.

All seven H-linear coefficients of each full cubic are literally zero,
both in this complete91-column matrix and in the independently reviewed
raw16-by-7 H matrix. These are full H-free cubics, not merely H-free
representatives after restricting Z=W=0. The previous functional-input
receipt also records that all three G11 subtraction coefficients are zero.
Their delta7 source indices are unchanged.

The accepted frame adapter identifies each saved full row with

```
Gi(v)=F_actual(T34*diag(Tden)*v)/row_unit_i.
```

The three units are respectively `189a/F`, `-189a/F`, and `-4779a/(5F)`.
The saved252 coefficients already include this normalization: no additional
unit multiplication or division is applied in this preflight.

The full saved coefficients have231 nonzero entries and4996 QQ[a,e,d]
terms, at most40 per entry, degree8, and24-bit rational coefficients.
Their source input is116190 bytes. It is complete, so a later curve job can
specialize these existing polynomials directly, without any new oracle
query or cubic graph reconstruction.

## Source and ideal transport

Each of the three cubics is multiplied by each of the seven physical
lights `[v0,v1,v2,v3,v4,v6,v7]`. These21 actual319 source rows have indices
174..180,195..201,223..229. The full definitions and normalized63 row
labels are saved individually in `source-input.json`.

They have no H1 or H2 terms. Thus neither the separate H2 unit deletion nor
an H1 lift changes these source polynomials. To compute their actual pure98
images, transport each full cubic and its physical light through the
accepted inverse light frame into `(L,ell,x0,x3,r,Z,W)`, form joint210, and
apply the fixed NF117 operator.

The complete122-output ancestor audit of that operator shows that its
typed linear map preserves H-degree. In particular, a pure joint210 input
has only pure98 outputs. Each normal form differs from its full actual
source by certified quadratic source rows. This is the precise source-span
argument. The old three-by26 ZW-restricted input is **not** used to replace
the full cubics; its discarded terms have not been assumed to lie in the
full quadratic ideal.

## Bounds and next operation

Only exponent and existing-operator bounds are propagated on the whole
P100 ring. No polynomial is expanded after curve substitution. The252
original coefficient values have a conservative4403-term live bound. A
common coefficient bound after the four inverse-light factors is numerator
degree57 and five-localizer denominator exponents `[34,12,4,0,0]`, hence109
terms per joint210 coefficient. The fixed normal form gives14972 terms as
a conservative bound for one complete pure98 row.

Retaining all21 output rows would give314412 terms, beyond a150000-term
retained-output limit. Seven rows from one cubic have the lower bound
104804 for that same conservative output estimate. A future bounded job
must supply a last-use release/streaming plan or select a smaller batch;
the live budget must include the retained inputs, inverse-frame values,
normal-form scalar values, current210 input and current98 output. This
receipt does not authorize or launch that arithmetic.

No rank, independence, unit pivot, matrix completion, full57 or full239
claim is made. The earlier finite-degree functional negatives do not imply
any negative rank conclusion for these21 actual source rows.

## Evidence and frontier

Code: `code/s61_nonsquare_P100_pure21_source_preflight_20261104_v63.py`.
Reports: `S61_NONSQUARE_P100_PURE21_SOURCE_PREFLIGHT_20261104_v63`, containing
the full source input, bounds certificate and pinned-input manifest. The
successful source-only pass took0.358 seconds. No value job was executed.

A first metadata audit visited unrelated earlier variable nodes in the
parent linear DAG and raised `KeyError0`; it evaluated no values. Its source
and frontier are preserved separately in the `_FAILED_FIRST` folder. The
repair restricts that audit to the actual122-output ancestor closure. An
additional empty-directory refusal is also recorded. Total attempts stayed
below2 seconds, without changing a mathematical scope or resource limit.

Report SHA-256:
`564f74f12deec7b76779a9e12e2736c706c071c90b9096b2573898bd58ccb349`.
Full source input SHA-256:
`b9d34f452bd27dcac3117c021c803b157404fedafe841d8f46c4400ad99f8fda`.
Certificate SHA-256:
`1c51039f459affba56bb7b12b1c15852a19f62706e6054956b771332c3d22b68`.
All primary bytes are immutable.

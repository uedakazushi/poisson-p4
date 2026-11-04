# Source-only budget for one P100 pure98 row

Date: 2026-11-04. Frozen input preflight. **No cubic coefficient values,
curve-specialized source values or pure98 values have been evaluated.**

The proposed single full source is the remaining q5H row minus its lift
through the already verified P100 H1 right inverse. q5H has zero pure part,
so the resulting pure quartic is minus the pure tail of that cubic lift.
Both curve right inverses use exactly the same six genuine cubic sources,
with indices `[1,6,8,9,11,14]` and original equation labels
`[20,61,85,86,99,129]`. No other genuine quartic is needed for this first row.

## Exact input and normalization

The frozen frame adapter is used as

```
FrameAdapter.basis57_coefficient(41+i, monomial)
```

for each selected i and all 84 pure light monomials of degree three. The API
extracts `F_actual(D*v)` and **divides by** the frozen source row unit. Its
source starts in T34 and is transported to the actual Tn frame; there is no
extra quartic normalization for these cubic tails.

All 504 scalar coefficient roots have been prepared as arithmetic circuits.
Both `core_value` entry points were replaced by a forbidden-call guard, so
no root was expanded to its parameter rational function. The fixed source
circuit remains immutable; the new queried roots and only their required
node/atom closure are saved in `cubic-coefficient-input.json`.

The actual cache indexing is q5H row0, cubic products rows1 through112, and
quartics rows113 through133. To combine the right-inverse coefficients with
the original six tails, undo the actual delta7 source change:

```
kappa_old6=delta*kappa_new6,
kappa_old1=kappa_new1-delta*kappa_new6.
```

Every other selected coefficient is unchanged. The q5H H1 row has only
526 numerator terms and 23 denominator terms. Its multiplication by the
already saved right inverse is a *future* computation, not performed here.

## Measured source closure and bit sizes

| Item | Count |
| --- | ---: |
| Prepared normalized cubic coefficient roots | 504 |
| Required cubic scalar nodes | 19,960 |
| Required cubic atoms | 604 |
| Literal cubic atom numerator/denominator terms | 3491 |
| Maximum cubic atom coefficient bits | 22 |
| Required pure98 normal-form linear nodes | 547 |
| Required normal-form scalar nodes | 244 |
| Required normal-form atoms | 122 |
| Literal normal-form atom terms | 4404 |
| Maximum normal-form atom coefficient bits | 46 |

The saved q5H coefficients have maximum coefficient size30 bits and the
already computed right inverse99 bits. These latter figures are read from
frozen values, not obtained by a new specialization.

The 98 target labels and the complete scalar/linear ancestor lists are
bound to the accepted NF117 operator. The projection is applied to joint
coordinates, so a future calculation must first use the saved old-physical
to joint light inverse. The actual variables and all row-unit scalings are
retained in the manifest.

## Conservative whole-curve size bounds

Every denominator is expressed using the known P100 linear factors
`[z,z+20,z+60,203z+2700,z+12]`. Bounds use their componentwise least common
multiple when adding circuits; they do not depend on unobserved coefficient
cancellation. The known source denominators are factored only against
`[a,d,F,B,C,J]`, reusing the accepted atomic unit certificates where present.

| Stage | Conservative bound |
| --- | ---: |
| Cubic DAG live rational-polynomial terms | 99,432 |
| All 504 retained coefficient outputs | 24,314 terms |
| Joint pure210 common numerator degree | 109 |
| Joint pure210 common denominator powers | `[61,23,10,2,0]` |
| All 98 final expanded numerator/denominator outputs | 24,576 terms |
| Pure98 linear normal-form live values | 26,157 terms |

The native evaluation plan must release last-used cubic nodes, consume the
504 tails while forming the joint pure210 vector, and release them before
normal-form evaluation. It must retain the hard 150,000-live-term,
2-GB/32-MiB guards. The measured input graph is small enough to select the
**full single98 row** as the next bounded job, rather than expanding all
quartic sources. Runtime remains an estimate until that new job is
authorized and executed. No such evaluation job has started.

## Assessment of a structured 89-row kernel

The known H1 right inverse gives an abstract split kernel of rank89 for
the cubic112-to-H23 projection. On this localized one-variable ring that
projective module is free. This does not make an unconstructed list of89
rows into an audited basis.

A candidate structured construction is:

- Split the 16 cubic sources using the actual six-row H basis, leaving ten
  cubics with zero H-linear part; their seven multiples give70 rows.
- The six independent linear H parts have15 degree-one Koszul relations.
  These cancel only the H1 part of the full cubic sources; their pure tails
  are generally nonzero and must remain as actual pure quartic generators.
- The five quadratic H relations intersect the quadratic-evaluation kernel
  in rank4 because q1 has unit evaluation. Four explicit lifts can use the
  same universal formula as the right inverse.

There is a valid split-exact module explanation of the count70+15+4=89.
To promote those particular vectors to a basis one must still construct
their actual source cofactor matrices, verify the invertible16-source split,
the15 Koszul and four lifted relations, and independence in that source
frame. None of those89 vectors has been constructed or claimed independent
by this preflight. Their evaluation would require all16 cubic pure tails,
or1344 coefficient roots, rather than this first row's504. That larger
query has not been made.

After such a construction, a possible redundant pure98 presentation has
43 pure quadratic rows,89 cubic-kernel rows, one q5H clearing row and21
H2-free/H1-cleared quartics:154 generators. This is a prospective source
recipe, not a computed154 matrix, a new free source quotient or a rank98
certificate. Full quartic239 remains unresolved.

## Reproduction

Producer: `code/s61_nonsquare_P100_pure98_source_preflight_20261104_v63.py`.
Evidence: `reports/S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63`.
The source preflight took4.396 seconds. It created no point, field-value,
GB, matrix-rank or native coefficient job. Input hashes are pinned in
`report.json`; supplemental bit bounds are in
`additional-budget-receipt.json`.
Manifest SHA-256:
`6f60670d78e895c3220e31a8cac42148020ae7d69ebc8638604508f0edda5fec`.
Cubic coefficient-input SHA-256:
`903d7e023eefaaf43f753213d81d592b4e728b96d2c63b7e7e410b3a0171eea4`.

# Deeper binary: core PID exploration (WIP, 2026-11-01)

This archive is a continuation checkpoint, **not a new exhaustion theorem**.
The v49 eight-factor result is unchanged. Unknown factors of the fourth-order
minor remain unresolved over QQ(e,d).

## Algebraic reduction and verified observations

On the seven-parameter deeper chart, set the five nuisance parameters to zero.
The quadratic ideal contains all products of either v5 or v7 with a normal
variable. In degrees at least two, projection to the other seven variables
therefore reduces the fourth-degree coefficient space from495 to210 dimensions.
The projected quadratic rows are the seven primitive QQ[e,d] rows stored here.
Their common contents before primitivization are rational constants.

The quadratic-linear multiplication matrix has generic rank129 on216 rows.
The existing114 constant pivots give86 exact linear syzygies. The additional
syzygy in `extra_syzygy.txt` is verified against the saved quadratic rref over
QQ(e,d); its addition raises the rank of the86 coefficient vectors to87 at
(e,d)=(2,3), modulo32003. Thus all87 generic linear syzygies are accounted for.
The additional syzygy gives **zero** projected quartic at nuisance zero.
This does not say its full nuisance-dependent quartic is zero.

The86 other syzygies give56 nonzero projected quartic rows. With7 quadratics
multiplied by28 quadratic monomials and2 cubic rows multiplied by7 variables,
this gives a266-by210 polynomial matrix. Its discarded numerator contents and
rational denominators have been checked overQQ: all d-dependent factors belong
to the seven already treated factors in `horizontal_factors.txt`. Factors in e
remain coefficient-field units here and must be covered in a final theorem.

At **e=2, over F_32003[d]**, polynomial row reduction gives full rank210 and
cokernel length430. Its determinant support consists exclusively of the known
seven factors; their removed degree totals are0,42,42,1,2,28,315, leaving degree0.
The chosen28 quartics alone leave extra support in this slice. The other
syzygies are therefore useful for removing artifacts of that chosen basis.
A separate QQ reconstruction of the entire matrix and the first55 pivots with
e-only numerators completed. The subsequent QQ(e)[d] computation is unfinished.

These facts **do not prove** the absence of extra factors overQQ(e,d).
Affine components may disappear at infinity on reduction, or merge with the
localized factors. A mod-p computation after inverting D does not permit a
Nakayama argument for the localized module.

## Reconstructing the finite-field slice

Run from the repository root:

```sh
.venv/bin/python research/archive/binary_deeper_pid_20261101/run_probe.py \
  --output work/deeper-pid-archive-replay.json \
  --work-dir /tmp/poisson_deeper_pid_archive_replay
```

Use an absent output and work directory. The script uses the frozen v47
preparation functions, reconstructs all inputs from repository certificates,
compiles native FLINT helpers, and regenerates the one-variable module and its
support calculation. It needs the existing local FLINT headers/library and
severalGB of temporary space. Add `--exact-inputs` to reconstruct the large QQ
matrix and its55 e-only pivots as well. Large caches are intentionally omitted.
`observations.json` records their hashes and the independently executed stage
results. The new aggregate driver has been syntax checked; its entire pipeline
has not been replayed as one command. Each main mathematical stage was run in
scratch before archiving. No independent reviewer is claimed.

The saved seven quadrics come from `primitive_quadratic_rows` in the uniform
quadratic report, by nuisance-zero specialization and removal of v5,v7, followed
by row content division. The seven horizontal factors are exactly its minor
factors that depend on d. `extra_syzygy.json` records the coefficient solution
before common-denominator clearing; `extra_syzygy.txt` is the primitive
polynomial version in the full24-by9 source convention.

## Current continuation

`generic_division_popov.cpp` and `exact_division_popov.cpp` implement experimental
polynomial row reduction overF_32003(e)[d] andQQ(e)[d]. Row contents may be divided
only by e-only units and by the listed known factors. They are unfinished
exploration, not frozen proof producers. Full coefficient-field arithmetic
currently causes substantial coefficient growth.

The valid lifting option is to construct a **new polynomial row matrix** by
QQ-verified divisions by known factors, then prove its cokernel finite over
V=Z_(p)[e]_(p), with residue fieldF_p(e), using an annihilating polynomial whose
highest d coefficient is a p-unit. If its unrelocalized mod-p cokernel is zero,
Nakayama gives zero overV and henceQQ(e). This cannot be applied directly to
M[D^-1]. An exact QQ(e)[d] row certificate is another route.

If every horizontal core curve has a nuisance-zero upper65 witness, openness
and the same actual relative binaryC3 lower65 give a nonempty removed open on
the corresponding16-dimensional stratum. Source generic smoothness and a
complement of dimension at most15 then remove its general-beta fiber. The
nuisance-independence of all86 quartics is unnecessary for that restricted use.
Vertical core factors introduced by e-only denominators still need a finite
cover, or an independently proved second PID argument with the core roles
exchanged. No ambient or global statement is presently inferred from this WIP.

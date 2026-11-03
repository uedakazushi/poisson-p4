# V6: full cubic 27×61 coker reduction (frozen checkpoint)

This packet contains an exact reduction, not a new all-point rank theorem. On
`a F B C d J != 0`, where `J=5ad-6e-4`, the full cubic rank is 30 plus the
rank of the saved 27×61 matrix. Its generic rank is 27 because the independently
audited original cubic rank is 57. Pointwise lower rank, free cokernel, the
quartic length65 cover and ambient nonsquare exhaustion remain open.

The starting 57 rows are the independently audited alternative ad basis from
`S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62`. We delete original rows
48,55,59,57,61,22 and retain the previous 30 pivot rows. The first 21 pivots
use only a,F,B,C. Eight further pivots are constant column combinations and
use only the six named localization factors. The last degree1 witness is
**a column replacement**, not an append: at relative column14 both the old
witness coefficient and target-row entry vanish; adding that existing column
gives new coefficient1, hence an invertible column operation and a true
28×62 → 27×61 reduction. The target original row is23. Every pivot is checked
as a product of declared units, and every subsequent division is exact.

`primitive.cpp` verifies every rational row/column scalar equivalence and
reverses every monomial of the P=ad chart conversion. `column-clean.cpp`
verifies every entry after common known-unit column division by multiplying
back its recorded factor. In its serialized three-variable format the third
variable is named d but represents P. Thus the final localization is
`a F B0 C P J != 0`, with `B0=(50a+18e²−9e−2)P+(60e−40)a` and
`J=5P−6e−4`. The final matrix has 1,310,260 terms; its large regenerated
outputs are omitted with hashes in the manifest.

To resume, extract the archive at repository root. Existing public sources
and reports are dependencies; the archive includes the two pinned old work
inputs needed for byte equality checks. Restore old `fresh-cubic/schur` from
`reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/matrices.tar.gz`.
The replay sequence, with `w=work/nonsquare-v62/cubic27-v6`, is:

```
source /workspace/cloud-setup/activate.sh
python $w/prepare.py
# Compile each named .cpp with c++ -O2 -std=c++17 and pkg-config flint.
$w/forced21-bin $w/raw-alternative $w/after21 $w/original-row-ids $w/forced21
python $w/prepare-span.py
$w/forced-span-bin $w/after21-ordered $w/after29 $w/old-units 1 $w/forced8
python $w/prepare-last.py
$w/last-pivot-bin $w/after29-ordered $w/full-row27 $w/old-units 1 $w/replacement-witness 1
$w/primitive-bin $w/full-row27 $w/integer27 $w/P-integer27 $w/scalar-equivalence $w/degrees.txt
$w/column-clean-bin $w/P-integer27 $w/P-units $w/P-clean27 $w/column-unit-multipliers $w/column-degrees.txt
```

The QQ(e) M21.24.11 GB attempt was capped at240s/3.3GB virtual memory. It
returned timeout124 after input-ready, without any GB, inverse or annihilator
result. No mathematical assertion follows. Its exact data generator and input
hash are saved. A tiny QQ(e) ChangeMatrix/inverse fixture passed. A new thin
Singular wrapper adds its existing MOD directory to SINGULARPATH so the
optimized finite-field procedures load; its direct QQ(e) liftstd fixture also
passed. Singular alternatives after this freeze are mutable new work and are
not evidence for a coker theorem.

If exact coker zero is proved after localization to each QQ(a), QQ(e), QQ(P),
finite generation gives a nonzero univariate annihilator in each coordinate,
so the actual rank-defect support has dimension0. A single coordinate-field
unit result would not suffice. No modular sample can replace these identities
without an additional proper model or accepted homogeneous certificate.

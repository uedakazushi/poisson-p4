# Frozen V24 scalar blocks and coefficient-unit exposure

The definitive accepted results are the immutable public code/notes/reports:

- S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v62 and its independent review:
  original quadratic-product rank at least35 on adFBCJ≠0, not41.
- S61_BINARY_NONSQUARE_SCALAR_PRODUCT_ROWS_20261103_v62 and independent review:
  original full26 rows0,7 admit M T=aI2; preferred scalar step32.
- S61_BINARY_NONSQUARE_SECOND_SCALAR_PIVOT_20261103_v62 and independent review:
  transported138-term cofactor and constant cofactor-span1 give explicit step33.

Current primary generated source is
`work/nonsquare-v62/main-scalar-unit-v21/full-after-two-rows.txt`:24×58,
2,755,382 terms, eight remaining product rows and sixteen genuine cubic rows.
It is an expanded cache omitted from Git; its hash/size are in manifest.json.
All-point product41, cubic57, quartic65 and ambient exhaustion remain open.

Extract the archive at repository root. Reuse the existing activation script.
Restore the pinned full26 from V6 then V9. The preceding public compressed
certificates contain both original source rows, scalar vectors and the2-row
right inverse, so all small mathematical identities can be checked without
rebuilding growing caches. The exact public verifiers currently also pin the
available native source/output caches; rebuilding those is deterministic.

Native step32 uses the V9 one-witness helper:

    witness-stats-bin FULL26 AFTER32 P-units 1 row7.witness 7

It must report replacement43,25×59 and1801170 terms. Then run the included
`transport_second.py`, which reads the saved2-row right inverse and row7
witness, and compile the PUBLIC utility
`code/polynomial_cofactor_constant_span_pivot_20261103_v62.cpp` with FLINT.

    constant-span-pivot AFTER32 AFTER33 P-units 1 after-row7-second.witness 6 second-cofactor-unit-combo.txt

It must report five constants exposing coefficient50, replacement50,24×58,
2755382 terms and PASS1. The optional eighth argument is a column-count plus
rational constants file. Every full polynomial image is calculated before
and after this invertible QQ column change and asserted identical.

The private2×2 fixture and independent review's different2×3 fixture both
validate the new branch. No coefficient is accepted solely from modular rank.
The saved original3-row test uses ALL original coefficients. An earlier
per-entry stripped3-row module test is explicitly INVALID and excluded as
mathematical evidence: entrywise unit stripping is safe for scalar ideals,
not independently for different rows of a matrix.

Separate five-light sufficient20×35 criterion progress is also preserved:
row6 has163-term degree3 scalar lift of−a²P²J; append+pivot gives10×26,
124869 terms (ten unit pivots). This does not prove criterion rank20 or
original product41. Its growing output is omitted with SHA and replay data.

Completed bounded module attempts in STATUS.json produced no global GB or
all-point claim. Later experiments continue separately; this archive and
its source/status/manifest are immutable.

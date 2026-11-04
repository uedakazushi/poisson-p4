# Explicit unit split of the P100 pure target

Date: 2026-11-04. Exact finite target/source transport; independent review
pending. No Gaussian, rank or additional coefficient query is performed.

Let r be the accepted actual q5H-cleared pure98 row, and v one of the seven
saved G18-times-light pure98 rows. The accepted literal Bezout certificate
gives two polynomial coefficients c0,c1 in QQ[z] with

```
r0*c0+r1*c1=1.
B2=[[c0,-r1],[c1,r0]],  B2^-1=[[r0,r1],[-c1,c0]].
```

Its determinant is1. Both inverse products, the determinant and rB2=(1,0)
are checked exactly over QQ(z). Apply B2 to the first two target columns,
then clear each later column j by subtracting r_j times the first column.
This is an explicit invertible98-target transformation. There are no
unknown-factor inversions or projective-complement assumptions.

For each saved source row, define

```
a=v0*c0+v1*c1,
w1=-v0*r1+v1*r0,
w_j=v_j-a*r_j, j>=2.
```

The whole98-entry source reconstruction is checked:

```
v0=a*r0-w1*c1,
v1=a*r1+w1*c0,
v_j=w_j+a*r_j, j>=2.
```

Thus subtracting a times the actual q5 pure source clears the pivot and
leaves a genuine97-component pure source w. All98 q5 identities and all686
G18 identities pass. All saved coefficient denominators remain products of
the five known whole-P100 curve units. The target order is the transformed
old column1 followed by old columns2..97.

The result is an actual unit direction plus seven residual rows. The same
target formula applies to every other source row if later needed; those
other source coefficients remain lazy. The seven-by97 block alone is not
claimed to present the whole quartic cokernel. Its row independence,
generic rank and whole-curve invariant factors are uncomputed. The known
pointwise full-quartic lower bound142 is unchanged.

## Bounds and evidence

The source-only known-unit denominator profiles gave a73662-term live RF
bound and107066-term serialization bound. Each residual row is streamed
directly as gzip JSON and released; no uncompressed giant output matrix is
saved. A28MiB compressed-output guard reserves4MiB for all other files.
The actual peak was49776 live terms including the literal inputs. Native
arithmetic took6.643 seconds and the whole job12.112 seconds. The seven
rows contain40106 serialized terms; compressed rows occupy22612816 bytes
and the full packet23705660 bytes, below32MiB.

Code: `code/s61_nonsquare_P100_pure97_split_20261104_v63.py` and `.cpp`.
Evidence: `reports/S61_NONSQUARE_P100_PURE97_SPLIT_20261104_v63`, with
`preflight.json`, `certificate.json`, `audit.json`, and
`pure97-row0.json.gz` through `pure97-row6.json.gz`. Each row saves its
actual319 source index and its q5-source subtraction coefficient a.

Report SHA-256:
`9f86e680ef96d565cb75969af7f38b04b921e41be7563c389216b90932fce559`.
All primary bytes are immutable. No source/oracle values, point, prime,
minor, rank, Gaussian or GB were generated beyond this authorized finite
transport.

For a future generic-rank test, the original compact eight-by98 block
`[r;G18*v_j, j in (0,1,2,3,4,6,7)]` is preferable to expanding the large Bezout-derived
residues: the exact transport proves its rank is one plus the rank of the
seven-by97 block. This is only a proposed input choice, not an executed
field calculation or a rank conclusion.

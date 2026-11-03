# The34th constructive product pivot as a complete exact Schur circuit

2026-11-03. The full resulting23×58 matrix is saved as a finite exact arithmetic
circuit. Its1334 entries are complete polynomials, not a degree truncation or a
finite jet. The expanded native cache remains the33-pivot24×58 matrix; the saved
constructive circuit performs pivot34. Product41, full cubic57, quartic65 and
ambient exhaustion remain open.

## Inputs and all entries

The source A is the pinned actual24×58 QQ matrix from the independently accepted
33rd-pivot packet. The independently accepted58-entry vector q satisfies
`A_7 q=g=aP²J`, J=5P−6e−4. It has975 terms in33 active components. Source variable
d means P=ad; the old open is aPFB0CJ≠0. These SAME input files and their complete
SHA256 values are leaves of the saved graph.

For every old row i, define the exact polynomial circuit

    v_i = sum_(k=0..57) A_ik q_k.

Append this column. For each remaining old row and EVERY old column save

    B_ij = g A_ij − v_i A_7j.

The ordered old rows are `[0,1,2,3,4,5,6,8,9,...,23]`; all58 old columns remain
in their original order. The first seven resulting rows are products, followed
by sixteen genuine cubic rows. All1334 root formulas, all24 full-image formulas,
all1392 original entry references and all33 nonzero cofactor references are in
the6244-node acyclic graph. Zero cofactor components are explicitly recorded.
The verifier regenerates and compares EVERY node and root from the frozen inputs.

No giant constant cofactor-span transformation is used. Such a23-constant
combination can expose aPJ, but its constants reach about67000bits. The appended
column construction avoids those constants and introduces no new core divisor.
No expanded outer-product cache is attempted.

## Exact cokernel equivalence, including nonreduced bases

This is an identity over any commutative coefficient ring in which g is a unit.
Appending v=Aq leaves the column image unchanged. Keep pivot row7, and replace
each other row R_i by `g R_i−v_i R_7`. The row-operation determinant is g^23,
a unit. The appended column now has only the pivot entry g. For each old column j,
subtract `(A_7j/g)` times that appended column. This clears its pivot-row entry
and leaves all B_ij unchanged. After permutations the matrix is exactly
`diag(g,B)`. Its cokernel is the cokernel of B because g is a unit.

Only g is divided by, and g is already a product of the named units a,P,J.
There is no division by any unknown factor of q, by a cofactor, or by another
polynomial in the core. The argument proves module cokernel equivalence and
works over nonreduced rings; it is stronger than comparing pointwise ranks.

## A different full polynomial fixture

Over QQ[s,t,h], localize ONLY at s. Take

    A = [[1,0,h],[s,s,0],[0,1,0]],   q=(t,1−t,0),   pivot row1, g=s.

Neither nonzero cofactor is individually a named unit. The graph evaluator gives

    B = [[s(1−t),−st,sh],[−s(1−t),st,0]].

The verifier explicitly constructs the invertible row and column operations,
checks their determinants s² and1 and the exact unit block diag(s,B). The source
maximal minor is sh. The resulting maximal-minor ideal is (s²h): its two nonzero
minors are s²h(1−t) and−s²ht, whose difference is s²h. Thus the unknown factor h
survives after localizing at s. All of these are full polynomial identities.

## Finite replay

`circuit.json.gz` is the entire new matrix, with referenced source/cofactor files
and hashes. Restore the33-pivot input by V6/V9/V24 recipes if its expanded cache
has been removed. `build` regenerates the exact graph; `evaluate` is a generic
polynomial/rational backend evaluator and is exercised on the symbolic fixture.
It is unnecessary to expand the large actual output to verify the block proof.

    .venv/bin/python code/s61_binary_nonsquare_schur34_circuit_20261103_v62.py --output work/schur34-circuit.json

# One actual global unit pivot from the full pair dual

Over FULLU, the actual saved projected cubic presentation admits one
additional unit pivot. The eight-step residual presentation
`17x24` is transformed to an explicit finite-DAG presentation
`16x23` with the same cokernel. There are now nine global unit pivots
in this transformed projected presentation. This is not a claim that
an original fixed nine-by-nine minor is a unit, and it does not prove
whole-FULLU cubic rank 57 or full quartic rank 239.

The code and report stem is
`S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63`.
The successful exact calculation took 0.265 seconds. It expanded only
the six-term scalar identity and the necessary row-4 coefficients.
The other sixteen rows remain lazy arithmetic expressions in pinned
checkpoint entries. The final certificate has 2,138 scalar DAG nodes.
Independent review is pending.

## The fixed scalar identity and its actual source

Use the accepted eight-step checkpoint, whose remaining source rows
start with 4, 5 and 9. Row 4 is the transformed `1600*q5*x0` relation.
Its only previous elimination was the `q0*x0` pivot at old column 0,
whose pivot was `25a/4`. In columns 4 through 9, which are
`x0*Sym2(y1,y2,y3)`, its six coefficients are literally

\[
 a_j=(C_L)_{1,j}-\frac{c}{25a/4}(C_L)_{0,j}.
\]

Here `c` is the saved pre-elimination coefficient, and `C_L` is the
accepted L-frame pair `(q0,1600*q5)` on `L=x0=0`. The certificate
checks all six identities directly against both the original augmented
matrix and the normalized checkpoint row 0.

Let `v` be the full second column of the accepted pair dual, after its
monomial reversal and q5 scaling. All four pair-dual identities are
checked again. Consequently the fixed six-term dot is exactly

\[
 \sum_{j=0}^5 a_jv_j=1.
\]

Unlike the earlier failed restriction to three L columns, the whole
six-component dual is used and is transported to a different set of
six target monomials. Its component `v[1]` is literally zero. The
actual coefficients have only powers of `a` (at most two) and
rational constants in their denominators. All twelve denominator
factorizations for the dual and row coefficients are saved; no new
unknown factor is inverted.

## Explicit in-place target completion

Order these six original target columns as
`[4,6,7,8,9,5]`, placing the zero dual component last. Write the row
as `(b,k)` and its dual as `(vbar,0)`, with `b*vbar=1`. Define

\[
 U_0=\begin{pmatrix}I_5-\bar v b&\bar v\\b&0\end{pmatrix},
 \qquad
 E=\begin{pmatrix}I_5&0\\-kb&1\end{pmatrix},
 \qquad U=U_0E.
\]

These are actual finite operators in the original six-dimensional
target, not a claim that an unspecified projective complement is
free. The code checks universal integer polynomial identities in
eleven independent variables. For `t=b*vbar`, they include

\[
 U_0^2-I_6=(t-1)
   \begin{pmatrix}\bar v b&-\bar v\\-b&1\end{pmatrix},
 \qquad
 (b,k)U=((1-t)(1+k)b,t).
\]

Both products with the proposed inverse `E^{-1}U0` are checked with
their exact `(t-1)` remainders. A fixed universal 720-term Leibniz
calculation gives `det(U0)=-t`; this is an integer identity, not a
new source-matrix minor or parameter search. At the actual dot `t=1`,
`U0` is an involution, `det(U)=-1`, and `(b,k)U=(0,1)`.

Row 4 has three further nonzero coefficients outside this block, at
original target columns 1, 2 and 3. After embedding `U` in all 24
targets, perform the three elementary operations

`column_j <- column_j - row4_j * column_5`, for `j=1,2,3`.

They have explicit inverses and require no division. The complete
row 4 is now the unit vector in original target column 5. The target
change has determinant -1. Its full 24-by-24 finite operator and
inverse are saved. Original target columns 27 through 31, namely the
four remaining projected heavy coordinates and psi, are unchanged.
Together they are the five heavy quotient coordinates left after
the old G1 unit pivot; no heavy information is dropped here.

## The actual sixteen-row presentation

Apply the target operator to all seventeen saved rows. For every
other row subtract its transformed column-5 entry times row 4. This
is an explicit source change with determinant one. Since row 4 is
a unit vector, deleting that source row and target column gives
the exact sixteen-by-23 presentation with the same cokernel.

The certificate saves the source operator, the old and new index
lists, and every resulting matrix entry as a finite DAG. The retained
source indices are labels for these transformed relations, not a
claim that a literal old submatrix is being selected. Their original
25-row indices are

`[5,9,10,12,13,14,15,16,17,18,19,20,21,22,23,24]`.

The retained original 32-target indices are

`[1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,21,22,25,27,28,29,30,31]`.

The two remaining quadratic rows still have all-point rank two after
this valid unit elimination. No explicit right inverse or next
completion for those rows is asserted. There is no automatic
Gaussian retry, additional minor, new point or GB in this packet.
The known quartic kernel count 38 and gap 42 are unchanged.

## Preserved implementation failure

The first attempt stopped after the scalar and universal identities
because a metadata assertion incorrectly required all dual denominators
to be rational constants. The accepted coefficients actually have
`a` and `a^2` denominators. The failed source, log and receipt are
preserved under `failed-attempt/`; no final certificate or pivot was
promoted from that attempt. Replacing the assertion by explicit
factor stripping against the already certified FULLU units produced
the successful receipt. Combined execution time stayed below two
seconds, within the authorized single 60-second job budget.

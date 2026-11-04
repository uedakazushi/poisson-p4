# A second global cubic pivot from a three-term Bezout identity

Over FULLU, the actual sixteen-by-23 projected cubic presentation
admits another global unit pivot. An explicit GL3 target change,
all required pure-target tail corrections, and a source clearing
give a finite-DAG `15x22` presentation with the same cokernel.
There are now ten global unit pivots in the transformed projected
presentation. Full cubic rank 57 on all FULLU is still open.

The code and report stem is
`S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63`.
The exact calculation took 0.116 seconds. It expands only the three
actual coefficients and their fixed dot; all remaining row and
tail coefficients stay lazy. The certificate extends the preceding
2,138 scalar nodes unchanged to 3,914 nodes. Independent review is
pending. No Gaussian, minor search, new point, GB or next automatic
operation is part of this packet.

## Actual source identity

The source label 5 identifies the transformed relation originating
from `1600*q5*y1`. Its actual L coefficients at original target
columns 21, 22 and 25 are `(b0,b1,b2)`. They are literal old checkpoint
references even after the preceding GL6 and source operations,
because those operations involve no-L targets and a source row
whose L coefficients are zero. Put

`T = 900a - 126e - 121`.

The code checks the following two exact identities directly from
the saved physical-source coefficients:

\[
 b_0=\frac{aT}{150},\qquad
 2b_1+b_2=\frac{(18e+7)T}{1125}-\frac19.
\]

Consequently

\[
 (b_0,b_1,b_2)
 \begin{pmatrix}6(18e+7)/(5a)\\-18\\-9\end{pmatrix}=1.
\]

This is a fixed literal three-term identity. Its dual has the
constant unit component `-9`, so completion needs only rational
constants and the already inverted factor `a`. It is not inferred
from pointwise rank, a projective-complement assertion or a modular
calculation.

## Finite GL3 and tail operators

Write `v=6(18e+7)/(5a)` and define, in target order `(21,22,25)`,

\[
 U_0=\begin{pmatrix}v&1&0\\-18&0&1\\-9&0&0\end{pmatrix},
 \quad
 E=\begin{pmatrix}1&-b_0&-b_1\\0&1&0\\0&0&1\end{pmatrix}.
\]

The inverse of `U0` is explicitly

\[
 \begin{pmatrix}0&0&-1/9\\1&0&v/9\\0&1&-2\end{pmatrix}.
\]

All entries of both inverse products for `U0`, `E` and `U0*E` are
checked as universal polynomial identities over
`QQ[b0,b1,b2,v]`. If `t=b0*v-18*b1-9*b2`, then

`(b0,b1,b2)*U0*E = (t,(1-t)*b0,(1-t)*b1)`.

The actual dot gives `t=1`. Thus the three-coordinate row becomes
`(1,0,0)`. The determinant of this target change is `-9`, a rational
unit; no unknown core factor is introduced.

The full source row has additional pure coefficients. After embedding
the GL3 operator in all 23 targets, the code clears every such tail
coordinate by subtracting its saved row coefficient times the new
pivot column 21. These fifteen elementary operations are saved as
finite DAG expressions; their coefficients are not expanded. The
five heavy quotient coordinates 27 through 31 have zero coefficients
in the pivot row and are unchanged. They consist of four projected
heavy coordinates plus psi, not six coordinates.

The full pivot row is now the unit vector in target column 21. Its
column is cleared from the other fifteen rows by an explicit source
operator with determinant one. Removing that source row and target
column preserves the cokernel. The certificate saves every target
operator and inverse, every source coefficient, and all `15x22`
matrix roots. Source labels are identifiers of transformed relations,
not a claim that the answer is a literal old submatrix.

The remaining source labels are

`[9,10,12,13,14,15,16,17,18,19,20,21,22,23,24]`.

The remaining original 32-target indices are

`[1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,22,25,27,28,29,30,31]`.

There is one remaining quadratic source row, originating from
`q5*L`, and fourteen genuine source rows. Its all-point rank-one
result after the valid unit eliminations does not supply a computed
unit functional or next target completion. No further pivot is
claimed. The accepted cubic57 result on D(Delta34), full quartic239,
and the quartic kernel count 38 with gap42 are unchanged.

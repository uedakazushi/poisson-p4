# The final quadratic-row completion with one spare pure target

The actual finite-operator completion passed in 0.092 seconds after the
independent POWER review gate. It gives a `14x21` presentation with the
same cokernel as the actual parent `15x22` presentation. The code takes
the frozen review of
`S61_NONSQUARE_CUBIC_LAST_Q9_POWER_BEZOUT_20261104_v63` as an explicit
input and checks its binding to the actual POWER certificate. The
review report hash is
`8cfcd14f1c55202f36c1b9730e75763f2b2a33659ba0b4f80d1b7066aaab314c`.
Independent review of this new completion is pending.

The parent actual presentation is `15x22`. Its first row is the final
quadratic relation labelled 9. The nine coordinates at target labels
`[11,13,14,15,17,18,19,22,25]` are denoted h. The accepted POWER identity
gives the actual dual `c_i=H_i/a^29`, so `h*c=1`. The other thirteen
coordinates b of the full row are retained, not set to zero by assumption.

First apply the elementary target matrix

`T=[[I9,-c*b],[0,I13]]`.

Its inverse changes the minus sign to a plus sign, independently of the
unit identity: the off-diagonal block has square zero. The full pivot
row becomes `(h,0^13)` by the pinned actual inner product. Choose original
target label 1 as the spare zero coordinate. It is a pure target, so no
heavy coordinate is used as the spare.

On the block ordered as the nine h labels followed by label 1, apply

\[
 U=\begin{pmatrix}I_9-ch&c\\h&0\end{pmatrix}.
\]

For formal `t=h*c`, the universal polynomial identities are

`U²-I=(t-1)*[[c*h,-c],[-h,1]]`,

`[h,0]*U=[(1-t)*h,t]`.

The code verifies all 100 matrix coefficients and all ten
row coefficients in the universal integer polynomial ring. At `t=1`,
U is an involution and maps the full row to the unit vector in the
spare target coordinate.

The determinant argument does not assume a free basis of `ker(h)`.
For `B=I-ch`, column multilinearity leaves only the identity term and
nine singleton rank-one terms, giving `det(B)=1-t`. The adjugate is
`(1-t)*I+ch`; both matrix products and the bordered expression are
checked as polynomial identities. Hence

`det(U)=-h*adj(B)*c=-t`,

and the actual target determinant is -1. This avoids a factorial-size
determinant expansion.

Embed U into the full 22 targets and compose it after T. The inverse
is `U*T^-1`. All thirteen tails are kept as finite circuit references.
No degree-207 cofactor is expanded: the POWER DAG and preceding matrix
DAG remain external pinned inputs. Clear the unit column from each
other relation using an explicit source matrix of determinant one,
then remove that source row and target column. The resulting actual
presentation is `14x21` with the same cokernel.

All five remaining heavy quotient coordinates, labels 27 through 31,
are unchanged. The quadratic pivot row has zero coefficients there,
the tail clearing is therefore identity on those columns, and the
spare target is pure. The code compares their final
literal scalar roots to the parent roots.

The retained source labels are the fourteen genuine rows. There
are eleven global unit pivots in the transformed augmented
25-row presentation: ten quadratic pivots and the previously split
genuine heavy pivot. This is not a particular original eleven-minor,
and it is not full cubic rank 57. No further Gaussian, point, minor,
GB or automatic next operation is part of this recipe. Full quartic239
and kernel38/gap42 remain separate.

The complete new certificate has 6,880 scalar nodes and is 113 KB.
It refers to the external 25,572-node POWER circuit and the preceding
3,914-node matrix circuit without expanding their arithmetic. The
new inverse denominator is only the certified unit `a^29`. It saves
the full 22-target change and inverse, all 15 source clearings,
all 294 final `14x21` coefficient roots and all source/target labels.
The retained source labels are
`[10,12,13,14,15,16,17,18,19,20,21,22,23,24]`; the retained original
target labels are
`[2,3,4,6,7,8,9,11,13,14,15,17,18,19,22,25,27,28,29,30,31]`.
They identify the transformed relations and targets, not a literal
old submatrix.

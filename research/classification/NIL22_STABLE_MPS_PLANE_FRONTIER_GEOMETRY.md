# Stable sources and all listed rational MPS plane branches

2026-11-07. Hand proof from the frozen actual plane formulas. No new CAS or anchor is used. This is a statement about the listed uniform normal-cone branches, not the full MPS components or an exhaustive classification of all Lagrangian planes.

Write the stable divisors as `H²*C3`. For stable_line_conic, `C3=LQ`, Q is smooth, L is tangent to Q, H is a secant, and `H∩L∉Q`. For stable_three_lines, `C3=L1L2L3` and the four distinct lines H,L1,L2,L3 have no triple intersection. In both cases H meets C3 in three distinct points and H avoids every singular point of C3. The six magnetic parameters do not change these configurations.

## 1. T37 cell013 prime1: a repeated intersection with H

The existing actual factor formula is

\[
H=At^2-2B,
\]
\[
\begin{aligned}
C_3={}&t^2A^3-(6t^5+2)A^2B+12t^3(t^5+1)AB^2
       -8t^6(t^5+3)B^3\\
     &-144t^9AC^2-96t^9B^2C.
\end{aligned}
\]

Substituting `B=At²/2` gives the exact identity

\[
C_3|_H=-t^9 A\,(t^4A+12C)^2.                   \tag{1}
\]

Thus the binary cubic restriction has zero discriminant. This condition is closed in the proper projective factor space `(H,C3)`, including zero restrictions. Its multiplication image `H²*C3` is closed. At either stable target, the unique double H must be the limiting squared factor, so the remaining cubic is precisely its C3. Its restriction to H has three distinct roots, contradicting (1). This excludes the whole projective closure of this actual branch, including its rational-chart parameter boundaries.

## 2. T36 cell012: a common singular point of Q2 and C3

The complete actual plane formulas give

\[
Q_2=5tB^2-tuBC+25C^2,
\qquad C_3=3A Q_2+25B^3+5tC^3.
\]

At `p=[1:0:0]`, both Q2 and C3 are singular. Existence of such a point p is a closed incidence in the proper projective space of `(p,Q2,C3)`: every first derivative of each homogeneous form vanishes at p. Its multiplication image under `(Q2,C3)↦Q2*C3` is closed.

For the stable line-conic target, a degree-two factor Q2 can limit only to Q, HL or H². It cannot be Q because Q is smooth. If Q2=HL, its singular point is `H∩L`, but the remaining cubic HQ is smooth there since `H∩L∉Q`. If Q2=H², its singular locus is H, while the remaining cubic LQ is singular only at its tangency point, outside H. All assignments contradict the common-singular-point incidence.

For the stable three-lines target, Q2 is H², HLi or LiLj. In the first case every possible p lies on H, but the cubic L1L2L3 has no singular point on H. In the second, p=H∩Li and the remaining cubic HLjLk is smooth at p. In the third, p=Li∩Lj is outside H and Lk, so the remaining cubic H²Lk does not even vanish at p. Again every assignment is excluded. This proves the result for the whole closed actual branch, not just generic t,u.

The quadratic Q2 here is binary, hence geometrically reducible over C. Its degree-two grouping in the original QQ(t,u) factor report is not an assertion of geometric irreducibility; the present argument uses its exact singular-point structure.

## 3. T38 cell013 prime1: a cubic-cone vertex on Q2

The report gives

\[
\begin{aligned}
Q_2={}&A^2(t^4u-t^2)+AB(t^4v-2t^3u+4t)+ACt^5\\
 &+B^2(4t^5-2t^3v-4)-4BCt^4,
\end{aligned}
\]

and a cubic C3 depending only on A,B. Thus at `p=[0:0:1]`, C3 has multiplicity at least three and Q2(p)=0. Existence of such p is a proper closed flag incidence: `C3∈I_p³`, `Q2(p)=0`.

At a stable line-conic target, Q2 must be the irreducible Q and C3 must be H²L: alternative cubic groups contain the irreducible conic and cannot be a binary cubic cone. The unique multiplicity-three point of H²L is `H∩L`, which is outside Q, contradiction.

At a stable three-lines target, a cubic group of three distinct lines among the four has no multiplicity-three point, by the no-triple-intersection condition. The only remaining possibility is C3=H²Li and Q2=LjLk. Its vertex is `H∩Li`, which lies on neither Lj nor Lk. This again contradicts the incidence. All factor assignments and projective boundaries are included.

## 4. The other listed rational branches

`NIL22_STABLE_MPS_PLANE_INCIDENCE.md` handles T36 cell013 prime2 and T37 cell013 prime2 by their exact `Q=AL+γH²` identities: the simple line's conic tangency is forced to lie on the double H, unlike either stable target.

Every other rational branch listed in `CLASSIFICATION_GA9_MPS_PLANE_DIVISORS_20261107.json` has either a linear factor of multiplicity at least three, or two distinct squared linear factors with a fifth linear factor. These are proper closed multiplication conditions and neither stable target satisfies them: both have exactly one double linear factor and no triple linear factor. This assertion is only about the explicit listed branch families; it does not assume the list is an exhaustive Grassmann decomposition. The explicit T38 cell013 prime2 branch belongs to the two-squared-linear case. So does T38 cell012: its squared quadratic is the binary form `tB²+uBC+2C²`, which splits into two linear factors over C; coincident factors only increase the multiplicity.

Consequently **all rational plane families explicitly listed in that frozen report fail to supply either full stable divisor configuration by this uniform normal-cone construction**, including limits inside each family's projective factor incidence. This controls all magnetic coefficients because the obstructions concern their fixed actual divisor configurations.

The T37 nonrational generic lifted B branch has the separately proved `LQ²` identity and is likewise excluded from these stable configurations by its closed squared-quadratic factor condition. Its finite denominator complement is an actual20-dimensional QQ algebra and remains unresolved in `T37_OPEN_LAGRANGIAN_CURVE_FRONTIER.md`; it is not included by extrapolation.

No claim is made that stable sources fail to belong to T36/T37/T38 themselves. Other plane components, moving genuine periods, nonuniform contractions, or actual families unrelated to this plane construction remain possible. The full global T39/T40 graph theorem is a separate argument, not a consequence of the present divisor flags.

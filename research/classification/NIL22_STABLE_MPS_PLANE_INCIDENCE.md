# Two actual MPS plane branches have the wrong stable incidence

2026-11-07. Hand geometry from frozen actual factor identities. This excludes two specified uniform normal-cone source branches, not their entire MPS components or all moving-period contractions.

The complete stable sources have divisors `H²*L*Q`, where L is tangent to the smooth conic Q but `H∩L` is outside Q, and `H²*L1*L2*L3`, where the four distinct lines have no triple intersection. These configurations and their independence of all magnetic coefficients are already established in `NIL22_STABLE_PROJECTIVE_ACTION.md` and `NIL22_STABLE_CLOSED_FRONTIER.md`.

The frozen report `CLASSIFICATION_GA9_MPS_PLANE_DIVISORS_20261107.json` gives full actual Lagrangian-plane matrices and factor identities for the following two branches.

## The two whole factor identities

For T36, cell013 prime2, let

\[
H=B,\quad L=5Bt^2+C,\quad
Q=3(tu-5)A L+(-25t^5+5tu-25)H^2.
\]

The actual restricted quintic is a nonzero rational factor times `H²LQ`. The displayed Q is exactly the report's quadratic factor, since its A terms are `3(tu-5)A(5Bt²+C)`.

For T37, cell013 prime2, let

\[
H=At+2B,
\quad L=A(2t^6+t)+B(6t^5+2)-12Ct^4,
\quad Q=AtL+(1-t^5)H^2.
\]

Expanding the latter identity gives the report's quadratic factor

\[
A^2(t^7+2t^2)+AB(2t^6+6t)-12ACt^5+4(1-t^5)B^2.
\]

Again the actual restricted quintic is a rational factor times `H²LQ`. Thus in both entire generic branches

\[
Q|_{L=0}=\gamma\,(H|_{L=0})^2.                 \tag{1}
\]

Whenever Q is smooth and H,L are independent, L is tangent to Q precisely at `H∩L`. In particular the tangent point lies on the double hyperplane.

## Closed incidence, including parameter boundaries

Parameterize H,L,Q projectively. Impose that the restriction of Q to the line L has rank at most one and is proportional to the square of the restriction of H. Include zero restrictions: equivalently one may use the proper incidence of a point p with `p∈H∩L` and Q having intersection multiplicity at least two with L at p, allowing Q to contain L. This is closed in the proper flag/factor parameter space. Its multiplication image under `(H,L,Q)↦H²LQ` is closed. Formula (1) therefore holds as this necessary factor-incidence condition on the whole projective closure of both specified actual divisor families, including parameters excluded by their rational chart descriptions.

For the stable line-conic divisor, the irreducible Q must be the limiting quadratic factor. The squared H must limit to the unique double factor H, and L must limit to the remaining simple line: other assignments would introduce multiplicity at least three or lose the irreducible quadratic. The required tangency point is then `H∩L`, contradicting the established stable configuration `H∩L∉Q`.

For the stable three-lines divisor, the squared H again limits to the unique double line. L limits to one simple line, and Q to the product of the other two. Equation (1) forces both other lines to meet L at `H∩L`. This gives a triple intersection (in fact all four lines meet there), contradicting the stable four-line configuration. These cover all possible factor assignments at the targets.

Consequently neither full stable source can be the endpoint of the specified two uniform normal-cone branches. This is independent of all six magnetic coefficients because the divisor configurations do not change with those coefficients.

## Remaining source geometry

This does not show that T36 or T37 as whole components fail to contain a stable source. Nonuniform contractions, other Lagrangian-plane branches or varying genuine source parameters are not covered.

A different existing branch remains meaningful: T37 cell013 prime1 has actual divisor `H²*C3(t)` with `H=At²-2B`, and its cubic factor is explicitly retained in the same frozen report. A complete one-parameter factor/marked-incidence analysis of that cubic would determine whether its special fibers have either stable configuration. No such factor classification or new anchor computation is performed in this note.

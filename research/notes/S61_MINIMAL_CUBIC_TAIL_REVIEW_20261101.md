# Independent review of the M1 cubic-tail local theorem

Date: 2026-11-01. Accepted at the specified point and the justified open: the completed full135 homogenized local ring is C[[u,e]]/(e³), and the homogenizing parameter t is zero. No endpoint, intersection or whole-family saturation conclusion is added.

## Exact independent checks

The review script does not import the main local elimination or specialization routines. It reads the original134 equations and the saved universal G6 polynomial, independently collects their full homogenized coefficients, and checks all135 equations.

At the canonical tensor it independently obtains rank six for the original beta-linear equations. The saved16 by10 matrix has rank ten and is their full kernel. The extra G6 restriction vanishes on that kernel. The44-column leading Jacobian is reconstructed over QQ(b0,...,b9), and all entries agree with the certificate. Its two saved tangent vectors are independent, satisfy D*T=0, and an independent rational specialization has rank42. Thus the generic rank is exactly42. The tangent quadratic coefficient matrix H is independently reconstructed, and D*L=H is verified for all three quadratic monomials on the full ten-parameter beta kernel.

At the specified point, the full homogenized Jacobian is rebuilt, including the t column from the lower q-degree terms. Its rank is43 and its leading44-column rank is42. The recorded43-dimensional minor is recomputed and nonzero. The review then substitutes the saved normal tangent and second/third corrections directly into every original homogenized polynomial, as truncated one-variable series. Orders zero, one and two vanish in all135 equations; the third-order residual agrees entrywise, and original row66 has coefficient84035/389017. The eliminated43 rows vanish at order three, and the corrections leave the free coordinates fixed, as required for the implicit-function slice.

The same fixed-beta curve is checked against every leading equation, with its projective tangent and its nonzero chosen u coordinate. Finally the finite parameter-dependent normalization is checked as an actual symplectic tensor transformation for all a,b,d. Its denominator along the curve is3+7z, nonzero at the specified point. These are full coefficient identities, not numerical tests.

## Why the universal quadratic identity applies along this curve

The identity D*L=H is for the leading equations themselves, so it computes the quadratic obstruction of S/(t), not merely that of the full S. This distinction is essential to the lower bound used later.

On the rank42 locus in the irreducible ten-dimensional beta kernel, the kernel and cokernel of D are vector bundles. The quadratic obstruction is an intrinsic bundle map Sym² ker(D)→coker(D). Its vanishing on the dense open where the displayed rational basis is defined extends over the entire rank-constant open. This is a vanishing statement about a bundle map, so a pole in the chosen tangent basis does not introduce an exception to it.

The exact symplectic normalization transports an open neighborhood of the same actual curve to the canonical tensor. The previously proved full135 homogenized equivariance is a stated dependency. Rank42 holds at the specified point and hence on a nonempty open of that curve. Consequently the generic curve point has relative embedding dimension one and no quadratic relation in S/(t).

## Finite freedom and the equality t=0

The transverse full-homogenized slice has one tangent variable, no quadratic relation, and a nonzero cubic relation. Formal elimination by the invertible43-minor therefore gives length exactly three, not merely an upper bound inferred from a truncated computation: a nonzero cubic coefficient is a unit times v³ in the one-variable completed ring, and no relation of lower order exists.

Choose u to restrict to a uniformizer on the embedded smooth curve and put R=C[[u]]. The finite length of S/(u) makes the u-adic and maximal-ideal topologies cofinal. Complete Nakayama then gives a surjection R³→S, so S is finite over R with rank at most three.

Its quotient S_infinity=S/(t) has at the generic point of that same curve an Artin factor of length at least three, by the relative embedding dimension and universal absence of quadratic relations just proved. Finiteness is already available at this step. Thus

\[
3\geq\operatorname{rank}_R S\geq
\operatorname{rank}_R S_\infty\geq3.
\]

The kernel of R³→S has rank zero and lies in a free module over the domain R, so it is zero. Hence S is free of rank three. Its ideal (t) also has rank zero and lies in the torsion-free module S, giving (t)=0. This proves exact vanishing; the observed vanishing of the first few series coefficients of t is not used as a substitute.

The reduced curve provides a retraction S→R. Choose e in its kernel so that u,e generate the completed local algebra. The generic factor supported on that curve already uses the full rank three, so there is no other generic support. In the resulting length-three local algebra over Frac(R), e is nilpotent and e³=0. Torsion-freeness extends that equality over R. The surjection R[e]/(e³)→S is between free rank-three R-modules and is therefore an isomorphism. This argument uses the unique generic support; freeness alone would not justify discarding lower coefficients of a cubic relation.

Faithful flatness of completion descends t=0 to the original local ring. Saturation by t therefore removes the specified point. The same proof applies where its explicit rank, curve and cubic-coefficient conditions hold, while exceptional points of the rational curve remain outside this conclusion.

## Reproduction and dependencies

`code/s61_minimal_cubic_tail_review_20261101.py --output <new-path>` produces `reports/S61_MINIMAL_CUBIC_TAIL_REVIEW_20261101.json`. It pins the main code/report/note, original equations, universal G6 and homogenized-equivariance proof. It independently verifies the full-kernel obstruction identity and the direct order-three substitution; it does not reprove the earlier global equivariance theorem or the construction of G6.

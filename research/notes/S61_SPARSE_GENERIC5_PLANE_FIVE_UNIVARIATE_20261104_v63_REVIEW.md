# Independent review of the five-row formal-curve selector

The accepted sparse5 input/plane and full134 circuit reviews bind the original
134 equations, center C, original plane directions W0,W1 and
d=Q52-C. This audit reconstructs only rows 50,10,52,51,15 directly from that
source under the ring homomorphism
(H,L,s0,q0,t0)=(T,T²,T³,T⁵,T⁷). It does not use or rebuild the earlier
692953-term coefficient matrix or repeat its guarded field elimination.

For each selected original equation F, a new C++ checker computes, over QQ[T],

    [P F(C), J_F(C) W0, J_F(C) W1,
     polar_QF(W0,d), polar_QF(W1,d)].

Every linear and quadratic q monomial of these five original equations is
included. Diagonal q squares use both cross terms in the polarization. All
25 results are compared with the frozen selected matrix. Clearing one
integer denominator per row gives an integer polynomial 5-by5 matrix; its
exact polynomial-matrix determinant, computed independently of the author's
field elimination/Bareiss trace, agrees with the saved nonzero degree448
polynomial. A nonzero image under a ring homomorphism proves that the original
five-variable determinant was nonzero. No numerical beta, source point or
finite prime is used, and no assertion about every value of T is made.

Here is the separate five-functional argument. The exact trace corollary
Tr(c)=0 gives the centroid of the three generic toric points as
C=(2Q51+Q52)/3. In plane coordinates the other two points have coordinates
(0,0,-1/2) and (0,0,1). The toric triple is generically noncollinear. After
an invertible linear change in its two coordinates it is
(1,0,0),(0,1,0),(-1,-1,0). The audit checks the universal five-point
evaluation matrix has rank5 and exhibits five quadrics whose values of the
functionals (constant,linear_a,linear_b,a_sigma,b_sigma) form I5. Thus those
functionals are injective on the full five-dimensional vector space of
degree-at-most-two polynomials vanishing at the five points. Under any GL2
change of the triangle coordinates, the selected two linear and two mixed
functionals change by invertible blocks, preserving injectivity.

The actual directions are V0=W0/P,V1=W1/P,d, so the five displayed columns are
P times precisely these five functionals. The known generic five points
are actual solutions, so every restricted original equation belongs to this
five-dimensional space. Noncollinearity and distinctness are invariant open
conditions: the accepted generic configuration and the dominant Sp4 times
physical-Gm sparse5 cover transfer them to this source. This uses invariant
genericity, not an assertion that an arbitrary fixed non-invariant minor
remains nonzero on a slice. The independent formal-curve determinant proves
rank5 and identifies five actual independent restricted equations.

This establishes the plane-restriction rank5 only. It gives no containment
of the full affine S61 fiber in the plane, no cross126 rank, finiteness or
exhaustion. The preceding failed field job and its partial receipts remain
unchanged. The review executable stays in ignored work/.

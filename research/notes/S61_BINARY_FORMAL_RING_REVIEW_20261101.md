# Independent review of the rank65 binary formal-ring theorem

Date: 2026-11-01. Accepted at the specified point and on the justified constant-rank open. For R=C[[v0,v1,v3]], the full135 completed local ring S is R-free of rank65, its nilradical I satisfies I⁴=0, and its associated graded pieces are free of ranks(1,7,23,34). Moreover 0≠t∈I³, It=0, t²=0, and S/(t) is R-free of rank64. No assertion that the entire multiplication table is constant over R is needed or accepted.

## The same actual fixed-beta threefold

The review independently reconstructs the full135 Jacobian and tangent quadratic forms from the original equations and universal G6, without the main specialization or primitive-beta helpers. It confirms leading rank34, full rank35, ten tangent directions, five quadratic obstructions and a three-dimensional common Hessian radical.

It also reconstructs the old closed first-jet form directly from the projective S61 bivector in the affine chart. Its restriction to an isotropic two-plane agrees with the primitive form, so it gives the same two curve equations. Their derivative is independently obtained as

\[
\begin{pmatrix}22&-21&2\\12&-27&3\end{pmatrix},
\]

with kernel direction(1/38,7/57,1). Thus this fixed-beta curve is smooth at the stated plane. Adding the two binary quadratic-form parameters gives the actual three-dimensional family in the main proof. The review differentiates its tensor formula, normalizes the five directions in S61, makes the projective-chart correction and obtains exactly the saved three tangent columns. They lie in the full135 tangent kernel and span precisely its common Hessian radical. The minor on coordinates v0,v1,v3 is -1/19.

The family is embedded locally: its nonzero tensor coefficients recover the common support plane as their rank-two image, and then recover the binary form. On the non-square open these are regular reconstructions. This is stronger than merely locating a three-dimensional subspace of the tangent space. The prior universal R0/binary identities provide inclusion of the whole family in t=0; this review does not claim a separate recomputation of those full parameter identities.

## Uniform lower bound on that same family

The earlier universal rank, quadratic and cubic certificates are explicit dependencies. The review checks their stated bounds and the t-free support of the quadratic lift; it does not rerun the full universal Bianchi construction.

Their passage from the dominant slice to the constant-rank incidence is valid. On an open where the Jacobian, quadratic obstruction rank and the smooth threefold tangent dimension are fixed, kernels, cokernels and the relevant normal spaces are vector bundles. Choosing local invertible minors expresses the third-order initial relation matrix regularly. Its rank-at-most50 condition is closed, and the verified dense slice, together with full135 equivariance and irreducibility of the binary incidence, forces it throughout that open. The same reasoning applies to the vanishing of the first and second normal initial classes of t. This is not an inference from an arbitrary generic-beta sample to every special fiber.

Consequently the general point of the same fixed-beta threefold has seven normal tangent directions, at most five quadratic relations and at most fifty cubic initial relations. Its truncated normal length is at least1+7+23+34=65. The exact special transverse quotient has length65. These are bounds on the same local problem, with the actual parameter minor established above.

## Finite freedom over a three-dimensional power-series ring

The finite special quotient implies that the R-maximal-ideal topology and the S-maximal-ideal topology are cofinal. Complete Nakayama gives a surjection R^65→S. The smooth embedded threefold supplies a split R-algebra retraction S→R, so R injects into S.

After extension to K=Frac(R), the local factor corresponding to the threefold has length at least65; the generating bound gives total dimension at most65. Equality follows, and there are no other generic factors. The kernel of R^65→S is rank zero inside a free module over a domain, hence zero. This proves freeness over R without assuming that R is a discrete valuation ring; its dimension three causes no difficulty.

Let I be the kernel of the retraction. Then I is free of rank64. The generic length bound uses all65 dimensions through degree three, so the generic Hilbert function is exactly(1,7,23,34) and I_K⁴=0. Since I⁴ is a submodule of the torsion-free module S, it follows that I⁴=0 globally. Thus I is the nilradical and S_red=R.

## Why the graded pieces are free

This step needs more than flatness of S or smoothness of its reduced support. The retraction gives S=R⊕I as R-modules, hence the special fiber of I is exactly the maximal ideal of the transverse Artin algebra A. It follows that

\[
(I/I^2)\otimes_R\mathbb C\simeq\mathfrak m_A/\mathfrak m_A^2.
\]

This module has seven special generators and generic rank seven. A surjection R^7→I/I² has rank-zero kernel in a free module, so I/I² is free. The sequence splits, making I² a free direct summand of I. Its special fiber now injects into A and equals m_A².

Repeat for I²/I³: its special dimension and generic rank are both23, so it is free. The splitting makes I³ free of rank34. Together with I⁴=0 this gives the claimed graded free ranks. The initial equality of special fibers is valid without first assuming that powers of I commute injectively with base change; the successive splittings establish the later injectivity.

## The exact location and nonzero value of t

The universal lift has no t component in either first or second normal order. Thus t has zero classes in I_K/I_K² and I_K²/I_K³. Changing the implicit-function complement alters a second lift by tangent directions, on which t is already zero; equivalently its quadratic class is a section of the regular quadratic quotient bundle. The dense-slice identity therefore gives t∈I_K³ on the constant-rank open used here.

The graded-freeness argument makes S/I³ R-free. Its class of t is generically zero, so it is zero everywhere: t∈I³. Hence It=0 and t²=0. The exact special transverse algebra has a nonzero cubic socle class represented by t, so t is not zero in S.

Because It=0, the ideal(t) is Rt. In the free module S, the nonzero special class of t means that one of its coordinates is a unit. Therefore Rt is a free rank-one direct summand and S/(t) is free of rank64. This does not confuse nonzero square-zero t with t=0.

Faithful flatness of completion descends t²=0 to the original local ring. Localizing by t gives zero, so the specified point is removed by saturation. No conclusion is drawn about the unexamined exceptional binary strata or the full affine fiber.

## Reproduction

`code/s61_binary_formal_ring_review_20261101.py --output <new-path>` produces `reports/S61_BINARY_FORMAL_RING_REVIEW_20261101.json`. It independently reconstructs the point135 Jacobian, quadratic cokernel and actual fixed-beta family tangent. It pins the previously accepted transverse and universal certificates, whose full calculations remain explicit dependencies of this formal-ring assembly.

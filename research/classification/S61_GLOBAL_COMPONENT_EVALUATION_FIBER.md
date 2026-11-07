# Connected evaluation fibers identify the full normalized T39/T40 components

2026-11-07. PRIMARY E. New hand proof, independent acceptance pending. This concerns reduced supports and actual global components; it does not classify scheme thickness or the other S61 components.

## 1. The global evaluation fiber is irreducible

Work over C. Let `W=H^0(P4,Λ²TP4)`, an affine space of dimension 126. Let C be a global irreducible projective Poisson component of generic rank four, and let X be its affine cone. A connected algebraic group preserves each irreducible component: the orbit of a component in the finite set of components is constant. In particular X is stable under PGL5 and under the GL4 subgroup `diag(1,A)` fixing `p=[1:0:0:0:0]` and the hyperplane at infinity.

Evaluation at p is a linear GL4-equivariant map

\[
\mathrm{ev}_p:W\longrightarrow\wedge^2T_pP^4\simeq\wedge^2\mathbb C^4.
\]

Let U be the six-dimensional open set of nondegenerate skew tensors. The preimage `X_U` is a nonempty irreducible open subset of X. Nonemptiness follows by moving a rank-four point to p. GL4 is transitive on U, with stabilizer Sp4 of the fixed tensor J. Sp4 is connected.

For clarity, the connected stabilizer is essential to the next assertion. Equivariance identifies

\[
X_U\simeq\mathrm{GL}_4\times^{\mathrm{Sp}_4}Y,
\qquad Y=X\cap\{\mathrm{ev}_p=J\}.
\]

This is the associated-bundle isomorphism, checked after the faithfully flat locally trivial torsor `GL4→GL4/Sp4`. Sp4 preserves each irreducible component of Y. If Y had multiple components, each component would descend to a distinct closed component of the associated bundle, contradicting irreducibility of `X_U`. Hence **Y is irreducible**. This argument concerns the whole fiber, not only its generic smooth locus. Equivalently, it applies to reduced supports after base change to C and does not assume a flatness theorem for an unrelated family.

If `dim C=d`, then `dim X=d+1` and

\[
\dim Y=d+1-6=d-5.
\]

All fibers over U are isomorphic by GL4, so no generic-fiber dimension estimate is being substituted at the special value J.

## 2. The four perspective pivots give a global polynomial product

Use affine spatial coordinates `z=(z1,z2,z3,z4)` near p and standard `J12=J34=1`. For any global bivector of value J, its affine expansion is `J+L(z)+higher terms`. The full projective bivector can have affine cubic terms; these are retained.

For a row `a=(a1,a2,a3,a4)`, the actual projective matrix

\[
g_a=\begin{pmatrix}1&a\\0&I_4\end{pmatrix}
\]

fixes p and has identity tangent action there. Thus it preserves `ev_p=J`. Its affine coordinate action is `z↦z/(1+az)`. If `v=Ja^t`, the change in the affine linear jet is exactly

\[
L\longmapsto L+T_a,\qquad
T_a=-2(az)J+zv^t-vz^t.                         \tag{1}
\]

The higher terms of the original bivector cannot affect this first-jet formula. It is valid for every bivector of value J, without any nonvanishing condition on L or its invariants.

The saved Euler quotient kills every `P0i` constant and linear affine coefficient: its vertical-pivot list contains raw indices `0..4,15..19,30..34,45..49`. Hence the saved spatial `P12` linear coefficients are precisely the coefficients of affine `L12`. The four saved perspective pivots are raw indices `[61,62,63,64]`, i.e. the coefficients of `x0*x1,...,x0*x4` in `P12`. Formula (1) gives

\[
(T_a)_{12}=-3a_1z_1-3a_2z_2-2a_3z_3-2a_4z_4.
\]

Therefore their response matrix is the constant invertible matrix

\[
R=\operatorname{diag}(-3,-3,-2,-2),\qquad\det R=36.
\]

Given any point of the entire J fiber, choose uniquely `a=-R^{-1}ℓ`, where ℓ is its four pivot coefficients. Applying `g_a` sets all four pivots to zero. Matrix multiplication in this unipotent group is addition, so this normalization is constant on each orbit and has the actual inverse `g_{-a}`. Both the normalization and its inverse are polynomial maps on coefficient spaces: the projective representation of `g_a` is polynomial, and `R^{-1}` has rational constant entries. Passing to the fixed linear Euler quotient is also polynomial.

Consequently, for the entire J fiber and for each GL4/PGL5-stable component cone,

\[
Y\simeq\mathbb A^4\times Z,\qquad
Z=Y\cap\{\ell=0\}.                            \tag{2}
\]

In particular Z is irreducible and has dimension `d-9`. No β, covariant, modular discriminant, nilpotent-type or special-period open was removed. This is stronger than a generic or formal slice statement.

Formula (1), including the actual finite perspective action and its radial cubic cancellation, is independently derived in `code/s61_weight1_audit.py`. The saved four pivot indices and all Euler pivot indices are in the immutable v02 certificate. The current `S61.normalize` uses exactly `-R^{-1}ℓ`, with R computed by the actual projective pushforward and the Euler quotient.

## 3. The S61 equations describe the whole perspective section

After quotienting Euler24 and fixing the six constant coefficients and the four perspective pivots, the remaining affine coefficient space has dimension 116. The v02 construction eliminates 55 Jacobi variables in positive-weight order. Each selected equation has coefficient exactly one on its new pivot variable and has only previously solved lower-weight variables in the remaining terms. Thus these are **global polynomial triangular eliminations**, not divisions by a jet invariant or formal implicit solutions.

The saved eliminated-variable formulas have rational constant coefficients and degree at most three. All remaining Jacobi equations are substituted exactly; the 134-generator system uses exact Bianchi reductions with no parameter inversion. Hence the actual Poisson locus in this entire perspective section is globally isomorphic to the saved S61 scheme in A61. Nilpotent or zero first jets are included.

This is the explicit scope of `legacy/v02/code/build_symplectic_slice.py`: “This constructs an equivalent 61-variable scheme, not a prime decomposition.” Its weight-ordered pivot construction is in `build_universal_slice.py`. `audit_reconstruction.py` independently reconstructs the bivector first and recomputes Jacobi coefficientwise. `S61_REBUILD_20261030.json` records exact regeneration of all saved elimination/reconstruction files. The absence of a global *component exhaustiveness* claim there does not restrict the stated polynomial coordinate equivalence.

We identify Z in (2) with its closed reduced subset of S61. This identifies actual coordinate sections, rather than merely matching tangent spaces.

## 4. Application to the two known 25-dimensional components

The accepted universal-section results state that

\[
\Gamma_{51}=\{q=Q_{51}(\beta)\}\simeq\mathbb A^{16},
\qquad
\Gamma_{52}=\{q=Q_{52}(\beta)\}\simeq\mathbb A^{16}
\]

are the reduced S61 closures of T39 and T40 respectively. These are globally defined polynomial graphs for **every** β. In particular every point of each graph belongs to the corresponding global component: the graph is the closure of actual normalized points of that component, and the component is closed. No exceptional first jet has been discarded. These precise statements are in §1 of `Q51_UNIVERSAL_SECTION_20261101.md` and §1 of `Q52_UNIVERSAL_SECTION_20261101.md`, with the separate universal-section review reports.

For either component C, `d=25`. By §§1–3 its full normalized intersection Z with S61 is irreducible of dimension 16. Its corresponding graph Γ is closed, irreducible of dimension 16 and contained in Z. Therefore

\[
C_{39}\cap S61=\Gamma_{51},\qquad
C_{40}\cap S61=\Gamma_{52}                         \tag{3}
\]

as reduced supports in the actual normalized J section. In particular this includes every special/nilpotent/zero β and every component-boundary point having rank four at the selected evaluation point. We make no assertion about embedded structure or points whose tensor has generic rank at most two.

For an arbitrary rank-four tensor, choose any actual rank-four evaluation point, move it to p, Darboux-normalize its value to J and apply the global perspective normalization. Since all these are invertible actual projective transformations, (3) is an if-and-only-if test for membership in the two global components. One does not have to search for another evaluation point or frame after a negative graph test.

## 5. The existing E6 smooth-flex computation now excludes both parents

The frozen report `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_E6_FLEX_GRAPH_20261107.json` contains an actual canonical E6 smooth-flex tensor, the structural rank-four evaluation point with Pfaffian `-75/376`, the full invertible PGL matrix and its inverse, all150 reconstruction, all134 equations and all45 graph residuals for each graph. Each graph has 15 nonzero residuals.

Subject to independent acceptance of the new global-fiber proof above, that existing calculation proves that the canonical23 orbit is outside **both full global components T39 and T40**, independently of frame or evaluation-point choices. The whole dense nonzero-magnetic-invariant orbit is therefore outside both. Its boundary with zero invariant can still meet those components; no exclusion of every source-boundary point is claimed.

The older diagnostic report correctly limited its conclusion to a faithful graph entry because it did not yet contain §§1–4. Its bytes and scope remain immutable. The stronger parent exclusion is a new theorem with those additional hypotheses/proofs. It does not identify a different parent, prove a new component dimension or change the component count. The A5 source has not undergone an additional graph computation here.

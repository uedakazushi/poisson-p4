# Fixed T37, open Lagrangian planes: the omitted curve and its finite complement

2026-11-07. PRIMARY E. Producer exact identities plus explicit unfinished finite source. No all-parent exclusion or component count follows.

## Whole cell and actual source

The accepted complete T37 tensor is used at actual `p0=t0=t1=t2=t3=1`. In Schubert cell012 with row order `(0,3,4,1,2)`, write

```
x0=A, x3=B, x4=C,
x1=a0*A+a1*B+a2*C,
x2=a3*A+a4*B+a5*C.
```

The plane has rank three everywhere in this affine cell. Its intrinsic Lagrangian condition is all six coefficients of `{ell1,ell2}|plane=0`, with `ell1=x1-a0*x0-a1*x3-a2*x4` and similarly ell2. All six equations and the actual full/restricted quintic are saved in `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_T37_OPEN_LOCUS_20261107.json`. One bounded whole-cell lex Gröbner calculation completed in 13.684 seconds and gives seven polynomials. This is one whole six-variable incidence, not a selection of planes.

Set `u=a4,v=a5`. The last polynomial is a nonzero constant times

\[
(v+6u^2)\,B(u,v),
\]

where

\[
B=373248u^5+82944u^3v-139968u^2v^4+3888uv^2-31104v^5+1.
\]

This supplies the precise previously omitted nonrational-curve entrance. The factorization of a projection equation alone is not a primary decomposition of the full plane source. In particular no unproved assertion about all exceptional lifts is made.

## Exact generic identity on the lifted B branch

The existing fourth and fifth a3-linear Gröbner rows determine `a3`, and then the first three determine `a0,a1,a2`, on their nonzero-coefficient opens. Over the exact field `QQ(u)[v]/B`, every original Lagrangian equation is zero after substitution. A literal degree-five Euclidean gcd of the binary restriction `D|A=0` has degree two. Its square divides that binary quintic, leaving a linear factor. A three-coefficient linear lift and one additional scalar lift then give ternary L,Q, and **all21 degree-five coefficients verify exactly `D|plane=L*Q²`**. All coefficients and bindings are recorded in `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_T37_DOUBLE_CONIC_20261107.json`.

This exact field calculation took 47.271 seconds, maximum RSS 56,880 KiB. The projective multiplication locus `[LQ²]` is closed, so the identity extends over the entire closure of this actual generic lifted branch. A reduced H plus irreducible quartic is outside that locus. Thus this branch cannot supply the requested E6/A5 normal-cone endpoint.

The native Singular general-factor attempt was not used as a certificate: it rejected a multivariate parameter-dependent `minpoly`; the generated power notation also exposed a native precedence issue. Its report has `all_passed=false` and contains no factor result. The subsequent exact field gcd uses SymPy finite-extension elements directly. Their element division implements modular inversion correctly; a high-level `Domain.exquo` adapter instead attempts polynomial division, so the same degree-five algorithm is implemented explicitly. These are implementation/tool diagnoses, not mathematical parent exclusions.

## Source-chart complements: a genuine finite algebra remains

To avoid extrapolating generic field identities to unrelated exceptional lifts, let c4,c5 be the a3 coefficients of the two saved linear Gröbner rows. Eliminating u from c4 and computing the univariate gcd with c5 gives the nonunit polynomial

\[
r(v)=408146688v^{10}+42768v^5+1.
\]

Thus the two generic linear charts do not by themselves cover the whole actual lifted source. The exact patch calculation is in `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_T37_CURVE_PATCH_COVER_20261107.json`; it took 0.155 seconds. In the complement,

\[
u=v^3(1889568v^5+90),
\]

and the quadratic Gröbner row, using that v is a unit modulo r, becomes

\[
q(a_3,v)=a_3^2-(349920v^7+15v^2)a_3
       +288684v^9+\tfrac{99}{2}v^4.
\]

The first three Gröbner rows give polynomial values of a0,a1,a2 in a3,v with rational constant denominators. The 20-dimensional algebra

\[
R=\mathbb Q[v,a_3]/(r/408146688,q)
\]

uses basis `(1,v,...,v9,a3,a3*v,...,a3*v9)`. Every original one of the six Lagrangian equations and every saved Gröbner equation reduces to zero in R. The ideal of all those extra relations has rank zero: **the full actual finite source really has dimension20 as a QQ algebra**, not an extrapolated generic curve. This is a finite algebra; twenty distinct geometric points or reducedness have not been separately assumed.

The complete original rank-three plane chart and all base/plane parameter maps are retained. No point or P4 anchor has been chosen in this finite complement.

## The finite support computation and its real limit

A weak possible target exclusion is nonzero strict linear symmetry: both marked E6 smooth-flex and A5 smooth-hyperflex pairs have finite projective automorphism group. For the finite algebra R, all nine linear vector fields are mapped to all21 coefficients of `W(D)`, giving an exact QQ matrix of size420×180. A basis of its kernel whose coefficient ideal is the unit ideal would prove a nonzero strict linear symmetry at every geometric source point. Such a field cannot be Euler, because the coefficient of C5 in D is the nonzero constant1/3. It would therefore exclude these two marked targets, while saying nothing about every other reduced1+4 divisor.

The first implementation hit its 120-second limit after constructing the exact finite algebra; no symmetry result was obtained. A single diagnosed optimization evaluates the original low-degree D coefficients in the finite algebra *before* expanding the degree52 raw substitution polynomials. It completed the following exact stages:

- reconstructed the same full20 algebra and all original Lagrangian relations;
- evaluated all21 actual quintic coefficients as explicit20-vectors;
- formed the entire420×180 QQ matrix.

The optimized run again hit the 120-second limit, during QQ nullspace computation. Its immutable checkpoint is `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_T37_EXCEPTIONAL_ALGEBRA_FIRST_20261107.json`, phase `EXACT_FINITE_LINEAR_MATRIX_READY`, with `all_passed=false`. The original slower checkpoint is kept separately. No unit identity, symmetry, divisor factorization or target exclusion is claimed on this finite complement. The route stops rather than increasing the time bound or trying new anchors. The saved small algebra, all D coefficient vectors and matrix construction are sufficient to resume by a justified exact linear-algebra method.

## Actual-family implication if a remaining finite plane is relevant

There is a direct geometric bridge worth retaining, rather than assuming divisor coincidence. If an actual isolated reduced Lagrangian plane of a known parent has this target marked divisor and its six Lagrangian equations have invertible six-by-six plane derivative, move the plane by `M_e=M0+e*B`, scale both normal coordinates by e and scale the tensor by e. The only possible negative term is `e^-1{ell1,ell2}|M_e`, whose constant coefficient vanishes by the Lagrangian condition. The resulting full tensor is regular. Its magnetic endpoint is the full six-dimensional derivative image. For e nonzero it is an actual coordinate transform and scalar multiple of the fixed known parent. Thus full Saito-frame/marked normalization plus that derivative would yield an actual all-magnetic containment, rather than an inverse boundary inference. These hypotheses have **not** been verified on the finite complement here.

The present result excludes the generic lifted B branch and its closure from this particular uniform normal-cone route to E6/A5; the explicit finite complement remains unclassified. Other Schubert/source charts, moving periods, nonuniform coordinate contractions and unrelated parents remain outside the exclusion. The global T39/T40 graph theorem and their separately verified E6/A5 tests are independent results.

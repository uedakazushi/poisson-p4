# Generic five-parameter S61 input and normal42 circuit: independent review

This review checks the four frozen packets INPUT_PREFLIGHT,
CONSTANT_UNIT_PREFIX, SMALL_TAIL_FIELD, and PLANE_RANK3_CIRCUIT, all dated
20261104_v63. Its independent code reconstructs the historical weight-one
columns from the immutable S61 model, passes to the affine first jet, removes
the Lee part, and verifies the primitive 24-by-16 basis and adapter determinant
36/625. It substitutes the sparse five-parameter source into every term of all
134 original equations. It also reconstructs the universal Q51/Q52, center,
L51-by20, and N45-by20 specializations. All 150 reconstruction weights and the
70/64 distribution of equation weights 4/5 are checked.

The separate C++ checker is an interpreter for saved arithmetic certificates,
not an elimination program: it makes no pivot choices, searches, new points,
or modular substitutions. The Python review emits exact QQ rational-function
identities in H,L,s0,q0,t0. It checks the recorded 11 constant pivots, seven
field pivots and final zero block, all 51-by2 kernel equations, free-coordinate
identity, both old/new kernel spans, all N*T plane entries, and the shared
denominator. For the rank-three certificate it checks every row scaling and
every recorded operation, including all multipliers and normalized rows.
Nonzero numerator/denominator tests are exact QQ tests. No large 3-by3
determinant or adjugate is expanded after substitution.

Thus L has rank 18 and its kernel is the displayed two-dimensional space over
K=QQ(H,L,s0,q0,t0). The original plane matrix W has rank 3: its selected rows
3,6,13 correspond to u6,u9,u19, its column permutation is [2,1,0], and its
nonzero selected determinant is minus the product of the three row factors
and three saved pivots. Row-primitive matrices are used only to prove rank;
the coordinate map uses the original W.

An independent commutative-polynomial evaluator checks all 1608 circuit nodes
against the abstract entries of W, the center and 45 new coordinates. It
checks both the determinant recipe and adjugate recipe, and proves the
universal B adj(B)=det(B) I identity without substituting large polynomials.
Consequently E=W adj(B)/det(B) has selected block I3, and

    q = center + E z + (42 nonpivot coordinate directions)

is an affine coordinate isomorphism over K. Its inverse first reads the three
selected coordinates and then subtracts E z from the remaining 42. Every
original q coordinate is retained. All inverted quantities depend only on
the five parameters and are nonzero in K.

The source with lambda stripped is not asserted to dominate W16 under Sp4
alone. The accepted sparse6 dominance certificate, together with physical
Gm covariance, supplies the generic cover. For a term of equation weight d,
F(lambda*u,lambda^2*q)=lambda^d F(u,q); the full 150-coordinate reconstruction
has the corresponding physical weights. This is a generic coefficient-field
input, not a chart valid at every parameter value. No substitution into the
42 normal variables, elimination, finiteness, unit ideal, or exhaustion of
generic S61 solutions is established here.

The review report records current input/evidence hashes and the exact native
identity counts. Compiled executable files are kept only in ignored work/.

# A5 smooth hyperflex: one full rank-four normalization and the two global graphs

2026-11-07. PRIMARY E. Exact producer computation, independent replay pending. No new component is counted.

The complete μ5 source batch gives the canonical full-source orbit representative `R=h²`, with transformed quartic

\[
F=h^2y^2-2hw^2y-16hy^3+w^4.
\]

The full six-magnetic source has this dense projective orbit of dimension 23, with the zero magnetic-invariant boundary in its closure. Its full ambient derivative has rank 100 and its actual orbit rank is 48, so the component dimension remains between 23 and 25. No higher local calculation is made here.

Choose the structural rank-four point `(h,x,y,z,w)=(1,0,0,0,1)`. The quartic value is `F=1`, directly from its pure `w4` coefficient; this is not a point search. The actual translation `w_new=w_old-h` moves it to the affine origin. The constant Pfaffian is `25/36`. Exact symplectic Gram–Schmidt gives an invertible rational GL4 frame. The complete PGL matrix has determinant `36/25`, with its exact inverse recorded.

The actual perspective transformation and the saved S61 adapter reconstruct all 150 coefficients, verify all 134 equations and verify the full projective Jacobi equations. All linear and cubic affine terms are retained. The 45 graph residuals are saved for each graph:

- T39/Q51 has 11 nonzero residuals; the first, at zero-based free-coordinate index 2, is `6275/15552`.
- T40/Q52 has 12 nonzero residuals; the first at that index is `1175/1296`.

The computation took 2.317 seconds and maximum RSS 105,628 KiB. It is the only graph normalization of this canonical A5 orbit; no E6 anchor was repeated.

`S61_GLOBAL_COMPONENT_EVALUATION_FIBER.md` supplies the additional global theorem: connected Sp4 makes the full J-evaluation fiber of each global component irreducible; the four perspective pivots give a global polynomial Ga4 cross-section; the two closed A16 universal graphs therefore equal the full normalized intersections of T39/T40. Subject to independent acceptance of that theorem and the present arithmetic, these nonzero residuals exclude the dense A5 canonical orbit from both full global components, independently of the point/frame used for normalization. The zero magnetic-invariant source boundary can still meet either component. No other parent is identified or excluded here.

The report `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_A5_HYPERFLEX_GRAPH_20261107.json` contains exact producer bytes, source bindings, the full PGL/inverse/perspective matrices, all 61 coordinates and both full residual lists. Later isolated modes preserve the graph function AST; use the report's execution bytes for exact historical replay if current whole-file bindings have advanced. The finite job gate is a contract check, not a mathematical certificate.

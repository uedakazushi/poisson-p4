# E6 smooth-flex source: one faithful S61 graph entry

Status: exact finite negative diagnostic, not a parent exclusion. No component or membership is asserted here.

The complete marked-source batch and its magnetic normalization identify this full six-magnetic source with the closure of one projective orbit of dimension 23. Its full ambient derivative has rank 100, so the tangent upper bound is 25. This makes the two known 25-dimensional S61 graphs T39/T40 legitimate recognition targets, without implying containment in either.

The canonical source has transformed quartic

\[
F=w^3(h-8y+3w)-y^4+6y^2w^2,
\qquad R=y^2.
\]

The structural evaluation point is `(h,x,y,z,w)=(1,0,1,0,0)`: `F=-1`, so the affine tensor has rank four there. The actual translation `y_new=y_old-h` moves this point to the affine origin. Exact symplectic Gram–Schmidt on its constant value gives a rational invertible GL4 frame. The combined PGL matrix has determinant `-376/75`; its inverse and the constant-value identity with the standard symplectic matrix are stored in the report. The original constant Pfaffian is `-75/376`. No point search was used.

The existing S61 implementation then applies the actual perspective transformation. All 150 quotient coefficients reconstruct, all 134 defining equations vanish, and the full projective Jacobi equations vanish. Thus the comparison is of the full affine `J+L+Q` tensor through actual invertible coordinate changes.

The exact 45 graph equations have 15 nonzero residuals in each graph. The first nonzero residual, at zero-based free-coordinate index 2, is `-2999375/1484448` for T39 and `-496525/989632` for T40. Every residual is retained in the report.

This proves only that this faithfully normalized graph entry fails both graph equations. It does not exclude a different point/frame/PGL translate or a different parent component. The E6 smooth-flex source remains unassigned. This finite route stops here; no further anchor or higher-jet calculation is implied.

Reproduction:

```
python3 scripts/classification.py gate research/classification/jobs/plane_hyperplane_quartic_e6_flex_graph.json
python3 code/classification_ga9_cuspidal_quartic.py --e6-flex-graph --output NEW_REPORT.json
```

The report `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_E6_FLEX_GRAPH_20261107.json` contains the exact producer bytes, all coordinate matrices, the faithful S61 point and both full residual lists. Later isolated runner modes preserve every old function AST; the report provides the exact execution bytes if the runner has been extended. The job gate checks bindings and is not a mathematical certificate.

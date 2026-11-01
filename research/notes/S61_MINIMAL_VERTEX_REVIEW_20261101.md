# Independent review of reduced minimal-vertex exhaustion

Date: 2026-11-01. This review applies to the corrected minimal-vertex exhaustion with the four additional baseline free-coordinate lines described below. It does not accept the earlier eighteen-curve list as complete.

The reduced case split is accepted after that correction, subject to its explicitly stated historical fixed-trace and vertex-normalization dependencies. The conclusion concerns unsaturated infinity over a sufficiently general beta. It does not remove these points by homogenizing-variable saturation or establish an affine-fiber bound.

## 1. A substantive correction: four canceled baseline pivots

After g2=g4=g5=g6=g7=g8=0, the only surviving original restrictions are rows70,84,87. Row87 first determines g3. The actual coefficients of g0 and g1 in rows84 and70 then are

\[
\begin{array}{c|cc}
&g_0&g_1\\\hline
M_0&-\dfrac{2(3k+1)(9k+7)}{15(2k+1)}&
\dfrac{2(3k+2)(9k+2)}{15(2k+1)}\\[2mm]
M_1&\dfrac{(9k+4)(9k+7)}{30(2k+1)}&
-\dfrac{(9k+2)(9k+5)}{30(2k+1)}.
\end{array}
\]

The first draft inspected only the denominators of the canceled function-field solution. This missed a factor in each coefficient. At the four additional values, both the coefficient and its constant term vanish identically in the nuisance parameters:

| Plane | k | g3 | Free coordinate |
|---|---|---|---|
| M0 | -1/3 | -2 | g0 |
| M0 | -2/3 | 2 | g1 |
| M1 | -4/9 | 1 | g0 |
| M1 | -5/9 | -2 | g1 |

The other linear coordinate is determined, giving four genuine affine lines over base divisors. They cannot be discarded by a codimension-two exceptional-locus argument. The independent audit identified this omission and the main result was corrected before acceptance.

At k=-1/2,-7/9,-2/9 the full specialized baseline ideals are still generically unit ideals over QQ(K,L,U,V). The review recomputes these and the four free-line bases. Away from all five values for each plane, the three explicit pivots give uniform uniqueness for every nuisance value. This avoids the invalid inference that polynomial-looking function-field solutions by themselves rule out nuisance vertical components.

## 2. Complete algebraic case coverage

The highest-cubic radical is the union of the four coordinate axes: the five displayed original restrictions force the square of every pairwise product into their ideal. Therefore either all highest coefficients vanish or exactly one of g5,...,g8 is nonzero. The reverse reduced containment is immediate on each axis.

When all four vanish, the exact row combination forces g2*g4=0. The baseline is treated in §1; for g2 or g4 nonzero, the resultant has precisely the two reported constant k roots for each plane and axis. Its content is a nonzero rational number, so it hides no nuisance divisor.

For the g5 and g8 axes, the three equations used to force the other quadratic coordinate to zero genuinely have coefficients in QQ[k,g3,g_other], with no omitted nuisance parameter. Their saturated unit-ideal checks are therefore uniform. The remaining resultants have only the displayed k roots.

For g6 and g7, the review checks the actual determinants of both successive linear solves, not merely denominators after cancellation. In both planes the first determinant is

\[
-\frac3{25}g_i^2(3k-1)(3k+4).
\]

The second is respectively

\[
-\frac6{25}g_6^2(2k+1)(3k+4),\qquad
-\frac6{25}g_7^2(2k+1)(3k-1).
\]

Together with the forcing-row exception k=-3 or2, these are exactly the stated exceptional values. There is no hidden nuisance factor in either solve. After the solves, the two residual polynomial ideals are checked in both directions over QQ, with all nuisance parameters still polynomial variables. Thus their factorization into the rational-tail case and the isolated H/D-axis case holds uniformly on the specified localization.

All exceptional-axis unit ideals are recomputed over the stated nuisance function field. Every rational graph is recomputed from the independently reconstructed restrictions and compared by both ideal containments with the saved graph. The graph computations use the same efficient elimination order and CAS as the main calculation; the independent parts are the restriction reconstruction, uniform pivot analysis and case audit. All full restrictions vanish on its formulas, including at zero of the free curve parameter. The rank-one case is similarly checked by both ideal containments over QQ(U,V).

## 3. Independent coefficient calculation

The review reconstructs M0 and M1 directly from their displayed Hamiltonians, applies S61 normalization, substitutes the rank-two and rank-one beta charts, and extracts all terms of degree two in the quadratic coordinates. It does not import the main restriction helper or the main case-splitting code. All 134 entries of both restrictions agree with the report. The archived S61 reconstruction and equations remain an explicit shared dependency.

The original G6 identities are read as part of the main certificate: their reduced vanishing follows from the highest-axis radical and from the rank-one bases. This review does not independently reimplement the full universal G6 construction.

## 4. Intrinsic exceptional-locus argument

The vertex incidence I={(beta,[v]):tau_beta(v)=0} has dimension seventeen. Its rank-two induced-quadratic open has a smooth normalized presentation Sp4×A11→I of relative dimension four. The three normalizations of the nondegenerate binary quadratic remove three dimensions from the parabolic stabilizer of the line; the remaining four-dimensional stabilizer accounts for that relative dimension. Discrete root ordering does not affect it.

For rank one with nonzero additional scalar, the previously verified shears remove K,L and the remaining torus fixes the scalar. The normalized beta slice has dimension eight, and Sp4×A8 presents the rank-one divisor of I, of dimension sixteen, with relative dimension two. These are presentations of the whole indicated opens, not merely dominant maps of the correct dimension. The explicit normalizing operations are smooth in characteristic zero after the finite root choices.

After retaining the four free baseline lines, every remaining omitted rank-two case lies in a fixed-k divisor together with a proper additional nuisance condition, or in H=D=0. The latter has codimension two for every k: H is a nonzero linear form in K,V and D is a nonzero linear form in L,U, and neither coefficient pair ever vanishes simultaneously. Function-field unit ideals show that the exceptional set inside each fixed-k divisor is proper; clearing finitely many denominators makes it closed after enlarging it. A similar proper condition inside the rank-one divisor again has codimension at least two in I.

The correct object is the intrinsic bad beta-vertex locus: a pair admitting a solution outside the symplectic saturations of all the listed charts. Every normalized presentation of such a pair is bad. Its full inverse image is therefore contained in the exceptional union. Smooth relative dimension four, or two on the rank-one divisor, implies that the intrinsic bad locus has dimension at most fifteen. Projecting to W16 and taking closure preserves that upper bound. This argument does not treat a choice of frame as intrinsic and does not depend on dimensions of the omitted q fibers.

The earlier exclusions of zero induced quadratic and rank-one zero scalar have the same explicit dimension-at-most-fifteen property. Removing the resulting finitely many proper image closures gives the claimed generic-base coverage.

## 5. Scope and reproduction

The review accepts point coverage using the full archived minimal fixed-trace normal forms, the independently reviewed vertex necessity/normalization results, and the corrected finite case split. It makes no assertion about the number of global components from the number of chart families. Exceptional points and intersections on the rational curves still require local saturation analysis.

`code/s61_minimal_vertex_review_20261101.py --output <new-path>` produces `reports/S61_MINIMAL_VERTEX_REVIEW_20261101.json`. The report records the uniform baseline pivots, all four released-coordinate specializations, saturated graph comparisons, exact two-step axis determinants, residual ideal equalities and exceptional unit-ideal checks. The frozen report must hash the corrected main files.

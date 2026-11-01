# Reduced minimal-nilpotent infinity over a sufficiently general first jet

Date: 2026-11-01. Exact case exhaustion in characteristic zero. Local saturation at the intersections and exceptional points of the listed curves is not determined here.

## Theorem and scope

For sufficiently general beta, every nonzero minimal-nilpotent modular point of the135-leading zero set is represented by one of the finite charts below, up to symplectic change of coordinates and projective scaling of Q. The same reduced list applies to the134-leading zero set within the two minimal planes M0,M1. The other open minimal-nilpotent stratum is already excluded generically by `S61_OTHER_NILPOTENT_INFINITY_EXCLUSION_20261101`.

The list consists of the generic elliptic-curve chart, rational curves over resonance divisors of the vertex incidence, two additional isolated-axis charts, and the rank-one induced-quadratic chart. The counts of normal-form charts are not counts of global irreducible components. Symplectic permutations identify some charts, and their closures can intersect.

Combined with the separately reviewed exclusions of nonnilpotent, regular-nilpotent and nil22 modular types and the R0/C0/L0 analysis, this supplies the remaining reduced normal forms at135 unsaturated infinity over a sufficiently general beta. It does not prove that this infinity is entirely removed by saturation, or that the full affine fiber is finite or has only the five known points.

## 1. Reduction to the vertex incidence

The all-beta necessity theorem `S61_MINIMAL_INFINITY_VERTEX_NECESSITY_20261101` shows that the modular image line [v] lies on E_tau. The historical minimal fixed-trace classification leaves precisely the full nine-parameter planes

Q_0=x² dx wedge dy+dy wedge(-G_w dz+G_z dw),
Q_1=(x²/2)(dx wedge dy+dz wedge dw)-xz dy wedge dz
    +dy wedge(-G_w dz+G_z dw),

with G in Sym³(x,z,w) modulo x³ and coefficient order

(g0,...,g8)=(x²z,x²w,xz²,xzw,xw²,z³,z²w,zw²,w³).

Use the vertex normalizations already proved, including their generic-base incidence exclusions, in `S61_R0_GENERIC_EXHAUSTION_20261101`. The zero induced quadratic and the rank-one case with zero additional scalar have beta images of dimensions at most14 and15, respectively, so neither occurs for a sufficiently general beta.

For rank two, use

u4=u12=u5=u13=0, u14=k, u49=2k+1,
u26=U, u35=V, u50=K-5V, u51=5U-L.

All six omitted beta coordinates remain free. Put

H=(3k+1)K-(6k-2)V,
D=(3k+2)L-(6k+8)U.

## 2. The four highest-cubic axes

In both planes, constant multiples of original rows74,75,100,101,102 give

g5g6, g5g7, g6g8, g7g8, 9g5g8+g6g7.

The radical is the ideal of all pairwise products among g5,g6,g7,g8. The code checks that the square of each pairwise product belongs to the displayed ideal; the reverse reduced inclusion is immediate. Thus at most one of g5,g6,g7,g8 can be nonzero. This is a uniform statement, before any genericity in k,K,L,U,V.

The exact restriction of G6 to this entire vertex chart is

G6=9(2k+1)g6g7

in both planes. It therefore vanishes on the reduced support of the original134 restrictions. G6 can still change the nonreduced structure.

## 3. Vanishing highest cubic and the generic chart

When g5=g6=g7=g8=0, the combination g2*row71+g4*row86 equals (6/5)g2g4. Hence either g2 or g4 is zero.

If both vanish, the only remaining nonzero original rows are70,84,87. Row87 first solves g3. The actual coefficients of g0 in row84 and g1 in row70 after that substitution are

| Plane | g0 pivot | g1 pivot |
|---|---|---|
| M0 | -2(3k+1)(9k+7)/(15(2k+1)) | 2(3k+2)(9k+2)/(15(2k+1)) |
| M1 | (9k+4)(9k+7)/(30(2k+1)) | -(9k+2)(9k+5)/(30(2k+1)) |

Away from the zeros of these pivots and 2k+1, these three linear equations uniquely give the previous three-G solution, uniformly in all nuisance parameters. It is essential to use the pivots before canceling the generic solution: solution denominators alone miss free-coordinate specializations.

At k=-1/2,-7/9,-2/9, the specialized ideal is the unit ideal over QQ(K,L,U,V). Consequently any additional solution there requires a proper further nuisance condition; section7 applies. The other pivot zeros instead give four genuine baseline lines, valid for every nuisance parameter:

| Plane | k | Free coordinate | Other nonzero coordinates |
|---|---|---|---|
| M0 | -1/3 | g0 | g1=2L-6U, g3=-2 |
| M0 | -2/3 | g1 | g0=2K-6V, g3=2 |
| M1 | -4/9 | g0 | g1=-L+2U, g3=1 |
| M1 | -5/9 | g1 | g0=-K+2V, g3=-2 |

All g2,g4,g5,g6,g7,g8 are zero on these lines. Exact specialized Groebner bases and direct substitution certify the formulas. They are additional divisor charts, not codimension-two exceptions. Independent review identified their omission in the earlier draft; the current certificate explicitly enumerates every root of the three triangular pivots.

If g2 or g4 is nonzero, eliminating g3 from the two necessary linear equations gives precisely the k values in the first four rows of the table below. The code records the exact resultants. No unspecified resultant factor is dropped.

## 4. All rational resonance charts

Each row records the only potentially nonzero highest or extra quadratic coefficient t. The remaining eight g coordinates are explicit polynomials in t over QQ(K,L,U,V), stored in the report. A lexicographic Groebner calculation after adjoining h*t-1 gives exactly the graph equations together with h*t-1. Substitution verifies that each graph extends to t=0 and satisfies every original restriction. The t=0 point is already on another listed chart.

| Plane | Free coefficient t | k values |
|---|---|---|
| M0 | g2 | 1/6, -11/12 |
| M0 | g4 | -1/12, -7/6 |
| M1 | g2 | -1/3, -11/12 |
| M1 | g4 | -1/12, -2/3 |
| M0 | g5 | -1 |
| M0 | g8 | 0 |
| M1 | g5 | -1, 0 |
| M1 | g8 | -1, 0 |
| M0 and M1 | g6 | -2/3 |
| M0 and M1 | g7 | -1/3 |

There are eight M0 and ten M1 graph charts in this table, in addition to the four baseline lines in section3. Their images in q coordinates need not be affine lines: some coordinate formulas are quadratic in t. Their projective closures are rational curves, and no disjointness assertion is made.

For g5!=0, original rows force g4=0 universally: saturating the three nuisance-independent necessary equations by g4 gives the unit ideal over QQ. The remaining two linear equations in g3 have exactly the listed resultant roots. For g8!=0 the identical argument forces g2=0. Thus no further base divisor is hidden in either cubic-end case.

## 5. The g6 and g7 axis charts

For g6!=0, row72 forces g4=0 unless k=-3. Away from

k=-3, 1/3, -4/3, -1/2,

two linear rows determine g2,g3, and two further rows determine g0,g1. The certificate verifies that their rational denominators only involve g6 and (2k+1)(3k-1)(3k+4). The two remaining equations, after clearing these units, generate precisely

(k+2/3)H,
(k+2/3)(D*g6-N6),

where N6=2(3k+2) for M0 and N6=3k(3k+1)/2 for M1. Equality of these two polynomial ideals is independently recomputed over QQ, not inferred from common zeros.

The k=-2/3 case is the rational chart in the table. Otherwise H=0 and D*g6=N6 give the isolated-axis chart; the four remaining g coordinates are uniquely reconstructed by the saved rational formulas. Loci with H=D=0, or H=0 and an additional fixed k value, have codimension at least two in the vertex base and are excluded as explained below.

For g7!=0 the symmetric computation has exceptional k values

k=2, 1/3, -4/3, -1/2,

and leaves exactly

(k+1/3)D,
(k+1/3)(H*g7-N7),

with N7=2(3k+1) for M0 and N7=-3(k+1)(3k+2)/2 for M1. The k=-1/3 case gives the table chart; otherwise D=0 and H*g7=N7 give the second isolated-axis chart.

At all four exceptional k values for each axis, the full original axis ideal saturated by the nonzero axis coefficient is the unit ideal over QQ(K,L,U,V). All sixteen calculations (two axes and two planes) are saved. They imply that a remaining solution at such a k requires a further proper nuisance condition; they do not assert uniform emptiness at every nuisance value.

The M0 H=0 component is separately certified dominant, generically finite and locally reduced in `S61_MINIMAL_AXIS_INFINITY_20261101`; its full homogenized rank45 proves its generic-beta points are removed by saturation. That local assertion is not used to deduce the present reduced exhaustion.

## 6. Rank-one induced quadratic

For the remaining vertex stratum, the prior normalization gives

u4=u12=u13=0, u5=-1, u14=1, u49=2,
u26=U, u35=V, u50=-5V, u51=5U.

This uses the nonzero additional scalar and exact symplectic shears; the other six beta coordinates remain unrestricted. Over QQ(U,V), the complete original134 restriction has the following Groebner bases.

For M0, let t=g8. The solution is

```
g0=-2V/3,
g1=2U/3+20V²t/9-10V/27,
g2=g3=g5=g6=g7=0,
g4=1/3-5Vt/2,
g8=t.
```

For M1 the unique solution is

```
g0=4V/9, g1=-4U/9+20V/81,
g2=0, g3=-1/2, g4=-1/12,
g5=g6=g7=g8=0.
```

The exact G6 restriction is -9(5g5g6+6g5g7-4g6²), which reduces to zero in both bases. The M0 chart is a rational curve and the M1 chart a point. The computations are over the indicated function field: exceptional nuisance values are treated by the generic-base argument, not silently asserted identical.

## 7. Why the parameter-field exceptions can be removed for general beta

The incidence I={(beta,[v]):tau_beta(v)=0} is a vector-bundle incidence of dimension17. On its rank-two induced-quadratic open, the normalized chart above together with Sp4 is a smooth surjective presentation of relative dimension4. This is the same normalization used and independently reviewed in the R0 exhaustion. Rank-one induced quadratic is a divisor in I; after the stated normalizations its presentation has relative dimension2.

Every omitted rank-two base case has codimension at least two in I: it is either an intersection of the independent divisors H,D, a specified k value plus a further nonzero nuisance equation, or an intersection of a specified k value with H or D. Exact unit-ideal calculations over a nuisance function field imply emptiness on a nonempty open of that divisor; denominators can be cleared to give a proper additional closed condition. Likewise the function-field exceptional set in the rank-one divisor has codimension at least one within that divisor, hence at least two in I.

This dimension assertion is applied to the intrinsic bad locus, meaning a beta-vertex pair that admits a solution outside the symplectic saturations of the listed charts. Every normalized presentation of such a pair must lie in the exceptional sets of the exhaustive algebraic case split. For a smooth presentation of relative dimension r, the full inverse image of an intrinsic bad locus has dimension dim(bad)+r. Containment in a subset of codimension at least two therefore forces dim(bad)<=15. This avoids treating a single arbitrary choice of frame as intrinsic.

Projection to W16 cannot increase that dimension. Removing the closures of the finitely many bad beta images, together with the previously excluded zero-quadratic and rank-one-zero-scalar loci, gives a nonempty base open on which every point is covered. The argument concerns beta images of vertex-base exceptions, regardless of how large their q fibers might be.

## 8. Reproduction and remaining work

```
.venv/bin/python code/s61_minimal_vertex_exhaustion_20261101.py --output work/minimal_vertex_exhaustion_check.json
```

The report `reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json` stores the complete restrictions, radical certificates for highest terms, all actual baseline pivots and their roots, all four additional baseline lines, all resultants, all18 further rational graphs, generic isolated-axis formulas, all unit-ideal exception calculations, and rank-one Groebner bases. The six unused beta coordinates are never specialized.

This completes the generic reduced point coverage for the minimal planes. It does not determine full homogenized local rings at resonance-curve intersections, at their projective endpoints, or at the rank-one charts. In particular generic-point saturation results on a positive-dimensional family cannot remove its whole generic-beta fiber: exceptional divisors of a seventeen-dimensional total family may still dominate W16. These local strata remain necessary before deriving any affine-finiteness conclusion.

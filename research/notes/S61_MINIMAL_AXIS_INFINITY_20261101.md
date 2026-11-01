# A dominant component of isolated minimal-nilpotent infinity points

Date: 2026-11-01. Exact characteristic-zero construction and smoothness certificate.

There is a sixteen-dimensional irreducible component of the unsaturated leading incidence which dominates W16 and is generically finite over it. Its general points have nonzero minimal-nilpotent modular matrix and belong to M0. Its generic fixed-beta points are reduced and isolated. This is a different family from the elliptic infinity curve in M0; existence of isolated infinity points does not imply an affine extension.

## Formula

In the normalized vertex chart put

u4=u12=u5=u13=0, u14=k, u49=2k+1,
u26=U, u35=V, u50=-(9k+7)V/(3k+1), u51=5U-L.

The other six beta coordinates are free. Let

D=(3k+2)L-(6k+8)U,

and work on (3k+1)(3k+2)D != 0. Define

```
g0=-2V(D-8U)/[(3k+1)D],
g1=2U(L-2U)/D,
g2=-4V(3k+2)/[(3k+1)D],
g3=-8U/D,
g6=2(3k+2)/D,
g4=g5=g7=g8=0.
```

With

G=g0 x²z+g1 x²w+g2 xz²+g3 xzw+g6 z²w,
Q=x² dx wedge dy+dy wedge(-G_w dz+G_z dw),

the certificate verifies all 135 leading equations as rational identities with all ten source beta parameters free. The cubic coefficient g6 is nonzero on the chosen open, so this is not the previously constructed G=a x²z+b x²w+c xzw family.

For context, writing K=5u35+u50, the new source condition is

(3k+1)K-(6k-2)V=0.

The displayed formulas arise from the g6-axis case in the reduced highest binary-cubic support. The main certificate does not depend on a claim that this case split is exhaustive. It directly substitutes its final formulas.

The code also checks the exact restriction G6=9(2k+1)g6g7 on the entire nine-G M0 family after u4=u12=u5=u13=0 and u49=2u14+1. It vanishes on the displayed axis family because g7=0. It is not identically zero on the full nine-G vertex chart.

## Exact point and differential

Take k=1,U=2,V=3,L=5 and set the six remaining beta coordinates u25,u27,u34,u36,u37,u48 to 2,3,5,7,11,13. The sixteen beta values in the standard order are

(0,0,0,0,1,2,2,3,5,3,7,11,13,3,-12,5).

The G coefficients are

(-19/2,-4/3,5,16/3,0,0,-10/3,0,0).

The canonical discriminant is exactly 4267428033666687552, hence nonzero. The source beta parameters and Sp4 action have differential rank sixteen. The joint map to (beta,[v]) also has rank sixteen. All minors are stored.

At this point both the134 and135 leading Jacobians with respect to the forty-five q coordinates have rank44. Since the equations are homogeneous, the radial direction is in the kernel; after projectivization the relative tangent space is zero. The full nine-variable fixed-v Jacobian has rank9 as an additional check.

## Why this is a whole component

The leading incidence is a closed subscheme of W16 times P44, of ambient dimension60. Its relative q-Jacobian rank44 bounds its tangent dimension at the exact point by16. The irreducible image of the explicit source and Sp4 action has dimension at least16 because its beta differential is surjective. Therefore the local dimension and tangent dimension both equal16. The local ring is regular and the image closure is the unique irreducible component through the point.

This component dominates W16 and has the same dimension as W16. Because it is closed in W16 times P44, its projection is proper and generically finite. After restricting to a nonempty base open it is finite etale; no degree is claimed here. This argument establishes actual dominance, rather than extrapolating from the presence of one special infinity point.

The separate full-homogenized local calculation gives relative Jacobian rank45 at this point. The present code independently reconstructs its homogenizer column directly from the original monomials and verifies the saved 45 by 45 determinant

-124102513209883356037120000000000000000/59049.

The selected rows already belong to the original134 equations. Thus the fixed-beta full homogenized local ring is C and t=0, for both134 and135 systems. This conclusion uses the full rank45, not leading rank44 alone.

The same rank45 minor is nonzero on a dense open of the sixteen-dimensional component. At those points the full homogenized total incidence is smooth of dimension16 and the infinity component is its unique local component; hence saturation removes it locally. The complement of that open in the component has dimension at most15. Its image in W16 has proper closed closure (the component is projective over the base). Removing that image shows that the whole generic fiber of this component is absent from the saturated closure. This does not address any other infinity component.

## Reproduction and scope

```
.venv/bin/python code/s61_minimal_axis_infinity_20261101.py --output work/minimal_axis_check.json
```

Frozen report: `reports/S61_MINIMAL_AXIS_INFINITY_20261101.json`.

The generic degree of this infinity component, all other M0/M1 exceptional fibers, and the full affine S61 fiber remain separate questions. No additional projective Poisson component is counted from this unsaturated infinity construction.

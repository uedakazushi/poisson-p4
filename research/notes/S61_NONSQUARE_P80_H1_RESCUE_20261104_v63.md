# Whole-P80 H1 rescue by the cubic products

Date: 2026-11-04. Exact characteristic-zero primary certificate; independent
review pending. This is one bounded fixed pipeline on the whole P80 curve.

Use the source parametrization

```
a=-8(7z+80)/(45z^2), e=2(z+40)/(9z), d=-3z(z+20)/80
```

over

```
C80=QQ[z,1/(z(z+20)(7z+300)(7z+80)(17z+216)(7z+60))].
```

The actual raw H map has rank six everywhere on this curve by the accepted
two-curve cokernel result. The first part of this job constructs its
normalized seven-component dual through the same exact source-unit
operators as on P100. Its physical v3 component is one. All sixteen
actual raw H-linear rows annihilate it by full QQ(z) coefficient checks.
All six original source quadrics are then evaluated. All values are nonzero;
q0, q1 and q4 have certified whole-curve unit values. In particular

```
q1(eta)=41600/(81z^2),   q1(eta)^-1=81z^2/41600.
```

The polynomial coefficients and whole-curve unit tests are saved in
`dual.json`. No point, interpolation, new prime or generic rank solve is
used.

Because this eliminated quadratic value is a unit, the same fixed job
continues to the explicit source-cofactor construction. For
`li=vi-eta_i*v3`, its seven source lifts are obtained from the fixed six-row
unit basis. The identity

```
f*g-eta(f)*eta(g)*v3^2
 = (f-eta(f)*v3)*g+eta(f)*v3*(g-eta(g)*v3)
```

and its q1 version give all 23 retained target monomials as H1 projections
of the 112 cubic products, modulo q1H. The latter is an actual eliminated
H-quadratic source relation of NF117.

The resulting **112-by-23 right-inverse coefficients** are saved explicitly
in `right-inverse.json`. All 23 source identities are checked in all 28
physical quadratic coefficient slots, in both the original and actual
delta7 source frames: 1288 full rational coefficient checks. The five
eliminated H-quadratic relations are rebuilt from the full corrected joint
prime rows and checked separately. Every nonzero right-inverse coefficient
denominator is proved to be a unit in C80, so this is a whole-curve result.

The actual source indexing is `cache_row=1+7*genuine_index+light_position`:
the cubic rows are 1 through 112 and q5H is row0. The saved source change
rows43 through49 is exactly genuine6 replaced by `delta*(C6-C1)`, with
`delta=-64000C_core^2`. Its inverse cofactor transport is checked in the
source identities. No optional q5H or quartic cache row is queried.

Thus the H1 projection is surjective at every point on P80, despite the raw
H-linear rank being six. The 117 known quadratic source pivots, the separate
H2 unit and these 23 H1 directions give pointwise rank at least 141 for the
full quartic source matrix along P80. Appended source combinations exhibit
the unit block; no single fixed original 141-minor or new free source
complement is asserted.

All identities concern the H1 projection. Actual cubics have pure tails;
the same full-source combinations retain those tails in the lazy pure98
operator. No `H*li=0` or full `H*v3^2=0` ideal identity is inferred. All 21
quartics remain present, and full cubic57 on the whole curve and quartic239
are not proved by this packet.

Producer: `code/s61_nonsquare_P80_H1_rescue_20261104_v63.py` with the two
new fixed native kernels `...P80_H1_dual_kernel...cpp` and
`...P80_H1_right_inverse_kernel...cpp`.
Evidence: `reports/S61_NONSQUARE_P80_H1_RESCUE_20261104_v63`.
The combined dual-plus-right-inverse pipeline took 6.690 seconds under the
single 57-second/2-GB/32-MiB/150,000-term guard. It stored 17,982 polynomial
terms with native peak length31. There was no full134 Gaussian, POWER
expansion, pure98 or quartic coefficient query, GB, new point or prime.
Certificate SHA-256:
`aad1bab0292abe873ee0b3b1fd9770f9c4f75b79150999e9e44d68b781e25b2a`.

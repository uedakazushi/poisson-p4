# Sharp raw-H rank drop on two explicit rational curves

Date: 2026-11-04. Exact characteristic-zero primary certificate; independent
review pending. This concerns the actual raw H-linear coefficient map and
its heavy block, not the full cubic or quartic matrix.

Let T be the actual FULLU Theta ring, identified in the frozen coordinate
packet with

```
T=QQ[z,d,d^-1,1/(z*(z+20)*L*K*Fnum*Jnum)].
```

The permitted named-unit stripping of D14 left the polynomial

```
Q6=-9z^4-360z^3-3600z^2-540dz^2-10800dz-8000d^2.
```

Define `Pm=3z(z+20)+m*d`. The new fixed coefficient checks give

```
Q6=-P80*P100,
P100-P80=20d,
1=(P100-P80)/(20d).
```

Thus the two ideals are comaximal on T. Neither is declared to be a unit;
their whole zero schemes are retained.

## The whole two curves are nonempty rational opens

On `Pm=0`, `d=-3z(z+20)/m`. Since d is already a unit, z and z+20 are
units in either graph ring. Direct unary QQ coefficient checks give

| Localizer | m=80 | m=100 |
| --- | --- | --- |
| L | `-z(7z+300)/8` | `-z(z+60)/2` |
| K | `-(z+20)(7z+80)/8` | `-(z+20)^2/2` |
| Fnum | `-z(z+20)^2(17z+216)/16` | `-z(z+20)^2(203z+2700)/250` |
| Jnum | `-z(z+20)(7z+60)/8` | `-z(z+20)(z+12)/2` |

Their products are nonzero polynomials, so the full localized curves are
the nonempty rational open schemes

```
T/(P80)=QQ[z,1/(z(z+20)(7z+300)(7z+80)(17z+216)(7z+60))],
T/(P100)=QQ[z,1/(z(z+20)(z+60)(203z+2700)(z+12))].
```

No point specialization was used to establish nonemptiness. For an explicit
source-frame description the original coordinates and FULLU localizers are
also checked by rational-polynomial cross multiplication:

| Source quantity | m=80 | m=100 |
| --- | --- | --- |
| a | `-8(7z+80)/(45z^2)` | `-10(z+20)/(9z^2)` |
| e | `2(z+40)/(9z)` | `(z+100)/(9z)` |
| d | `-3z(z+20)/80` | `-3z(z+20)/100` |
| B | `z(7z+300)/60` | `z(z+60)/12` |
| C | `-4(7z+300)/(9z)` | `-25(z+60)/(9z)` |
| F | `-10(17z+216)/(9z)` | `-(203z+2700)/(9z)` |
| J | `(7z+60)/30` | `(z+12)/6` |

Thus no part of the claimed curves has escaped the original FULLU
localization. The two graph ideals are prime and comaximal, so their union
has the exact reduced quotient `T/(P80*P100)`, not just a set-theoretic
description.

## Actual coefficient-map cokernel

The orientation here is the row-generator presentation
`coker(rawH)=T^7/(the sixteen raw H rows)`. The typed finite operators are
saved in `operator-recipe.json`.

First, the previously independently reviewed original raw16-by-7 reduction
is pinned, rather than assuming a pointwise rank lower bound. Its constant
GL7 matrix and both inverse products are checked. Its literal unit formulas
are rechecked:

```
delta=-64000C^2,
eta=192aFB^2C^3/125,
gamma=-4FBC^2/625,
cZ*delta_z+cW*delta_w=eta.
```

Replacing source row 6 by `delta*(C6-C1)`, applying the pair GL2 of
determinant `1/gamma`, clearing source rows by the saved eta-functional, and
removing its unit row/target preserve the cokernel. The second fixed pivot
is `-BC/14062500` after the saved known-unit row clearings. All initial and
subsequent row-unit specifications are checked to involve only nonzero QQ
and FULLU units. The accepted peer receipt supplies the full coefficient
replay of these operations. Thus the original raw16-by-7 map and the actual
checkpoint14-by-5 map have the same cokernel, including after Theta
base change.

Second, the exact Theta raw8 unit pivot has inverse `37500C^2*g1/B^4`.
Its middle three entries are literally zero. The known B^4 target scaling
leaves the seven-by-four matrix `[B|S]`, with the other six source rows zero.

Third, the known-unit stripping certificate proves the selected B minor
Delta is a unit on all T. Use `adj(B_pivot)/Delta` to normalize those three
source rows and a target shear to clear their S components. The four
remaining scalar entries are the four fixed defects divided by Delta.
Three are identically zero. For raw11 this uses domain cancellation by the
nonzero, explicitly computed V, not a global V inverse. The fourth is
`D14/Delta`, a known-unit multiple of `P80*P100` by the exact source and
localizer quotient certificates.

Consequently the actual cokernel is

```
coker(original raw16x7)
 = coker(heavy14x5)
 = T/(P80*P100)
 = T/(P80) x T/(P100).
```

This is an exact module statement with the typed unit operations; it is
stronger than a generic rank calculation. It gives, at every geometric
point on Theta,

```
heavy14x5 rank = 4 on the two curves, 5 off them;
original raw16x7 rank = 6 on the two curves, 7 off them.
```

In particular a uniform raw H rank-seven assertion on FULLU is false.
Pure cubic components can compensate for this heavy-coefficient defect.
Neither full cubic57 failure nor quartic239 failure is inferred. The full
Poisson normal algebra, generic-beta exhaustion and any new component
claim remain outside this packet.

Producer: `code/s61_nonsquare_Theta_rank_drop_curves_20261104_v63.py`.
Evidence: `reports/S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63`.
The fixed job took 0.046 seconds under 57-second/2-GB/32-MiB guards, with no
GB, new point, prime, solver or minor search. Source dependencies and the
original raw-H independent review are pinned in the report.
Certificate SHA-256:
`612faf65f2795c2b4ff4da64387b16496ddc6595f4682433d79027c11ea530fb`.
Operator-recipe SHA-256:
`0189663af6bf5d389b10dd5e4981d07004ba3b67f34dffb37b45d065fe016b40`.

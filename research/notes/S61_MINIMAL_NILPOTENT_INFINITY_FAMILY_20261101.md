# Two minimal-nilpotent curves at unsaturated infinity

Date: 2026-11-01. Field: characteristic zero. Status: exact universal construction and generic fixed-vertex calculation; no complete minimal-nilpotent exhaustion or saturation assertion.

## Statement

For sufficiently general primitive first jet beta, the zero set of the 135 leading forms has two distinct curves whose dense open points have a nonzero minimal-nilpotent modular matrix. Their normalizations are the elliptic null-correlation curve E_tau associated with beta. The first has generic quadratic bivector rank two; the second has rank four. This corrects any proposed exclusion of all nonzero-modular infinity over a general first jet.

Only the unsaturated leading-equation scheme is considered. Neither curve is asserted to be reached from the affine S61 fiber. The full generic affine fiber remains unclassified.

## Exact formulas

Write J=dx wedge dy+dz wedge dw, v=dy and

G=a x^2 z+b x^2 w+c xzw.

The two tensors are

Q_0=x^2 dx wedge dy+dy wedge(-G_w dz+G_z dw),

Q_1=(x^2/2)(dx wedge dy+dz wedge dw)-xz dy wedge dz
    +dy wedge(-G_w dz+G_z dw).

Both have modular matrix A with sole nonzero entry A_(y,x)=-2. It has rank one and A^2=0. For Q_0 the Pfaffian is zero; for Q_1 it is x^4/4. Thus their projective families are distinct.

Use the same dense normalized vertex chart as the R0 theorem:

u4=u12=u5=u13=0, u14=k, u49=2k+1,
K=5u35+u50, L=5u26-u51.

The other six beta coordinates are unrestricted. On
(2k+1)(9k+7)(9k+2) != 0, set

- Q_0: c=-2/[3(2k+1)],
  a=-[2K+18(2k+1)u35]/[3(2k+1)(9k+7)],
  b=-[2L-18(2k+1)u26]/[3(2k+1)(9k+2)].
- Q_1: c=-(3k+1)/[3(2k+1)],
  a=[-(3k+1)K+12(2k+1)u35]/[3(2k+1)(9k+7)],
  b=[(3k+2)L-12(2k+1)u26]/[3(2k+1)(9k+2)].

The certificate substitutes all sixteen beta coordinates and checks all 135 leading forms over QQ. For the displayed three-parameter G, only equations 70,84,87 survive before imposing these formulas. G6 is identically zero on this shape after imposing the normalized vertex conditions, without using the displayed a,b,c formulas. It is not asserted zero for arbitrary beta.

The computation also retains all nine coefficients of G and computes the Groebner basis of the first 134 restrictions over QQ(k,K,L,u26,u35). For each family the basis consists of the nine linear equations setting G to the displayed unique solution. This is a function-field uniqueness result; it does not classify the exceptional values of the chart parameters.

## Geometry and dominance

The vertex condition tau_beta(v)=0 defines E_tau in P(V). The prior R0 generic-exhaustion theorem establishes that the displayed beta chart and Sp4 cover a dense open of this incidence: the induced quadratic on v-perp/v has rank two generically and can be diagonalized to the cross term. The two constructions therefore transport to a dense open of the incidence. Stability of the 135 leading zero set under Sp4 follows from the reviewed universal G6 transformation argument; that argument is not restricted to flag tensors.

Uniqueness of the generic fixed-v solution gives descent independent of the choice of normalized symplectic coordinates. Conversely, the line v is recovered from a nonzero rank-one modular matrix as its image. Consequently, for a sufficiently general beta, each construction is birational to E_tau. A rational map from the smooth projective curve E_tau to projective Q45 extends across its finitely many undefined points. Its image is a closed irreducible curve with normalization E_tau. Boundary points may have zero modular matrix; no disjointness of the two closures is asserted.

The code independently recomputes the following exact differential checks at
k=1,u26=2,u35=3,u50=4,u51=5, all other unrestricted beta coordinates zero:

- The eleven beta parameters plus ten Sp4 directions have rank 16, with determinant 48.
- The joint map to (beta,[v]) has rank 17, with determinant 20.
- The fixed-v nine-variable leading-equation Jacobian has rank 9 for each family.

The joint rank and recovery of v prove a genuine positive-dimensional fixed-beta family. These assertions do not rely on counting numerical samples.

The Q_0 sample has (a,b,c)=(-25/18,98/99,-2/9); Q_1 has (2/9,-47/99,-4/9). The exact full q-Jacobian ranks and minors are stored separately in the report. They describe leading-equation tangent spaces only, not completed local rings or the homogenizer.

## Reproduction and dependencies

Run:

```
.venv/bin/python code/s61_minimal_nilpotent_infinity_family_20261101.py --output work/minimal_infinity_family_check.json
```

The report is `reports/S61_MINIMAL_NILPOTENT_INFINITY_FAMILY_20261101.json`; its root entry is Q_0 and `M1_family` is Q_1. Both store full S61 coordinates at the exact samples, all rational formulas, the generic fixed-vertex Groebner basis, rank witnesses and input hashes.

Geometric dependencies are `S61_R0_GENERIC_EXHAUSTION_20261101.md` for the dense vertex chart, `NULL_CORRELATION_DUALITY_20261101.md` for E_tau and its smooth elliptic geometry, and `S61_UNIVERSAL_WEIGHT6_REVIEW_20261031.md` for Sp4 stability. The construction itself does not require the historical full Q45 classification. Further minimal-nilpotent strata and exceptional vertex parameters remain to be treated separately.

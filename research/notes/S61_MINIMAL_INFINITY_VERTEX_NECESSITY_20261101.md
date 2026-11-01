# Minimal-nilpotent infinity necessarily lies over E_tau

Date: 2026-11-01. Exact rational, set-theoretic statement for all beta.

Let Q belong to either of the two minimal-nilpotent fixed-trace planes M0,M1, or a symplectic conjugate. If the 134 leading equations vanish at (beta,[Q]), then

tau_beta(v)=0, where [v]=P(image(modular Q)).

This assertion needs no genericity or discriminant assumption. It concerns the unsaturated leading equations. It does not classify the solutions on this vertex locus and makes no assertion about the affine S61 fiber.

## Uniform reduction

Normalize the nonzero minimal modular field to -2x dy. The planes have the full nine-parameter forms

Q_0=x² dx wedge dy+dy wedge(-G_w dz+G_z dw),

Q_1=(x²/2)(dx wedge dy+dz wedge dw)-xz dy wedge dz
    +dy wedge(-G_w dz+G_z dw),

where G runs through all cubics in x,z,w modulo x³, with coefficients
(g0,...,g8) in the order x²z,x²w,xz²,xzw,xw²,z³,z²w,zw²,w³.

At v=dy the null-correlation section is tau_beta(v)=(-u12,u4). The SL2 on (z,w) preserves the modular field and both full planes, and acts transitively on nonzero covectors. If tau were nonzero it could therefore be normalized to u4=1,u12=0. No other beta coordinate is set to a numerical value.

Write the remaining 64 leading equations as M(g)b, with b in the stored sixteen-coordinate order. The certificate constructs constant left combinations annihilating every column except the u4,u12 columns. On the normalized tau chart they give seventeen polynomials solely in g. Their exact rational Groebner bases are

- M0: g2,g3+6/5,g5,g6,g7;
- M1: g2,g3+1/5,g5,g6,g7.

Thus these five equalities hold at every point of the chart, not just at a generic point.

Substituting them into original row86 leaves respectively -14u13/5 and u13/5. Hence u13=0. Original row72 then becomes the nonzero constant -32/25 for M0 and -2/25 for M1. Both charts are empty. This proves the assertion on all points.

## Certificate

`code/s61_minimal_infinity_vertex_necessity_20261101.py` independently constructs both full nine-parameter tensors, normalizes them to the S61 model, and forms the complete beta-linear leading matrices. It saves every constant row combination, verifies all canceled columns identically, recomputes the two Groebner bases and verifies the final contradictory rows.

Run:

```
.venv/bin/python code/s61_minimal_infinity_vertex_necessity_20261101.py --output work/minimal_vertex_check.json
```

Frozen report: `reports/S61_MINIMAL_INFINITY_VERTEX_NECESSITY_20261101.json`.

For sufficiently general beta, E_tau is the smooth elliptic null-correlation curve. The separate `S61_MINIMAL_NILPOTENT_INFINITY_FAMILY_20261101` certificate supplies the unique generic fixed-vertex solution in each plane and hence two elliptic infinity curves. The remaining tasks are the exceptional points over E_tau and the nonnilpotent modular strata; this necessity theorem does not silently discard them.

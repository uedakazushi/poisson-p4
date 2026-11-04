# Source preflight for the ten-dimensional H-kernel on P100

This is a source-only receipt. It creates no new curve coefficient values.
The six accepted pivot sources have genuine indices `[1,6,8,9,11,14]`,
actual135 labels `[20,61,85,86,99,129]`. Their504 pure coefficients and the
original-source-frame seven linear lifts sigma are already evaluated.

The ten nonpivot indices are `[0,2,3,4,5,7,10,12,13,15]`. Literal pure
sources are0,3,7 (actual18,39,68). The remaining seven actual source labels
are23,40,43,92,115,120,134. Their588 pure coefficients are already present
as QQ[a,e,d] polynomials in the complete normalized57-by91 source matrix.
No lazy-oracle root, ancestor or coefficient query is needed. All112 full
H-linear coefficients agree exactly with the accepted raw16-by7 source.

Write H_i for its H-linear form and eta for the accepted dual with eta_W=1,
W=v3. The saved sigma lifts obey `sigma_j*H=v_j-eta_j*W`, with support only
on the six pivot sources. For each nonpivot i define

```
alpha_i,p = sum_j H_i[j]*sigma_j,p,
K_i = C_i-sum_p alpha_i,p*C_p.
```

The H-part is zero because H_i.eta=0. Restricting sigma to the six non-W
physical columns gives a square left inverse to the selected pivot H block.
Thus those six H rows are a basis of ker eta, and the ten K_i, with their
nonpivot identity block, form the whole **formal source kernel** of the
16-to-seven H map. The pivot-first source permutation has sign+1, and the
remaining block triangular source change has determinant1. This does not
assert independence of their polynomial images modulo quadratic products.

Sigma is in the original16 source frame. To use the actual delta7 frame,
transform arbitrary source cofactors by `new6=old6/delta` and
`new1=old1+old6`, with `delta=-64000*C^2` a whole-curve unit. The whole
source-basis determinant is then1/delta, also a known unit.

The direct missing588 source polynomials contain76421 core terms; all16
pure source rows contain139746 terms. Therefore a later job must stream
one nonpivot source rather than retain the full source JSON alongside its
curve values. The ten pure84 output bounds total39305 RF terms, at most5250
for one source. Cached pivots and sigma have6131 and341 RF terms. The
source-only pass took1.295 seconds.

Code and report stem:
`S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_PREFLIGHT_20261104_v63`.
Its source input saves all full polynomials and exact source recipes for
inspection; that large JSON is not intended as a live arithmetic input.
Report SHA-256:
`f9c4d7c24a868293ca066c0b8da4d8fb88c05cd5de72ab0b4f8ceb870ad8d626`.

No70 quartic products,15 mixed Koszul vectors, four quadratic-overlap lifts,
full89-vector basis, rank or full239 conclusion is supplied here. Those
require their own actual source coefficients and identities.

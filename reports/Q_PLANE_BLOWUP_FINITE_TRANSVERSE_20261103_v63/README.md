# Actual blowup removes generic directions of both15-planes

Full actual35-variable model, exactQQ.  In the z22 blowup chart write
x=r*(v(p14)+n20), and G_i=F_i/r^2.  At the explicit rational parameter
point shared by the two frozen15-planes, the full350x21 linear matrix
[dQ_v,C(v)] has rank19.  The pivot includesr, and its normalized2D kernel
has radial coordinate0.  Hence the kernel-restricted quadratic forms are
exactlyQ_i(kernel), and after actual implicit correction their coefficient
matrix has rank3.  All19/3 minor values and full matrices are retained.

Thus the transverse closed fiber has length3.  The completed UNSATURATED
weak transform is finite over QQ[[delta_p14]].  Its r-saturation quotient
remains finite, withr a nonzerodivisor if nonzero, so the completed actual
exceptional fiber has dimension at most13.  An entire15-plane would give
a14D projective chart through that point, a contradiction.

Therefore EACH15-plane meets the actual ordinary tangent cone in a PROPER
closed subset, of affine-cone dimension at most14.  Generic directions
are removed.  This does not exclude all points, claim the selected point
absent, bound the full model dimension14/12, or construct/count branches.
Puiseux possibilities are handled by saturation/dimension, not dismissed
using failure of an ordinary second-order arc lift.

Proof: research/notes/Q_PLANE_BLOWUP_FINITE_TRANSVERSE_20261103_v63.md.
Producer: code/q_plane_blowup_finite_transverse_20261103_v63.py, ~2seconds.
The diagnostics folder separately preserves the unsuccessful polynomial
kernel ansatz methods; the successful proof does not depend on them.

```sh
python code/q_plane_blowup_finite_transverse_20261103_v63.py \
  --output work/q-plane-blowup-v63/replay
```

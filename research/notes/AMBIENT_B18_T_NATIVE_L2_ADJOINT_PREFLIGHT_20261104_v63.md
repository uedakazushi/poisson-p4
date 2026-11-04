# NativeL2 source adjoint after the full lower correction

This source-only packet uses existing actual origins and first coefficients.
It computes no L2,N2, new source jet, or quartic value.

The full corrected source multiplier is
c1=-Xp*a1-Xstar*b1-XL*l1. All84 normalized unit-star columns and37 selected
linear columns have separate old/native constant part zero. Thus S0 times
their native source combination is exactly zero. The Xp*a2 term starts in
degree2. Consequently the entire native L2 sensitivity in degree3 is
L2*(S1*c0-S0*Xp*a1), independent of the unknown lower b1,l1,a2 values.

For each actual619 constant source, the proper first residual is recovered
by applying the saved invertible C4 matrix(N0^t)^-1 to
nativeXc1-nativeXp0*a1. The worker checks the forward identity and folds
the exact Gamma4 source adjoint against L2 times this residual. It records
both the actual48-by48 L2 coefficient adjoint and the sub-adjoint of a free
proper N2 direction using origin C4 multiplication.

A nonzero symbolic adjoint is a source dependency, not a proof that a new
N2 computation is necessary: the oldB scalar source functions may depend on
the same gauge and cancel it coherently. Other Res4/basis contributions to
L2 are not discarded. No unknown coefficient is assigned a value here, and
actualHF4 or totalT length does not follow from this packet.

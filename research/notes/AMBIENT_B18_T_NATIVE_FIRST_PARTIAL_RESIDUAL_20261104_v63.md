# The actual first source identity does not split into old and native parts

The accepted full56 identity Xc1=Xp0*a1 is an equality after summing all
original768 source columns. It does not imply the corresponding equality
on the old0..639 and native640..767 column subsets separately.

This diagnostic uses only the actual619 c0/a1 arrays, the frozen first27
Schur source combinations, full56-by768 X0, and the independently accepted
old640/native128 first coefficients. For every619 source it independently
forms the two partial residuals

    sum_partial X1_j*c0_j-sum_partial Xp0_p*a1_p.

Their sum is zero in every56-by8 slot. The native residual is nonzero for
82 of the619 sources, with634 exact QQ coefficients. At source619-index288
its entries are41:-1,49:1/4,257:-1. The slot order is8*actual fullV56 row
plus the external variable index, ordered zA0,...,zA6,rA. The old partial
residual is the exact negative, including at this witness.

In the actual scalar Res4(C4) frame the proper residual is
Omega1=L0^(-1)*native_partial_residual. The accepted L0 is a unit, so these
82 proper residuals do not vanish. This excludes the proposed automatic
cancellation from the full source identity alone. It does not establish a
nonzero final L2 sensitivity: Gamma4*(L2*Omega1) may still cancel, and that
source-adjoint fold remains the next gate. No nativeN2 coefficient or new
physical source value is computed here.

# Nine actual fourth-order source functionals: independent review

This review verifies the nine requested source coefficients without solving
for any source coefficient. The primary is a new restricted higher-order
calculation, not a cache-only calculation.

The common105 r-line coefficients through order3 are reconstructed directly
from the original bivariate joint-source certificate by setting both plane
parameters to T. They agree with the separately transported second-order
cache. Substitution into the original105 source equations, with every unknown
fourth coefficient zero, gives the stated fourth residual F4. A separate
sparse forward derivative calculation reconstructs the entire original
Jacobian J0. All nine supplied adjoints satisfy rho J0=e_target, so the
requested coefficient is exactly -rho F4. This is the coefficient equation
J0 u4+F4=0 and does not require constructing the full vector u4. The CA
target additionally agrees with the independently reviewed ownCA order6
readback.

The selected targets are also checked at their source: all619 fourth-order
sensitivity folds are reconstructed with separate physical and honest
regular-function derivative labels. The copiedCB57 inverse is its actual
zero7 section psi(0^7,YB7,YB8). Its derivative is transported through the
two actual YB functions, using the saved own inverse first coefficients;
it is not an additional independent inverse variable. The accepted nonlinear
source2/copy verification supplies the literal copy and unit-uniqueness
bridge. Every surviving fourth-order derivative is a coordinate of the
original105 source, and their union is precisely the nine requested targets.

No source solve, mixed third-order source or gauge calculation, complete619
quartic value, actualHF4, or exact totalT length is inferred here.

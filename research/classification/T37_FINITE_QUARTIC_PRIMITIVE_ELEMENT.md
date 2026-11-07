# A single primitive element for the T37 quartic field

This is an additive hand-algebra supplement to the frozen `T37_FINITE_QUARTIC_FIELD_DESCENT.md`. No earlier file or certificate is changed and no computation is launched.

From

\[
3U^2+11U+9=0,\qquad Z^2-(10U+15)Z+33U/4+99/2=0
\]

one obtains

\[
U=\frac{N(Z)}{D(Z)},\quad
N=4Z^2-60Z+198,\quad D=40Z-33.
\]

The denominator is nonzero in K: D=0 would make Z rational, contrary to the degree-four field presentation already proved. Equivalently substitution of Z=33/40 in the second relation leaves a nonzero rational constant.

Direct expansion gives

\[
3N^2+11ND+9D^2
=48Z^4+320Z^3+2100Z^2+13860Z+55539.
\]

Consequently K=Q(Z). Since `[K:Q]=4`, this polynomial is irreducible and its monic version is

\[
f(Z)=Z^4+\frac{20}3Z^3+\frac{175}4Z^2
          +\frac{1155}4Z+\frac{18513}{16}.
\]

The whole actual descended T37 parent and plane from the preceding note can therefore be expressed over this single constant number field. This is suitable for exact native algebraic-field arithmetic, including Singular's supported univariate constant `minpoly`. It is not the unsupported multivariate parameter-dependent extension attempted for the generic curve.

As arithmetic checks, the trace of Z from the quadratic tower is `10*(-11/3)+30=-20/3`, agreeing with f's cubic coefficient, and

\[
N_{K/Q}(Z)=N_{Q(U)/Q}\big((33/4)(U+6)\big)
 =(33/4)^2\big(3-22+36\big)=18513/16,
\]

agreeing with the constant coefficient. The use of this field does not select a numerical embedding or a new P4 anchor, and it does not by itself identify the marked divisor or establish the plane-Jacobian determinant.

# Actual linear G2 source profile from existing finite action matrices

Only the five active variables and ell are needed. For each of the42
chosen degree2 source lifts, the saved fixed67 source-basis receipts
give coefficients alpha in the56 d2 source columns with
d2(0)*E_j=d2(0)*E_selected*alpha. The actual source inclusion at the
origin is I2=E_j-E_selected*alpha, hence d2(0)*I2=0.

Differentiating the exact source-unit transport gives

    G2_linear=Pi17*d2_linear*I2(0).

Derivatives of the target projection multiply d2(0)*I2(0)=0. Derivatives
of I2 multiply d2(0), killed by Pi17. Thus all coefficients are known
from the constant unit matrices and the literal periodic/Koszul linear
terms; no physical source coefficient is requested. The resulting
102-by42 QQ matrix is saved with fixed-minor/column-factorization
receipts and exact kernel vectors if needed.

In the full tensor resolution, the active degree2 block contributes
this rank. The51 active-degree1/zero-degree1 constants are separated by
their distinct zero variables in the active17 target; the9 zero-wedge2
constants are then separated by the ordinary zero-variable Koszul
linear differential. Therefore the full102-source linear rank is the
active rank plus60, and injectivity is equivalent to active rank42.

This only profiles the known resolution. It computes no X1 value,
source implicit equation, T2 origin relation, finite length or H0.

# Actual constant comparison X and the first-jet minimal-generator gate

This packet evaluates only the full source X:768→56 at the closed Rbar point.
The first640 columns retain all chain B48 and honest P8 components. The last128
are the actual Res4(N^t)*S columns, including the sixteen free-zero source
syzygies. Neither B0 nor those syzygies are assumed zero.

At this point the native coefficient algebra is exactly
QQ[tauA,tauB]/(tauA²,tauB²), of length4. Individual normalizer sources need only
their saved first directions: higher powers of their own tau vanish. The CA
source gauge is the saved tauA times its odd center, with no fresh implicit
solve. Its action on the CB frame retains the mixed tauA*tauB coefficient.
All flags, fixed19 readbacks, R17, first5 graph elimination, R12 and the native
H12-to-G8 adapter are evaluated by their finite original formulas. Polynomial
Hadamards use finite telescoping, and all inverses have nonzero rational
augmentations. In particular d and any vanishing source function are never
inverted. The resulting actual N^t*S is restricted in the original basis
(1,tauA,tauB,tauA*tauB).

The fullV56 unit comparison supplies m0:56→24, with kernel dimension32.
All m0*X0 columns are checked to vanish. Projection to the32 unselected source
coordinates is an isomorphism on ker(m0); its fixed column elimination gives
rank r of X0 and explicit constant combinations for every remaining cycle.
It follows already that mu_Rbar(T) is at least32-r. If r≤23, this would give
length(T)≥9 and mu_D(K)=length(T)-4≥5, hence nonzero DVR torsion because the
separately proved generic rank of K is4. The report records whether this
actual constant-only criterion succeeds; it does not assume it.

There is a precise smaller next gate if this bound is insufficient. Use the
constant pivot combinations to replace each remaining cycle X_j by X'_j with
X'_j(0)=0. Then

    (m*X'_j)_linear = m0*(X'_j)_linear.

All unknown first derivatives of m cancel. In each of the three literal
Koszul sectors, d2's linear part is injective on its28 constant wedge
coefficients: the21 wedges among seven linear coordinates are independent,
and each of the seven wedges with the last r^n generator has its own linear
last-component entry. Thus the84 constants w0 are uniquely determined by
d2_linear*w0=m0*X'_1. No V second coefficient or m first derivative is needed.
The resulting complete formula is

    mu_Rbar(T)=116-r-rank(w0).

This formula still requires source-bound X1 for all retained columns. The
native source adapter's first coefficient array has not been serialized in
the prior packets; its exact length36 restriction is a separate coefficient
gate. The present job neither evaluates X1 nor builds a348-ring matrix.
No origin depth conclusion is made unless its constant criterion actually
succeeds and is independently reviewed.

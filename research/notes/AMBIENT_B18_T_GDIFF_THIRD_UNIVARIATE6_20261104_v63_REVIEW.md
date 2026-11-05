# Independent restricted own-source order6 verification

This review postverifies the114 new order6 coefficients from the authorized
ownCA57 and ownCB57 calculations. It performs no new unknown solve. The
primary is a NEW higher physical computation, not a cache-only query.

For both families the checker pins the same original57 unknowns, equations,
signed numerator/denominator nodes, source DAG and unit Jacobian as the prior
independently reviewed order5 calculation. Every provided coefficient0..5
matches that previous calculation. It substitutes the provided order6
coefficients and checks all57 equations through degree6, all physical x9
readbacks, the saved residual with order6 set to zero, and J0*u6=-residual.
The fixed J0 inverse is checked on both sides. These equations and the fixed
lower coefficients prove uniqueness without resolving the114 unknowns.

For each of the eight honest regular quotients, the checker reads the
literal original numerator and denominator. Since h0=0 and g4=0, ordinary
power-series coefficients satisfy b6=g3*h3+g5*h1. The scalar g3 is nonzero;
there is no source-parameter division by the vanishing conductor. All seven
coefficients of b=g*h are checked. The full original CB ratio times f is
retained by the literal source-DAG evaluation.

The newCA target curve is exactly r=s. The oldCB target has linear slopes
(17/12,7/12), zero quadratic coefficient and the independently checked cubic
bend YB3. Thus its actual common third coefficient is the own-slope third
coefficient plus the saved linear h gradient applied to YB3; the quadratic
h term has no bending contribution because YB2=0.

The actual common values are

    oldCB=(23809/1728,0,0,1093/256),
    newCA=(-761/128,0,0,-171/128),
    oldCB-newCA=(68165/3456,0,0,1435/256).

No common105 sixth coefficient, mixed third, gauge3, complete619 quartic
value or HF4 is computed or inferred by this review.

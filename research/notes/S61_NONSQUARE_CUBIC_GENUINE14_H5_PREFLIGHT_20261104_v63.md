# Actual fourteen-source heavy block and its known-unit preconditioner

The heavy five-column block of the actual final `14x21` presentation
is literally the corresponding block of the old eight-step checkpoint.
All 70 scalar roots are traced to those checkpoint entries without
visiting a POWER-circuit node. The remaining original source labels
are

`[10,12,13,14,15,16,17,18,19,20,21,22,23,24]`.

The heavy target labels are `[27,28,29,30,31]`, at final target
positions `[16,17,18,19,20]`. They are four projected heavy
coordinates and psi. Every recent quadratic pivot source has zero
coefficients there, and every associated target tail operation is
identity on these five columns. No POWER coefficient of degree207
or full `14x21` row is expanded in this packet.

The 70 old entries have 43 nonzeros and 13,854 numerator terms.
Their exact known-unit row/column preconditioning is

`Hfinal = diag(rho14) * Ncompact * diag(1,1,1,1,gamma)`,

where `Ncompact` is precisely the independently reviewed A0 raw
H-linear checkpoint after its single completed pivot, and

`gamma = -4*F*B*C²/625`.

All 70 identities are verified, including the zeros. The compact
block has 2,407 terms, 43 nonzeros and at most 148 terms per entry;
its column term counts are `415,304,302,302,1084`. The fourteen
row units are all explicit nonzero rational multiples of `C`,
`a*F*C²`, or `F*C²`. No unknown p or determinant factor is inverted.
Exact rational constants and exponents are saved in the certificate.

## Why the L-coordinate corrections cancel

Let `Xi` be the actual raw oldX H-coefficient row and let
`s=X1_x1=F*B²*C/12500000`. The L-frame relation has unit coefficient
`L_x1=-40C`. Before the old genuine pivot, source row i has the known
clearing unit `lambda_i`. Eliminating the L coordinate using source1
leaves, for each other oldX coordinate h,

`lambda_i*(s*Xi_h-Xi_x1*X1_h)/s`.

The coefficients of L in the inverse coordinate change cancel
exactly. For psi, write

`li = hcore*Xi_Z - a*g1*g2*Xi_W`.

The same expression holds with `Xi_h` replaced by `li`, with one
additional known column factor gamma. The earlier raw H-linear
checkpoint performed only named-unit initial row cleanings and
one completed pivot. If the corresponding cleanup product is
`tau_i`, its row is this cross-product numerator divided by tau_i.
Consequently `rho_i=lambda_i*tau_i/s`. The code derives these
products from the fixed saved transcript and verifies all entries;
it does not rerun a unit scan or Gaussian.

This also provides a legitimate full-presentation preconditioner:
rescale each full source row by `rho_i^-1` and the psi target by
`gamma^-1`, retaining every pure coefficient lazily. The present
packet verifies the heavy block and saves the operator coefficients;
it does not expand or discard those pure tails.

## Parameter scope and next algebraic question

All row and column scalings are units on FULLU. The already accepted
raw H-linear generic-rank theorem therefore gives generic rank five
for this actual heavy block. The accepted explicit chart is
`FULLU intersect D(p*Q104)`, where p and Q104 are the previously
saved polynomials. Neither is declared a unit on all FULLU.
Failure of that one chart does not imply an actual rank drop.

The compact block has three literal zero-heavy rows, labelled
`10,13,16`; these are pure source rows. The accepted sparse pair
has raw checkpoint rows `[2,4]`, columns `[0,5]`, and determinant
`-(6/5)*B*p²`. Its middle three columns are zero. The other three
column directions are supported on the accepted seven-by-three
B block. On D(p), rank failure of the five-column block reduces
to the maximal-minor ideal of B. On p=0 other rows may still
restore the rank, so that reduction is not asserted globally.

The next B-block support/content/minor-cost preflight is owned by
the other agent. No B content stripping, 35-minor cost computation,
actual minor, GB or further pivot was performed here. The compact
equivalence avoids repeating the already-stalled individual-unit
Gaussian. Any next Bezout combination must keep p, Q104 and other
unknown factors as nonunits unless separately proved invertible.

Global rank five of the heavy block and global raw H7 generation
remain open and may fail. Quartic rows can restore directions even
where this heavy-only map drops rank. No full cubic57, quartic239,
kernel count or exhaustion claim follows from this preflight.

## Receipt and preserved implementation failure

The successful fixed verification took 0.837 seconds. The first
attempt stopped before the 70 identities because a polynomial was
compared to a fraction-field element without explicit coercion.
Its source, log and receipt are preserved under `failed-attempt/`.
Replacing only that comparison by its correct field coercion
completed the authorized job. No mathematical method or search
scope changed; combined execution remained below two seconds.

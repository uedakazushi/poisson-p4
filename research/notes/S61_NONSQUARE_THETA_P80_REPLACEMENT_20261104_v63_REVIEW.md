# Independent audit: cubic57 at every permitted P80 point

The original actual P80 source17-by-24 and its selected determinant are
reused through their independent accepted source/coefficient review. This
new checker imports no producer or solver code. It uses dense univariate
integer arrays and exact convolution to verify all408 row/column unit
transports, every specified row and column unit power, and both A/B blocks
in their literal source order.

All six saved old-minor power/quotient identities are independently
multiplied back. The removed powers are242,0,34,27,0,0; the remaining
known powers are54,2,7,26,9,0. Their exact product gives the saved expected
degree144 determinant after the certified-unit preconditioning. Its
relation to the original degree447 determinant is checked in every
coefficient. Thus the solver denominator is the actual prepared
determinant times the saved nonzero rational unit.

The checker independently multiplies all17 equations of the saved selected
vector A X=delta b and the complete integer identity

    S*X00+T*R46=C,  C a nonzero integer.

It does not reconstruct the solution or query another minor. It does not
claim to recheck all119 solution cofactors, only the selected vector that
is sufficient for the theorem.

Over QQ(zeta), A is invertible because its actual determinant is nonzero.
Uniqueness and Cramer's identity identify the selected X00 with the actual
column-replacement determinant times the same nonzero QQ unit. Equality
holds as a polynomial before specializing to a root of R46; no inverse
of the old determinant is used at those roots. The integer Bezout identity
shows this replacement determinant is nonzero at every geometric R46 root.
Off R46 the old determinant is a certified unit times R46 and is already
nonzero. These two actual minors therefore cover the entire permitted
P80 open curve. No irreducibility assumption on R46 is needed here.

All source and target operations are invertible by previously certified
P80 chart units. The actual17-by-24 checkpoint has rank17 everywhere on
P80, hence the actual augmented matrix has rank25 everywhere. The accepted
pointwise source-projection criterion then proves full cubic rank57 at
every permitted P80 point. P100 and full quartic239 remain separate.

No new solve, field rank search, point, prime, GB or adjugate expansion is
performed. The saved large Bezout constant is handled exactly, not by
sampling or modular inference.

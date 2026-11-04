# The actual D4 second Betti number is 112

This independent fixed coefficient check uses the accepted exact literal
horseshoe over the regular local coefficient ring
R0=k[[z1,...,z8,v]], with residue field k=QQ (and hence also over C).
The extractor has a separate independent ROOT review. No new source solve,
point, prime, jet or off-diagonal coefficient expansion is performed.

The complete sequence is R0^112 -> R0^44 -> R0^9 -> D4 -> 0. The first 28
second-syzygy columns are the Koszul relations of
P=(z2-psi2(z1),...,z8-psi8(z1),v-psiv(z1)). The last 84 project onto the
accepted N second-syzygy columns S. Each coefficient of S is linear in
z1,...,z8,v and contains no nu. The whole 324-by84 scalar coefficient matrix
has exactly 252 nonzero entries, all +1 or -1.

There is an explicit identity minor, so no rank search is needed. The first
56 S columns are ordered by triples i<j<k in {0,...,7}; select the coefficient
of z_(i+1) in relation row (j,k). For the last 28 columns select the coefficient
of v in relation row (i,j). The checker reads every original coefficient and
verifies that the selected 84-by84 matrix is exactly I84. It also verifies
that the N projection of the actual horseshoe columns is literally this S,
and that the first 28 columns have zero N projection.

Put x=(z2,...,z8,v). Since psi depends only on the parameter z1, the
coefficient of x_j in the linear part of P_i is delta_ij. For the nil Koszul
column labelled i<j, select target row j and coefficient x_i. These selected
28 coefficient rows give exactly I28. The full selected scalar coefficient
minor of the 112 actual first jets is therefore block upper triangular:

    [ I28   *   ]
    [  0   I84  ]

Its determinant is 1, regardless of the uncomputed off-diagonal derivatives.
This is a minor of the coefficient matrix of first jets, not a 112-column
ring minor of the 44-by112 differential. It makes no impossible ring-rank
assertion. The checker separately verifies that all actual d2 entries vanish
at the origin, using the exact c36(0)=psi8(0)=0 inputs.
It also checks origin zero for every entry of d1, so the accepted nine free
generators and the 44 relation generators form a minimal prefix.

If a relation between these 112 columns had a nonzero constant coefficient
vector, reducing the relation modulo the square of the maximal ideal would
contradict this scalar coefficient minor. Hence ker(d2) is contained in
m R0^112. Since the accepted horseshoe says im(d2)=ker(d1), these 112 columns
are a minimal generating set of ker(d1), and beta2_R0(D4)=112. No knowledge
of the next differential is required for this minimality argument.

No honest9 transition, old relation lifts, finite-Z Tor comparison, T0, depth,
or higher Betti number is determined by this result. All selected source
coefficient rows and both identity minors are saved in the dedicated report.

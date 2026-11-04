# Independent review of the fixed Theta coordinate defects

This review binds the three frozen DELTA_LAMBDA_PROFILES, RAW12_DEFECT and
REMAINING_DEFECTS packets to the already independently accepted raw source
and rational Theta scheme isomorphism. Source rows are [7,9,10,11,12,13,14]
and the fixed pivot rows are [7,13,10].

The independent checker rebuilds all nine original B-coordinate polynomial
identities for raw9,12,14 with reversed-variable Fraction arithmetic. It
then clears only known powers of d and uses generic multivariate polynomial
composition in reversed DEGREVLEX contexts to check all coefficients of
Delta, nine lambda entries, six actual Schur entries and the old V. This is
a separate implementation from the producer's termwise cached powers.
All source hashes and the previous five-packet review are bound.

All three raw12 product polynomials and all eight raw9/raw14 product
polynomials are independently multiplied in exact rational Laurent
arithmetic. Raw12 and raw9 defects are zero. Raw14 has 269 nonzero terms.
The old V has eight nonzero terms; the original raw11 identity and integral
Theta domain therefore imply its Delta-coordinate defect is zero, without
inverting V globally.

Delta has 56 nonzero terms. The actual Theta ring is the previously checked
nonempty localization of QQ[z,d], hence a domain. Delta and the nonzero
raw14 defect have nonzero product there. Over its fraction field B has
rank three and the Schur column adds a fourth independent direction. Thus
the augmented block has generic rank four, heavy block rank five, and the
original raw H has generic rank seven after the accepted unit pivots.

Pointwise, the globally defined row combination kappa annihilates B and
has Schur entry D14. At every residue field where B has rank three and
D14 is nonzero, the augmented block has rank four. This uses no Delta
inverse and bounds the possible heavy defect by the B-defect locus union
V(D14). Neither this locus nor D14 is claimed empty, a unit, or factored.
There is no all-Theta rank drop. No full cubic57/quartic239 conclusion or
negative Poisson classification statement follows.

No points, primes, new minors, rank search, solver or GB are used.
Code: `code/s61_nonsquare_Theta_defects_review_20261104_v63.py` and `.cpp`.
Evidence/report: `S61_NONSQUARE_THETA_FIXED_DEFECTS_20261104_v63_REVIEW`.

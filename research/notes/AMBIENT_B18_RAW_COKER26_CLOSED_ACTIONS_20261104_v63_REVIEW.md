# Independent review: canonical closed quotients and raw26 action tables

2026-11-04. This checker imports no producer and runs no source solve, jet,
CAS, field inverse, new point, prime or rank search. It uses independent
sparse-vector Fraction arithmetic and univariate multiplication modulo r³.

The exact input is the separately independently accepted actual H20 ring
in CA-F9 coordinates (zA7,r_CA), with s_CA=r_CA only after ell=0. The two
canonical maps send zA to0 or -c*r² and use the SAME actual r_CA on both
targets. Each map is recomputed directly on every monomial of the fixed20
basis, not copied from its supplied matrix. The common lifts (1,r,r²)
at indices(0,8,16) prove both maps surjective.

The actual unit-normalized conductor generators show the maps factor
through their own conductor quotients. The beta3 terms, and every omitted
term of ordinary order at least3, vanish after zA has weight2 and r weight1
in QQ[r]/r³. Both final source cubics vanish there as well. The accepted
actual conductor quotients each have length3; hence these surjections are
isomorphisms, with no change of the common parameter. IA*IB is contained
in both conductors, yielding the canonical H20 quotient maps.

All16 generator intertwining identities are checked on every basis column.
An independent sparse-action computation reconstructs the20 basis-scalar
multiplication operators and verifies all800 quotient multiplicativity
identities. The supplied eight IA generators annihilate the A target and
the eight normalized IB generators annihilate the B target. The checker
then reconstructs every block of each26-dimensional coordinate action,
verifies their pairwise commutativity and all64 IA*IB product annihilations,
and verifies the20 supplied basis-scalar actions separately.

The accepted raw normalizer-image formula, followed by the right-exact
quotient by ell, gives the module decomposition (3)+(3)+(20). This is not
a new ring structure on a cokernel, and no kernel is computed by reduction
modulo ell. The callback for an arbitrary scalar is linear in its ALREADY
EVALUATED actual H20 coefficient vector; its action follows from the20
basis-scalar identities. No raw physical rho label is substituted for a
CA-F9 coordinate. Values of delta,nu,a,b, generic DVR actions, actual comparison
cokernel/kernel, torsion, socle and D0 depth remain outside this review.

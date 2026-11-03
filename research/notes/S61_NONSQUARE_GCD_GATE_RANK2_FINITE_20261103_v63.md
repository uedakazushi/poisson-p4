# The rank-two branch of the B,W gcd failure is finite

Use the raw four quadrics and the exact A,B,W=x5 factorization of the separate
structural packet. Work on the known core a*P*F*B0*C*J!=0. Write barQi for the
restriction to B=W=0. Then the locus where barQ0 has rank at most2 and
barQ0,barQ3 have a common nonconstant factor is zero-dimensional in the core.
This includes proportional quadrics of rank2. Proportional barQ3/barQ0 of
rank3 or4 remain a separate unhandled case. This result proves neither the
whole gcd gate nor original rank24/41.

The proof is exact QQ algebra, with no GB or finite-field calculation.

In the verified split coordinates,
4Q0=U R−5 V T, with radical coordinates r=x0 and W. On W=0 the B coefficient
of r is −20a*DB, where

DB=20P²+3Pe+20P+180a.

If DB!=0, eliminating r leaves a nondegenerate rank4 quadric. Thus restricted
rank at most2 requires DB=0. B remains nonzero for P!=0, since its x1,x3
coefficients use Sc,Vc with 3Sc−7Vc=20P. Restriction of a rank4 quadric by one
hyperplane has rank at least2. On DB=0, rank2 requires the remaining B covector
to be isotropic for the dual split quadric. The producer computes that norm
exactly. Eliminating a=−P*(20P+3e+20)/180 gives

norm = −4P²*N(P,e)/10125,

where N is the saved quartic and its e^4 coefficient is the CONSTANT72442188.
Consequently N is monic in e up to a nonzero QQ scalar, and this branch is
contained in DB=N=0.

A rank2 quadric over C is a product of two linear forms. If it shares a factor
with another quadric, that other quadric has rank at most2, including the zero
case. Hence every3-by3 minor of the restricted Q3 matrix must vanish. The B=W
coordinate chart uses only P denominators. Multiplying Q3 by P^4 and substituting
DB=0 gives a symmetric4-by4 matrix H3 with entries verified to be polynomials
in QQ[P,e]. Let M(P,e) be its principal minor with indices0,1,2.

The exact univariate calculation at P=1 produces a saved Bezout identity

u(e)*N(1,e)+v(e)*M(1,e)=1

in QQ[e]. N(1,e) has degree4 and M(1,e) degree18. All coefficients, the global
H3 matrix and both Bezout cofactors are retained and the identity is checked
literally.

Why this specialization proves FINITENESS: if V(N,M) had a curve component,
its common irreducible factor could be chosen monic in e, since N is monic in e
and QQ[P] (or C[P]) is integrally closed. A vertical P-only factor is impossible
because N's leading e coefficient is a nonzero constant. Specializing P=1
therefore preserves that factor's positive e degree, contradicting the saved
Bezout identity. Thus V(N,M) is zero-dimensional over C. The linear equation DB
then determines a uniquely, proving the stated finite core locus. This is a
monic-family argument; a lone sample without monicity would not justify it.

The proof does not invert DB, N, their factors or the minor. It retains their
zeros as case equations. The unhandled high-rank proportional branch, uniform
pair gate and five-quadric rank gate remain open. Independent review is pending.

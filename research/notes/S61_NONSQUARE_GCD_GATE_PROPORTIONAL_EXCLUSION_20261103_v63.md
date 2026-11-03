# Uniform exclusion of proportional quadrics on B=W=0

Use the raw four quadrics, Q2' and A,B,W=x5 of the frozen structural packet.
On a*P*J!=0, the restrictions of Q3 and Q0 to B=W=0 are NEVER proportional.
Together with the separate rank-two finite certificate, this makes the entire
gcd-gate failure locus zero-dimensional. Together with the finite-V Nakayama
certificate, those finite exceptions already have raw rank24. The uniform pair
and five-quadric rank gates remain open, so no whole-core product41 is claimed.

All identities here are exact QQ polynomials; no GB, field promotion, CRT or
unknown-factor inversion is used.

## DB!=0

Use the verified split coordinates(r,U,V,T,R,W), with r=x0 and
4Q0=UR−5VT. Write

B=br*r+bU*U+bV*V+bT*T+bR*R,
br=−20a*DB,  DB=20P²+3Pe+20P+180a.

Eliminating r on DB!=0 leaves the same split Q0 in U,V,T,R. In these coordinates
Q3 has U² coefficient EXACTLY −a³J/12, a known unit. Thus Q3 cannot be a scalar
multiple of Q0, whose U² coefficient is zero. The producer checks this after
clearing the row factor a²*br², so no unexplained denominator is present.

## DB=0

Eliminate a=−P*L/180, L=20P+3e+20. Since a,P are units, L is a known unit.
Here Q0 is independent of the radical coordinate r, and B has no r term.
B is still nonzero, from 3Sc−7Vc=20P in the structural factorization.

Let H be the symmetric matrix of a²Q3 in(r,U,V,T,R). Its rr entry vanishes
identically on DB=0. If Q3|B,W were proportional to Q0|B,W, then its r-linear
part would vanish on the main-coordinate hyperplane B=0. Hence the row
h=(H[r,U],H[r,V],H[r,T],H[r,R]) would be proportional to b=(bU,bV,bT,bR).
All six minors h_i b_j−h_j b_i would be zero.

Divide these minors ONLY by their recorded common powers of P,L,J and nonzero
QQ constants. Five of the resulting normalized polynomials G1,...,G5 have a
saved literal polynomial Bezout identity

sum C_i(P,e)*G_i(P,e)=1.

This alone proves the needed contradiction uniformly. For clarity, its small
elimination structure is as follows. Put delta=−5P+6e−5=−J−9. The first three
normalized polynomials are delta*f1,delta*f2,delta*f3, where

f1=−100P²+120Pe−25P−36e−60,
f2=−50P²+60Pe−25P−18e−15,
f3=−60Pe−40P+72e²+6e−29.

From f1−2f2=25P−30, the simultaneous equations f1=f2=0 force P=6/5,e=13/6,
where f3=118. The exact polynomial identity expressing118 in(f1,f2,f3) is
retained. Thus all three delta*f_i can vanish only if delta=0. On that line,
G4,G5 reduce to q1=25P²+35P−6 and q2=675P²+2445P+1594, with

−(10125P+24822)*q1+(375P+86)*q2=286016.

The producer combines these two elementary identities into the stated literal
five-polynomial Bezout identity, verifies it in QQ[P,e], and saves every
cofactor plus the exact factors stripped from each actual minor. Therefore no
proportionality point exists on the DB=0 core either.

## Consequence and scope

Restricted Q0 has rank at least2. At rank3 or4 it is irreducible, so a common
factor with Q3 would imply proportionality, now excluded. Rank2 common-factor
failures lie in the finite cover handled by
S61_NONSQUARE_GCD_GATE_RANK2_FINITE_20261103_v63; the later finite-V certificate
proves raw rank24 on the entire closed cover. Thus outside that raw-good cover
the gcd premise of the conditional linear-syzygy lemma always holds.

The remaining raw rank-defect locus is therefore contained in the union of the
pair-gate and five-quadric rank-gate failures. Neither is removed here. This is
an exact reduction of the outstanding task, not whole-core exhaustion. Author
replay of the proportionality certificate takes2.97seconds; independent review
is pending.

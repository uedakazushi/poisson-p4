# Five-target T: actual first relations and sufficient quadratic closure

The newly proved84 wedge units change the coefficient requirement. This packet
uses only existing exact X0/X1, their27 constant cancellation combinations,
the full-rank84 wedge constant block and literal d3 first terms. No m1, X2,
new source jet, scalar solve or348 matrix is needed.

First choose a constant32-row frame with the27 pivot columns equal to(I27,0).
The full unit cancellation of those cycles has lower Schur part
b_j-B_p*Y^-1*a_j. After the saved constant cycle combinations, both a_j(0)
and b_j(0) vanish, and B_p(0)=0. Its first part is exactly the constant
five-row projection of X'_j1. Unknown top pivot coefficients do not enter.

Write the remaining relation columns on84+5 targets as(W_j,Q_j), with Q_j(0)=0.
The84 selected W columns form a formal unit matrix. In the minimal five-target
presentation their linear terms are

    R_j1=Q_j1-Q_star1*W_star0^-1*W_j0.

Every derivative of W disappears because Q_star(0)=0. The program implements
the same formula by replaying the frozen exact W0 elimination on Q1; it does
not solve a new84 system or recompute the earlier rank. All657 remaining
cycle relations and their complete linear span are retained.

Before the first27 cancellation, the literal d3 relations have form(d3,0).
Their lower block stays exactly zero during that cancellation. The next84
cancellation makes their lower relation -Q_star*W_star^-1*d3. Since both
Q_star and d3 have zero constant terms, their leading degree-two terms are
exactly

    -Q_star1*W_star0^-1*d3_linear.

No unknown derivative of W, higher cycle coefficient or first pivot lift is
present in this formula. The unit-wedge coordinates of literal d3_linear
are solved in the already fixed normalized84 basis, with its exact signs.

Let rho be the rank of the40×657 actual linear relation matrix. Then the
exact first quotient has length5+(40-rho). A value above8 proves nonzero
DVR torsion in K using its independently accepted generic rank4 and
mu_D(K)=length(T)-4. If rho=37, the first quotient has length8. This alone
does not give total length8.

For that case take the37 selected actual linear relations times all8 variables,
and all168 certified d3 quadratic terms. These464 columns belong to the
degree-two initial relation space in five targets, whose ambient dimension
is5*36=180. If their exact rank is180, m²T is contained in m³T, hence they
are equal. Nakayama applied to the finite module m²T gives m²T=0. The total
length is then exactly5+3=8, and the accepted source-bound depth criterion
gives original cone depth20. This is a sufficient upper certificate: a
partial quadratic span proves no extra lower bound, since the uncomputed
quadratic cycle relations could supply further equations.

The code persists the linear input and rank phase before building only this
certified quadratic family. All source arrays are pinned to independent
origin, first640 and native128 reviews. The scientific outcome, including
whether the sufficient180 condition holds, is recorded in report.json.

# Full minimal Betti vector of the actual same-R0 D4 module

The actual D4 module in the accepted horseshoe has minimal Betti vector

    (9,44,112,182,196,140,64,17,2)

over its regular nine-dimensional coefficient ring R0. This corollary extends
the separately frozen beta2=112 result and does not alter that certificate.
It concerns the same-R0 D4 module, not the full B18 origin ring or its
remaining T0 torsion question.

The accepted exact extension is 0 -> L -> D4 -> N -> 0, where
L=R0/P, P=(x_i-psi_i(Z))_(i=1..8), and R0=k[[Z,x1,...,x8]]. Thus
R0=k[[Z,P1,...,P8]] as well. The eight P_i are part of a regular parameter
system. The nil resolution K is their Koszul complex. The independently
accepted N resolution is minimal, with ranks
(8,36,84,126,126,84,36,9,1). The actual extension's first off-diagonal map
is -c36, and c36(0)=0. The accepted second off-diagonal map also lies in the
maximal ideal m, being a constant linear combination of those residues.

Only the following weak strictness property is needed; a stronger statement
about all filtration levels is unnecessary:

    im(dK_j) intersect m^2 K_(j-1) = m im(dK_j),  1<=j<=8.

For each j-subset J, select in dK_j the coefficient of P_min(J) in target
row J minus its minimum. The selected scalar coefficient matrix is exactly
an identity matrix of size binom(8,j). This is checked in every degree by
the finite checker. If y=dK_j(u) has entries in m², reduction modulo m² and
this identity minor force every constant coefficient of u to vanish. Hence
u belongs to m K_j. Conversely dK_j has entries in m, so dK_j(m K_j)
belongs to m² K_(j-1). This proves the displayed equality exactly, including
coefficients involving the extra parameter Z.

Write a horseshoe differential as [dK_i,h_i;0,dN_i]. Its compatibility
condition is dK_(i-1) h_i = -h_(i-1) dN_i. Suppose h_(i-1) has entries in
m. Since the N resolution is minimal, the right side has entries in m².
It is a Koszul cycle by the preceding compatibility and dN²=0. In positive
degrees the exact Koszul resolution makes it a boundary. In the initial
degree its augmentation vanishes by the actual extension; the accepted
literal h2 already supplies that boundary. Weak strictness then permits a
lift h_i with every entry in m. This induction applies column by column
through the top degree eight, preserving the already constructed h1 and h2.
No value of any higher coefficient or any higher source jet is assumed.

The horseshoe lemma makes these block maps an exact resolution. Above degree
eight both component resolutions vanish, and injectivity at their top terms
also gives injectivity at the top horseshoe term. All entries of every block
map are in m, so the resolution is minimal. Its ranks are the termwise sums
beta_i(N)+binom(8,i), yielding the vector above. The finite certificate also
checks the accepted N rank count from its v-action cone and single unit
cancellation, and the alternating rank sum zero.

This is an existence proof for the higher minimal differentials; their
off-diagonal functions are not serialized. It proves no old honest9
transition or old relation lifts, no finite-Z Tor comparison, and no new
T0/full-origin-depth assertion. No CAS, source solve, jets, points, primes,
or additional residue calculations are used.

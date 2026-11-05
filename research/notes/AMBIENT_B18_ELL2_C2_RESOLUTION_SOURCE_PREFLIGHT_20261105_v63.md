# Exact C2 Loewy data and a complete hypersurface-base resolution

This calculation uses only the nine actual commuting14-dimensional C2
action matrices. It retains exact rank receipts for the Loewy filtration
and introduces no physical coefficient, implicit solve or source jet.

Let A=k[ell]/ell^2, B=A[[zA0,...,zA6,rA]], and L be the actual ell
matrix. The periodic A-free resolution has alternating matrices
ell*I-L and ell*I+L. Its augmentation sends u+ell*w to u+L*w. Direct
coefficient comparison, using L^2=0, proves exactness at every position.
All eight other action matrices commute with L and each other. Tensor
this periodic resolution with the Koszul complex of x_i*I-A_i, using
total differential d_periodic+(-1)^p*d_Koszul on periodic degree p.
Completed translation by the nilpotent commuting A_i proves exactness.
Its first ranks are14,126,518,1302. This is a free B resolution of the
actual C2, not a kernel over a finite conductor quotient.

The fourteen actual C2 basis lifts are rA/ell monomials times the three
physical targets tauA,tauB,tauA*tauB. The three literal target positions
are0,5,10. The polynomial monomial map m0:B14->B3 has m0*eta=I3.
For each of the eleven other basis vectors e_j, the source vector a_j
has the exact identity d1*a_j=e_j-eta*m0(e_j). Its only nonzero source
blocks are the periodic block and the rA Koszul block. The explicit
formula follows from the polynomial identity for r^a-A_r^a, and the
ell block when the monomial contains ell.

The eleven constant source anchors are distinct -1 entries. Their full
polynomial minor is -I plus a nilpotent matrix, and its two-sided inverse
is verified by a finite series. The target basis and source complement
are therefore explicit unit changes. Eleven contractible pairs are
removed, giving first ranks3,115,518,1302. The first3-by115 polynomial
differential is saved; the higher changes are exact finite typed
operators, without expanding518-by1302 scalar matrices.

This closes the C2 resolution/source-unit step only. For an honest T2,
the actual original R0 comparison V3-by56 and X56-by768 must be mapped
to this resolution using source identities d1*m=V, V*n=d1 and
d2*h=I-m*n, with compatible cycle lifts for all768 actual X columns.
No such m/n/h source coefficient is inferred from a finite C2 equality.
Alternatively [eta*V,K11] has target14 and source67; the target split
proves that its kernel is the same actual kernel of V and its cokernel
is C2. This is a method gate, not a value of T2 or H0.

The Loewy calculation bounds scalar evaluation orders on C2. It does
not bound the order lost by source Hadamard division: each division
requires its own numerator coefficient and literal regular-function
identity. Allorder lifts must remain exact operators when unavailable.

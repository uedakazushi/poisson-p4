# Independent finite R-gluing bounds

The bounds length(B_R/J_R)<=16 and
length(D4/image(S_F))<=19 are accepted. Here B_R is the independently
accepted canonical double DVR, J_R is the image of T0 in it, and
S_F=0:_(D0)F^infinity. This is a pure proof and a fixed2-by2 polynomial
calculation, with no repeated R source computation.

The scalar F=v²-nu*z1² is the N projection norm, not the old coordinate
used as a generator of the nilradical. It maps to zero in D4_red=N. The
accepted whole W1-tail nilradical is square-zero; this is an exact completed
ring statement, preserved by the regular tail cuts, not a generic-length
inference. Thus the image of F belongs to that nilradical and F²=0 in D4.
Since T0=ker(D0->D4), F²D0 is contained in T0. Projection to B_R gives
F²B_R contained in J_R.

Write A=k[[U]], B_R=A plus A*epsilon with epsilon²=0, and
F=f(U)+g(U)*epsilon. The fixed polynomial checker verifies that multiplication
by F² in the basis(1,epsilon) is [[f²,0],[2fg,f²]], with determinant f⁴.
The actual projected norm has ord_U(f)=4, so this determinant has order16.
The DVR determinant-length formula gives exactly length(B_R/F²B_R)=16,
and the ideal containment gives the first upper bound. The nil coefficient
g is arbitrary; no physical epsilon embedding is needed.

For the lower reduced-order constraint, T0 maps to zero in the N0 subalgebra
because N0 injects into its w-unit localization. Its reduced R restriction
therefore vanishes on the actual N/R meeting. R is fully torus-fixed and all
positive-weight functions restrict to zero on it. The degree0 meeting ideal
is consequently the accepted contact ideal(U⁴). This proves only
image(J_R->A) subset(U⁴), not equality.

Any finite-colength ideal of the double DVR has Hermite form
(U^m+b(U)*epsilon,U^n*epsilon), with b defined modulo U^n. Indeed choose
a generator of its projection ideal in A and of its intersection with
A*epsilon; subtract the first lift to reduce every other element to the
second. Closure under multiplication by epsilon forces n<=m. Contact gives
m>=4, whereas F² containment gives m<=8. Its colength is m+n. Hence
4<=m<=8, 0<=n<=m and m+n<=16. The checker enumerates the35 necessary pairs;
none is asserted to be the actual pair or automatically realizable with all
the remaining coefficient constraints.

The exact cokernel sequence needs no additional finiteness assumption. The
map D0/S_F -> D4/image(S_F) has kernel(T0+S_F)/S_F=J_R and image
B0/image(S_F), where B0=image(D0). Its cokernel is the already accepted
Z=D4/B0 of length3. Therefore

    0 -> B_R/J_R -> D4/image(S_F) -> Z -> 0

is exact, and the second length bound is16+3=19. It follows that image(S_F)
has finite colength in D4. The kernel T0 intersect S_F is not determined by
this argument. S_F or that kernel need not be finite length or curve-supported,
and their hidden socle and the full-origin depth remain open. No actual m,n,b,
saturation exponent, quotient matrix or nil embedding is computed.

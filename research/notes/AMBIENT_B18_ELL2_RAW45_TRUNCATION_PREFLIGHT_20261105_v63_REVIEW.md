# Tor-aware ell² truncation: independent proof and source-input review

The withdrawn tensor-injection argument is not used. The actual union is
the fiber product of two D-free cubic algebras over
A3=k[r,ell]/(r²,r*ell,ell²). Matching their Weierstrass coordinates gives
the explicit six-element D basis e0,e1,u,v,p,q recorded in the primary.
It retains the finite Tor attachment after reducing moduloell².

The independent checker verifies the universal multiplication matrices
in this basis, with symbolic Weierstrass and beta coefficients. The
important z*v calculation is also valid in this basis: reducing beta*r²
in the B branch gives ell²,ell³,ell⁴ divisibility of its r²,r,constant
coefficients, so rewriting the pair in the six-element basis has an
overall ell². Thus it vanishes in U/ell²U itself, not only branchwise.
The matrices verify r4=ell*r3=r²*zi=zi*zj=0 and r*zi=zi*r. In particular
m4(U/ell²U)=0.

Because U is ell-regular, the sequence from the intersection-product
module into the product ring remains left exact after modelli². This
is the legitimate use of regularity. The intersection-product module is
annihilated by J=IA+IB containingm². Therefore J*m4S=0 andm6S=0. In both
the actual and truncated rings, replacing the two factors f,g by their
cubic truncations changes a product by f*u4+g*t4+t4*u4, which is zero.
Applying this argument in each ring gives the reciprocal ideal inclusions.

The source input is checked independently: beta2=A0^-1*b2 and
beta3=A0^-1*(b3-A1*beta2), including the ordered A1 correction. The plane
quadratics have a nonzero fixed three-by-three minor, and both conductor
cubics have nonzero common-r coefficients. All48 z-dependent cubic terms
factor as M2*z. The unit matrixI+M2 removes them and changes pure-plane
coefficients only from degree4 onward, which the proved truncation allows.

The plane8-generator presentation is closed under every z action. Pair
and gA rows map to zero, and gB rows map to plane multiples of other gB
rows. The degree6 product gA*gB is zero by the already establishedm6
bound. The conditional numerical input therefore has at most147 rows
and74 columns. The accepted DVR decomposition predicts product length33
and raw cokernel length45. This review performs no quotient rank job,
physical action evaluation, new source coefficient, or comparison kernel
calculation, and draws no K-torsion or H0 exponent conclusion.

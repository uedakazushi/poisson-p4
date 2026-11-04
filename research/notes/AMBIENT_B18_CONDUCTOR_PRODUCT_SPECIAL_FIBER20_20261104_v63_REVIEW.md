# Independent review: the actual conductor-product special fiber of length20

2026-11-04. Scope is H/(ell), H=R0/(IA*IB), ell=5(s_CA-r_CA), in the actual
CA F9 coordinates (z7,r,s). No raw physical-rho readback, generic DVR action,
comparison-kernel torsion or D0 depth is asserted.

## Source normalization

No author producer or source worker is imported/run. The independent checker
reads the original saved two57 quadratic sources, the joint CA-plane cubic
restriction, and the CB own-plane cubic. Every conductor row/node is matched
to the original actual source. Both physical maps x_A(YA),x_B(YB) are read
directly from those sources and compared coefficientwise. The saved linear
inverse is checked on both sides; the exact quadratic composition
x_B(YB(YA))=x_A(YA) is replayed. Both off-diagonal linear blocks vanish.

For z_B=A z_A+b(r,s), the checker independently splits every quadratic
monomial according to the specified maximal z-index (ordered-prefix rule).
It verifies A0 and its inverse, every A1 coefficient, b2 and the actual
plane3 source b3. Thus beta2=A0^-1*b2 and
beta3=A0^-1*(b3-A1*A0^-1*b2) are exact through the required degree, including
their z-dependent cubic terms. A2 is not assumed zero: A2*b2 starts degree4
and is unnecessary for products with an IA linear generator modulo m5.

The eighth conductor generator on each side is legally replaced by its
own zero-seven-plane restriction using the exact coordinate Hadamard
identity. Their cubic terms come from their own actual source caches.
The CB tail-plane coordinates have no z-linear terms, so only the already
known linear two-coordinate substitution is required at cubic order.
After the exact ell=0 substitution s=r, the checker obtains
gA3=125*r³/24 and gB3=-15125*r³/144 from those separate sources. No type-only
choice gB=r³ is made. O(m4) errors in either cubic generator times a linear
generator, and O(m4) errors in the normalized z_B generators times z_A,
all lie in m5. Hence the computed 64 products are exactly the actual product
ideal modulo m5. They have 49 quadratic leading terms,14 quartic terms and
one zero product in that truncation.

## Fixed finite certificate

An independent sparse Fraction engine reconstructs the 64 products. Each
of the 579 saved relation rows and every duplicate alias is checked against
its literal product and monomial multiplier, and the complete set of all
eligible multipliers is exhausted. All 475 supplied pivot receipts are
replayed as exact row identities; this is not a fresh rank search. Distinct
unit pivots prove the lower rank, and the checker reduces every relation row
to zero using the given echelon basis, proving the upper rank. All 495
monomial normal forms and the twenty claimed independent basis monomials
are verified.

The eight multiplication tables are reconstructed from those normal forms.
All table pairs commute, and an independent matrix-action test verifies that
each of the 64 actual products kills every one of the twenty basis vectors.
Small20-dimensional image computations give m-power dimensions
(20,19,11,3,1,0), hence HF(1,8,8,2,1).

Finally, the accepted independent conductor-DVR theorem gives the ACTUAL
length of H/(ell) as20. There is always a surjection from this actual ring
to the verified quotient by m5, also length20. The length equality proves
the map is an isomorphism; thus m5=0 and these are the actual special-fiber
tables. This use of equal finite lengths is not a truncated-kernel or Tor
base-change inference. The finite special fiber still does not determine
the generic DVR actions, K torsion, socle of D0, or D0 depth.

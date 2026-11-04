# Actual closed conductor square: exact pair normalization and length348

This packet computes the finite base ring R0/(ell,(IA*IB)²). It constructs
no T or K comparison matrix, specializes no kernel, and makes no torsion or
D0-depth claim. The only new arithmetic is fixed polynomial division, fixed
S-pair verification and sparse finite-ring action checks; no GB search runs.

## All-order normalization of the actual pair

Work only after ell=5(s_CA−r_CA)=0 in the SAME actual CA-F9 coordinates
(zA_0,...,zA_6,r). The actual IA generators are its seven zero coordinates
and an own-plane function with leading (125/24)r³. Thus

    IA=(zA_0,...,zA_6,r³).

The actual first seven IB generators have unit zA-Jacobian, determinant
−1/1024. The accepted quadratic source calculation gives their normalized
quadratic difference c*r², where

    c=(9835/864,0,0,0,29875/3456,−4565/864,0).

Their formal implicit graph is therefore zA=phi(r)=−c*r²+O(r³). Use the actual
eighth IB generator after subtracting its own zero7 Hadamard differences.
Its own-plane initial cubic, pulled through the source-bound linear adapter,
has coefficient −15125/144 on r³. The graph has zA=O(r²), so its higher
corrections cannot change this leading cubic. Hence the eighth generator
restricted to the exact graph is r³ times a unit. It follows that

    IB=(zA−phi(r),r³)=(zA+c*r²,r³)

exactly: all unknown higher graph terms lie in(r³). This is an all-order
implicit-function and ideal argument. It is not inferred from the truncated
length20 multiplication table or from a guessed coordinate normal form.
The two source unit gates and quadratic c are hash-bound to their independently
accepted actual source evidence. No new source coefficient or jet is required.

Since c0 is nonzero, use the explicit rational constant GL7 transformation

    Z=zA_0/c0,  Wi=zA_i−(ci/c0)zA_0  (i=1,...,6).

Its determinant is1/c0 and its inverse is zA_0=c0*Z, zA_i=Wi+ci*Z. It sends
the actual pair to

    A=(Z,W1,...,W6,r³),  B=(Z+r²,W1,...,W6,r³).

The r and ell conventions have not changed. This normalization is asserted
after ell=0 only; raw physical rho coordinates are never substituted.

## Four exact two-variable coefficient ideals

Write A0=(Z,r³), B0=(Z+r²,r³), M=A0+B0=(Z,r²). Expanding A²B² by total
W-degree q gives coefficient ideals

| q | coefficient ideal in QQ[[Z,r]] | quotient length |
|---|---|---:|
|0|A0² B0²|26|
|1|A0 B0 M|14|
|2|M²|6|
|3|M|2|
|at least4|whole ring|0|

The packet checks these fixed lexicographic bases, with Z preceding r:

    q0: Z⁴+2r²Z³+r⁴Z², r³Z³−r⁷Z, r⁵Z²+r⁷Z, r⁸Z, r¹⁰;
    q1: Z³−r⁴Z, r²Z²+r⁴Z, r⁵Z, r⁷;
    q2: Z², r²Z, r⁴;
    q3: Z, r².

Every displayed polynomial has a saved constant cofactor vector in the original
coefficient-ideal generators. Conversely every original generator divides to
zero by its fixed basis. All20 fixed S-pairs divide to zero, with literal
division/remultiplication certificates. Buchberger's criterion therefore
proves equality and the stated standard monomials, without any basis search.

There are1,6,21,56 monomials in six W variables of degrees0,1,2,3. Thus the
full standard monomial basis has

    26+6*14+21*6+56*2=348

elements. The decomposition is an exact coefficient-ideal decomposition, not
just a dimension estimate or a truncated quotient.

## Sparse action tables and source annihilation

The evidence lists every348 basis monomial. Each of the eight canonical
coordinate actions is stored by sparse columns, using the fixed normal forms
in its W-degree block. Dense348×348 matrices are never built. The checker
verifies pairwise commutation and checks that every basis monomial applied to
the unit gives exactly its named basis vector.

All36 pairs of A generators and36 pairs of B generators are included. Every
one of the1296 fourfold source products kills the unit. Commutation and cyclic
generation then prove that these products annihilate the entire348-dimensional
ring. The constant inverse GL7 readback gives and checks the eight actual
CA-F9 coordinate action tables. This is the exact action of

    R0/(ell,(IA*IB)²).

The fixed bases are homogeneous for weights(Z,Wi)=2 and weight(r)=1. All
standard monomials have weight at most9, so the maximal ideal's tenth power
vanishes. The standard monomial r⁹ is nonzero, proving exact nilpotence index10
for this finite base ring. This statement concerns the finite base ring,
not the nilradical of the original intersection or any particular T module.

No actual T presentation, comparison matrix, Tor image, kernel, saturation,
torsion, socle or D0-depth conclusion is calculated. The new input is a finite
source-bound ring and sparse multiplication callbacks for later honest module
presentations. Bounds are60s/2GB/32MiB/150000 terms.

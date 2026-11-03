# Uniform product rank at least37 from an irreducible quadric

2026-11-03. On U={adFBCJ≠0}, J=5ad−6e−4, the ORIGINAL42 quadratic-linear
product rows have rank at least37 at every geometric point. This strengthens
the accepted35 bound; it proves neither41 nor full cubic57/quartic/ambient removal.
No generic3-plane is claimed to work uniformly over the entire core.

The accepted exact six-light coefficient certificate gives a four-dimensional
space W=(q0,q1,q3,q5) modulo ell at every point. In its saved cleared q0,
the symmetric coefficient matrix on coordinates(0,1,2) has determinant
`125a²/256`. q0 contains no original d and needs no a-clear shift in the
P=ad chart. Thus q0 has rank at least3 on U, hence is irreducible over C.

## Lemma: restriction of W to some3-plane is injective

Let V have dimension6 over an algebraically closed field of characteristic0,
and let W⊂Sym²V* have dimension4 and contain a quadric of rank at least3.
In P(W), the rank≤2 locus is a proper closed subset, of dimension at most2.
Consider pairs ([q],L) with L∈Gr(3,V) and q|L=0. The variety of such L for
fixed q has dimension at most6 when rank(q)≤2, and at most4 for rank(q)≥3.
Consequently this incidence has dimension at mostmax(2+6,3+4)=8, whereas
Gr(3,6) has dimension9. Its image is proper, so some L avoids every nonzero
q vanishing on L. Equivalently, W→Sym²L* is injective.

For completeness, the fiber dimensions for ranks1,…,6 are `(6,6,4,4,3,3)`.
If q has rank r, radical dimension m=6−r, and k=dim(L∩rad q), put s=3−k.
The corresponding stratum has dimension

    k(m−k) + [s(r−s)−s(s+1)/2] + s(m−k).

The three terms choose the radical intersection, its isotropic quotient, and
its graph lift. Require0≤k≤m and0≤s≤floor(r/2), then maximize. This also
covers ranks1,2; equivalently their isotropic3-planes lie in a hyperplane
factor of q. The incidence argument is pointwise and requires no fixed L.

## Cubic shadow and lifting

In three variables, four independent quadrics have at least8 independent
linear products: after row reduction their four leading degree2 monomials
are distinct, and their degree3 shadow has at least8 elements. There are
only15 choices of four among the six degree2 monomials; the verifier checks
all of them. Distinct leading monomials give independent products.

Choose L from the lemma, and let I_L have three independent linear generators
ell1,ell2,ell3. The twelve products ellj qi are independent: reducing a
constant relation modulo I_L² gives three separate constant relations among
the independent restrictions qi|L, so all coefficients vanish. These twelve
products lie in the kernel of restriction to L. The image has dimension at
least8. Therefore dim(V*W)≥12+8=20.

Apply this to the four quadrics modulo the original ell. The accepted pair
cover supplies17 independent original cubic products killed by restriction
modulo ell:13 from ell·(Z,W)·S1 and four from ell qi, independent modulo the
first13 because qi remain independent modulo(Z,W). Hence

    rank(original six quadrics × seven light variables) ≥20+17=37

at every point of U. All normalizing row factors used for the saved quadrics
are units on this open. The bound is unaffected by those rescalings.

## Consequence after the concrete33rd pivot

The independently accepted primary24×58 residual has eight remaining product
rows after33 product pivots. Their pointwise rank is therefore at least4.
For their8-row matrix N, the incidence of (p,[lambda]) with lambda N(p)=0
inside U×P7 has dimension at most3+(8−4−1)=6<7. Thus a general constant QQ
row combination generates the unit ideal of QQ[U]. A further scalar unit row
exists globally. No particular lambda or new Bézout vector is certified here;
this is a rigorous existence result guiding the next constructive search.

# A genuine ell-squared criterion for the remaining finite torsion length

This is a source and method preflight. No ell-thick coefficient, kernel,
rank, implicit jet, or new scalar value is computed. The exact comparison
length11 and Loewy filtration concern the fiber ell=0. Its coordinate rA
is not the DVR parameter ell=5(sA-rA).

Put D=k[[ell]], K=ker(D0->E4), and C=coker(D0->E4). Both are finite over D
and have rank4. The independently accepted length11 fiber calculation
gives mu_D K=7, hence

    K = D^4 direct_sum (sum_{i=1}^3 D/(ell^{e_i})), e_i>=1,

noncanonically as a D-module. The three exponents and their sum are not
determined by that fiber. Their sum is the length of H0(D0).

For the actual pre-specialization comparison R^768 -> R^56 -> R^3, write
J=image(V), L=ker(V), K=L/image(X). For every n>=1, J is ell-regular
because it is a submodule of R^3. Tensor the two short exact sequences for
V by R/(ell^n), then divide by the actual image of X. This gives

    0 -> K/ell^n K -> T_n -> C[ell^n] -> 0,
    T_n=ker(V mod ell^n)/image(X mod ell^n).

Since rank_D C=rank_D K=4, cancellation of the two free ranks gives the
particularly small exact numerical invariant

    d_n = length(T_n)-length(C/ell^n C)
        = sum_{i=1}^3 min(n,e_i).

Thus d_1=11-8=3. One honest ell-squared pair with d_2=3 proves all e_i=1
and exact length H0(D0)=3. If d_2>3, it instead gives the exact number
d_2-3 of exponents at least2; it does not determine their total length.
More generally d_n-d_{n-1} counts exponents at least n. No splitting of
the comparison as rings or R-modules is used.

The right-exact C part has a small finite input. The accepted product
curve type is D^6+(D/ell^2)^7+(D/ell)^7, so

    length R/(IA*IB,ell^2)=33.

The raw cokernel is D^12+(D/ell^2)^7+(D/ell)^7, and its ell-squared
quotient has length45=6+6+33. The actual C/ell^2 C is its quotient by the
twelve physical difference actions, at most a45 by540 QQ coefficient
map (24,300 slots), once the actual ell-squared quotient basis and its
coordinate/source actions have been constructed and independently
checked. This size bound does not construct that missing table. The old
closed raw26 table is only ell=0 and cannot supply it unchanged.

The T part requires an honest ell-squared kernel or a complete compatible
comparison resolution and source lifts. The literal Cbar cyclic resolution
and minimal5-by825 Schur operator are ell=0 constructions and are not an
ell-squared deformation by declaration. Nor may a kernel be computed
after tensoring by IA*IB or its square: those operations add Tor terms.
For example R=k[[ell,x]], V=multiplication by x, X=0 has K=0 and C=D.
Its true T_n is zero. Reducing first by x^2 creates the nonzero kernel
x*(R/(x^2,ell^n)), despite x^2 annihilating the true T_n. Annihilation
of the desired homology is not a kernel base-change theorem.

Available source interfaces are the allorder original9/640 D0 map,
physical free4/Q48 normalizer map, its complete Res4 syzygies, and the
actual UP=QB lift. Existing FULL9 source2 and bivariate source3 cache
contains some ell-thick directions. There is currently no serialized
honest ell-squared minimal comparison or explicit K[ell] source basis.
The next coefficient job must first close one of those exact source gates
and measure its complete live support. No such job is authorized or run
by this packet; no full348/696-base homology shortcut is proposed.

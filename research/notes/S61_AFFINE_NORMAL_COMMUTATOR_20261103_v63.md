# Normal cubic block, reduced origin, and determinant splitter without z

2026-11-03. This continues the frozen finite-five-generator normal-base theorem. It supplies a concrete normal-only cubic block with invertible linear Jacobian and removes the embedded-origin ambiguity from the scalar normal presentation. It does not prove that the unknown residual is empty.

The exact QQ identities z_k*n=P_k(n) have degree(P_k)≤2. Factor each monomial of P_k by its first normal variable to write P_k=R_k(n)n, where R_k is a42×42 matrix affine-linear in n. Since n and z commute, in the actual original quotient

    z_i P_j=R_j z_i n=R_j P_i,
    z_j P_i=R_i z_j n=R_i P_j.

Thus all entries of [R_i(n),R_j(n)]n belong to the original ideal. These are genuine polynomial relations of degree≤3, with a fully specified finite product recipe from the native QQ cross inverse; they are not inferred from modular samples.

At n=0, the three constant commutators have determinants21866,23983,19574 modulo32003. Every coefficient of the exact QQ R_k is p-integral by the preceding native certificate. Consequently all three determinants are nonzero QQ p-units. In particular the first42-cubic block has full normal Jacobian at origin. The verifier independently reconstructs the modular cross inverse, matches all119070 coefficients with the hashed exact QQ output, and saves all three constant matrices and commutators.

Let J be the prior normal ideal with2796 explicit polynomial recipes, and let

    J0=J+(42entries of [R0,R1]n).

Away from origin J already equals the actual residual ideal in every D(n_i) chart. The added cubics vanish on that residual. At origin, the full Jacobian implies (J0)_m=mR_m by Nakayama applied to m/J0. Therefore the global ideal is exactly

    J0 = m intersect I_B,

where B is the actual residual, closed in normal affine42 and disjoint from origin. This is an equality of ideals, preserving the full unknown nonreduced structure. There are2838 generators of degrees≤4 in finite recipe form. The only normal point added to B is a reduced origin.

A scalar normal-only residual equation is also available. Set

    D84(n)=det([R0(n),R1(n)])/det([R0(0),R1(0)]).

The denominator is beta-only and nonzero on the recorded open; degree_n(D84)≤84. Adjugate multiplication gives D84*n=0 in the full actual algebra. At the known five points n=0 it equals1, while on the residual n generates the unit ideal, so it equals0. Hence D84=e, the original global known-five idempotent, and

    I_B = J0+(D84)

is an unsaturated scalar normal-only presentation of the entire unknown residual. Degree exactly84 and an expanded QQ determinant are not claimed. This is a normal-only counterpart to the earlier degree42 determinant in(z,n), not an emptiness conclusion.

The cubic block has invertible linear terms but nonzero higher terms; it does not give42 linear consequences or contradict the complete degree3 Macaulay obstruction. No unchanged degree1 left-inverse search was repeated. Neither properness over beta nor a generic-unit inference from modp emptiness is used.

Reproduction first generates the exact cross output with `code/s61_affine_normal_finite5_20261103_v63.py`, then runs `code/s61_affine_normal_commutator_20261103_v63.py --output <freshdir> --cross-output <thatworkdir>/cross-output.txt`. The107MB exact expansion remains reproducible ignoredwork; its hash is pinned by the finite5 receipt.

# Actual closed raw cokernel26: quotient maps and coordinate actions

This packet constructs the finite input after ell=5(s_CA−r_CA)=0. It uses
the independently accepted actual conductor-product20 ring and its SAME
CA-F9 coordinates (zA_0,...,zA_6,r_CA,s_CA). These are actual functions of R0;
raw physical rho coordinates are not substituted for them.

The previously accepted raw normalizer-image formula gives, by right exactness,

    C_raw/ell = (R0/(IA,ell))*tauA
             ⊕ (R0/(IB,ell))*tauB
             ⊕ (R0/(IA*IB,ell))*tauA*tauB.

Its dimensions are3+3+20. This is a module decomposition; no ring structure on
the raw cokernel is asserted. Both tau symbols retain their own actual trace
conventions. The packet does not evaluate nu, delta, a or b, and does not
construct the actual comparison cokernel or kernel.

## Two canonical quotients of the actual20 ring

Let r mean the actual CA-F9 function r_CA on both quotients. The accepted
coordinate-normalized ideals, to the needed order, are

    IA=(zA_0,...,zA_6, (125/24)r³+O(m⁴)),
    IB=(zA_i+beta2_i+beta3_i+O(m⁴), (−15125/144)r³+O(m⁴)).

All these source arrays and their original physical generator IDs are retained.
The accepted actual20 theorem makes m⁵ zero, so its eight multiplication tables
are exact, rather than merely truncated approximations. The quadratic differences
are beta2_i=c_i r², where

    c=(9835/864,0,0,0,29875/3456,−4565/864,0).

The exact two target algebras are QQ[r]/r³. Their maps are

    piA(zA_i)=0,        piA(r)=r,
    piB(zA_i)=−c_i r², piB(r)=r.

Every beta3 and later term vanishes after these substitutions into QQ[r]/r³.
The source-bound cubic final generators vanish as well. Thus the maps factor
through the actual conductor quotients, whose independently accepted lengths
are3. They are surjective since r is retained, so they are isomorphisms. The
maps from H20=R0/(IA*IB,ell) are the canonical quotient maps, with source
IA*IB contained in both IA and IB.

The fixed20 basis is

    1; zA_0,...,zA_6,r; r*zA_0,...,r*zA_6,r²;
    r²*zA_5,r³; r⁴.

The checker constructs both3×20 matrices from these literal substitutions. The
common chosen lifts of (1,r,r²) are basis positions(0,8,16). Each quotient
composed with this lift is the identity. All eight generator intertwining
identities are checked on each side, as are all400 basis products on each
side. The eight actual IA generators annihilate the A target, and the eight
normalized actual IB generators annihilate the B target. No field solve or
new source coefficient query is used.

## Closed26 action table and callbacks

Order the raw module basis as

    (1,r,r²)tauA; (1,r,r²)tauB; (the accepted20 basis)tauA*tauB.

For zA_i the first3 block is zero, the second3 block is −c_i times
multiplication by r², and the third block is the accepted H20 action. For r,
the two3 blocks are the ordinary length3 shift and the20 block is its accepted
action. The ninth s_CA action equals r only in this closed fiber. All eight
26×26 matrices are saved in column-sparse form and their commutativity is
verified. Every actual IA*IB product annihilates the whole26 module.

The evidence also stores the20 multiplication operators for the20 monomial
basis of H20 and their corresponding26 scalar actions. The code exposes
`multiplication20(vector,evidence)` and `scalar_action26(vector,evidence)`.
Their inputs must already be evaluated actual H20 coefficient vectors. They
do not interpret physical coordinate names or evaluate source functions.
The quotient matrices and chosen lifts are exposed for later cross-sector
actions; an arbitrary chosen lift is not claimed to define a descended map
without its required source annihilator identity.

These finite tables provide exact multiplication callbacks for the next actual
source-function evaluation. They do not provide generic DVR actions, derived
specialization, Tor corrections, a complete comparison presentation, torsion,
socle or D0 depth. The entire calculation is fixed finite rational arithmetic,
under60s/2GB/32MiB/150000 terms, with no new point, prime, jet, GB or rank search.

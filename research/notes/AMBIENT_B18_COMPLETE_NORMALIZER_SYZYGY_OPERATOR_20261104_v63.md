# Complete actual normalizer relation module over the same R0

2026-11-04. This is a finite component-typed unit-implicit representation,
not an expansion of power-series coefficients or a new source solve. It
retains every physical root, row and unit in the accepted source interface.
It closes the representation of complete Syz(Q); it does not decide the
comparison kernel's ell torsion, the D0 socle, or D0 depth.

## 1. Actual coefficient algebra and source row map

The actual algebra is C=R0[tauA,tauB]/(tauA²-nuA,tauB²-nuB), finite free over
the SAME physical nine-coordinate R0, with basis (1,tauA,tauB,tauA*tauB).
Its twelve functions F are the actual source-bound differences in rows
18,19,28,29,60,67,72,63,64,68,73,74. Its accepted Q is the four by 48 matrix
whose column 4i+a is F_i times basis_a. This entire column identity is checked
again, using formal independent coefficient atoms and all sixteen basis
products. No old regular quotient or source image is substituted implicitly.

The preceding physical-interface packet pins the finite actual map kappa9
from C to the proper source: its Jacobian is 5/48. Its exact target row
identity is dH12(kappa)=R12 F, with det R12(0)=687194767360000.
The source adapter T is the factored invertible matrix from the original H12
to (G8,0^4), where G=(u5,Q0,Q1,sr²). All five d*u relations remain. The
current evidence embeds the actual H12 terms, core roots, source cofactor
arrays and every factor definition, not merely the resulting dimensions.
Consequently N=T(kappa)R12 is an actual unit C matrix satisfying

    N F=(dG0,...,dG7,0,0,0,0)^t.

Its inverse is R12^-1 T^-1. Only previously pinned local-unit inverses are
used. In particular d is never inverted. Both own trace roots retain their
original signs and means.

## 2. Literal finite four-component lowering

Every C scalar node x is represented by four R0 components. Constants embed
as (q,0,0,0); sums and negatives act componentwise; the sixteen literal
basis products use only tauA²=nuA and tauB²=nuB. Multiplication by x has
matrix

    [x0  nuA*x1  nuB*x2  nuA*nuB*x3]
    [x1  x0      nuB*x3  nuB*x2    ]
    [x2  nuA*x3  x0      nuA*x1    ]
    [x3  x2      x1      x0        ].

The code verifies all 64 basis associativity identities and all 16 products
of these multiplication matrices. For a known local unit x its inverse is
the unique unit linear system mu(x)y=e0; its determinant's closed-point
value is x0(0)^4. A native C matrix M is restricted by blocks mu(M_ij).
Multiplication, inversion of unit matrices and matrix composition commute
with this restriction. No transposition occurs inside a multiplication block.

A pointed C unit-Hensel system is lowered by replacing each unknown by four
R0 unknowns and applying the same multiplication table to each equation.
Its scalar Jacobian is Res4 of the native Jacobian and is a unit. For the
group gauge, the original thirty-scalar equations arise by first using a
FORMAL auxiliary etaA²=nuA to extract even/odd coefficients of fifteen group
equations. That auxiliary symbol is distinct from the physical tauA while
forming the thirty equations. Base-changing this already defined system
to C and restricting gives 120 R0 unknowns/equations, with Jacobian origin
(4/9)^4=256/6561. The unique solution is the original R0 even/odd solution
embedded in C; substituting etaA=tauA afterwards gives the actual group.
Thus no ambiguous decomposition inside C is introduced. The original 96
pure source roots remain their known R0 roots. The conservative 120-root
representation is redundant but finite and exact.

The thirty scalar origins are not all zero. Their even fifteen values are
zero, whereas the odd fifteen are the frozen first-derivative tauA column
(-1,0,0,2,0,0,-2,0,-2,0,0,-1,0,-4,0). The code pins these and the resulting
120 component origins. This is exact at etaA²=0, using the already verified
source/gauge first derivative. The native group15 is centered at zero modulo
the maximal ideal of C, which does not set every R0 component to zero.

The resulting component values are defined recursively by literal finite
algebraic equations. They are not unspecified coefficient placeholders. The
serialized output gives the typed inputs, component rules and all output
block expressions. It does not claim expanded R0 coefficient polynomials.

## 3. Complete syzygies and the next differential

C is regular of dimension nine: the two monic unit minor is -20/9, and the
actual source coordinate minor 5/48 identifies its chosen source coordinates.
The independently accepted full-core membership theorem makes G8 a regular
sequence in C. Since d is nonzero in this regular domain, the relations of
dG equal the relations of G; this does not divide by d.

Order Koszul pairs (i,j) and triples (i,j,k) lexicographically. K2 has column
G_j e_i-G_i e_j. K3 has column G_i e_(j,k)-G_j e_(i,k)+G_k e_(i,j).
The code checks all 28 entries of G*K2 and all 448 entries of K2*K3 in the
free polynomial ring over eight independent scalars. Let

    S = diag_block(K2,I4): C^32 -> C^12,
    V = [K3;0]: C^56 -> C^32.

The four identity columns correspond to the four zero source rows, not four
discarded equations. Koszul exactness proves that S generates every
relation of (dG,0^4), and V generates every relation of S. Hence the actual
physical relation module is generated COMPLETELY by N^t S. For column F,
the row presentation obeys F^t N^t=(dG,0^4)^t_as_row.

Restriction of scalars is exact because C is R0-free. Therefore the actual
complete operators are

    Q: R0^48 -> R0^4,
    D2=Res4(N^t)Res4(S): R0^128 -> R0^48,
    D3=Res4(V): R0^224 -> R0^128.

Their chain identities are consequences of the ACTUAL source identity and
the checked universal identities: Q D2=Res4((dG,0)S)=0, and
D2 D3=Res4(N^t S V)=0. Completeness follows before restriction, not from
these zero products or shapes alone. Continuing Koszul gives the nonminimal
R0 resolution ranks (4,48,128,224,280,224,112,32,4). Cancelling the four zero
C-generator pairs gives (4,32,112,224,280,224,112,32,4); even this is not
claimed minimal over R0.

## 4. Scope

The exact input is below 32 MiB and retains factored matrix expressions.
No source series solve, new point, jet, coefficient-polynomial expansion or
Gröbner basis is run. This is a complete finite actual coefficient operator
for Syz(Q) and its next differential, which can be evaluated in later finite
coefficient algebras. It is not an evaluated DVR action or a complete
comparison-kernel computation. Those arithmetic gates, including old-relation
Hadamard lifts, remain explicit and are not inferred from H0(E4)=0.

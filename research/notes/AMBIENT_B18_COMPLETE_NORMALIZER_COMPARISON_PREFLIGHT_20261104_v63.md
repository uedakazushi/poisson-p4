# Complete normalizer syzygies and a smaller source-bound DVR route

This is a proof and operator-size preflight. The new actual proper-core
membership theorem and its independent review imply J is a regular
sequence of length8, with core CI degrees(2,2,3). Thus normalized E4 has
depth1 and H0(E4)=0. It follows immediately that H0(D0)=H0(K), where
K=ker(D0->E4). The literal physical adapter and DVR coefficient actions
listed below remain inputs, not computed results of this preflight.

Write R=R0, and let C=B_A tensor_R B_B be the actual regular normalizer
algebra, free of rank4 over R in basis(1,tauA,tauB,tauA*tauB). Denote the
comparison cokernel by M to distinguish it from this algebra. The actual
nine R coordinates and both companion matrices are pinned in the frozen
normalizer-tensor packet; no abstract type substitution is used.

## Complete syzygies of the actual normalizer relation map

Once the two physical source interfaces are accepted and restricted, write the
twelve actual generators as F_old and choose T with

    F_new=T F_old=d*(f1,...,f8,0,0,0,0)^t,

where f1,...,f8 is regular in C. T is invertible and d is nonzero in the
regular domain C. A relation on d*f is a relation on f; d need not, and
does not, become a unit. Thus the complete C-syzygy module is generated
by the28 Koszul pair relations on f and four free redundant-generator
relations. Higher relations are the ordinary Koszul ones.

Orientation matters. For column generator lists as above, the row map
satisfies Q_new=Q_old Res(T^t). A block of Res(T^t) is multiplication by
the corresponding C coefficient T_ij; it is not the transpose of its
four-by-four multiplication matrix. Old coefficient syzygies are
Res(T^t) times the new Koszul/redundancy syzygies.

For c=c0+cA*tauA+cB*tauB+cAB*tauA*tauB, multiplication is
c0*I+cA*muA+cB*muB+cAB*muA*muB. The companions commute and square to
the actual nuA,nuB. Forgetting scalars is exact, so this gives complete
R-syzygies48-by-128 and the next matrix128-by-224. The free resolution
sizes are(4,48,128,224,280,224,112,32,4); these are size bounds, not
minimal Betti claims. Eliminating the four redundant C generators gives
Q_flat:R32->R4 and complete first syzygies32-by-112.

The source-only adapter already has an explicit unit determinant. That
alone is not the physical adapter: the regular-S9 to actual-C map and the
invertible row map from d*H12 to the actual twelve F rows must both be
serialized. MS's subsequent PROPER_PHYSICAL_INTERFACE_FIRSTJET primary
supplies a finite matrix-typed source interface, with Res4 lowering being
handled in its own packet. This preflight does not perform or claim that
coefficient computation. Functions use the same actual Cfree4
basis. A unit regular function maps to an invertible multiplication
matrix; a nonunit must not be inverted through its norm. This is MS's
separate literal-interface task.

## The honest low comparison and its limitation

Let P:R640->R9 present D0 and U:R9->R4 be the actual C0 map. A relation
lift B with UP=Q_flat B is constructed from the68 actual CB zero
polynomials times9 C0 generators. Their degree-at-most-two coordinate
telescoping differences give actual Hadamard multipliers for the twelve
F rows; the28 intrinsic C0 columns have zero lift. Existing old r/g
functions remain R-scalars. None is recomposed as a new0/0 quotient.

The literal first C0 generator maps to1. Eliminating this unit makes the
combined presentation [Q_flat,-U] a3-by-40 matrix V. The low complex is

    R752 -> R40 -> R3,

with112 normalizer syzygies and640 D0 relation columns. Its H1 is exactly
K and H0 is M. It is not a complete comparison resolution. Adding only
the known higher E4 resolution still leaves homology coming from Syz(P).
Its specialized H2 therefore cannot be called K[ell].

## Complete resolution of the small finite-DVR cokernel

There is a route through M that avoids computing Syz(P). Put
O=QQ[[ell]], with the actual regular parameter ell=5(s_CA-r_CA), and
R=O[[x1,...,x8]]. The accepted product-conductor theorem makes M finite
over O. The actual DVR basis/actions, still to be constructed, give a
presentation with at most26 generators and326 relations.

Replace the relation image by an O-free basis A:O^r->O^m, r,m<=26.
This is certified by both inclusions P_O=A B and A=P_O H and a nonzero
r-minor/injectivity certificate, not a claim about a specialized kernel.
It requires no large Smith normal form. Lift each actual x_i action to
T_i on O^m. Relation preservation gives unique S_i and H_ij with

    T_i A=A S_i,             [T_i,T_j]=A H_ij.

Since A is injective, these imply

    [S_i,S_j]=H_ij A,
    cyclic_sum(H_jk T_i-S_i H_jk)=0.

For the latter, multiply by A and use the matrix Jacobi identity, then
cancel A. Thus no unproved higher coherence is needed.

An explicit R-free resolution of M has

    F_n=R^m tensor wedge^n O8
        direct_sum R^r tensor wedge^(n-1) O8.

Let iota_i be exterior contraction. The differential is A from internal
degree1 to0; sum(x_i-T_i)*iota_i in degree0; minus
sum(x_i-S_i)*iota_i in degree1; and minus sum_(i<j) H_ij*iota_i*iota_j
from internal0 to1. The displayed identities make its square zero.
The amplitude is two rows, so two H factors compose to zero.

For exactness, filter by exterior degree. After A, only M remains; the
next differential is the completed Koszul complex for x_i minus the
commuting actual actions on M. These actions are topologically nilpotent
mod ell because M is a finite local O-module. Formal monic division
(equivalently completed translation) proves that Koszul complex resolves
M. This argument permits O-torsion in M.

The worst free sizes are26*binomial(9,n). The first four are at most
26,234,936,2184. They are factored Koszul blocks, not a proposed dense
matrix dump. The actual26 generators can include the three physical
cokernel generators. Their obvious23 redundant target directions can
then be cancelled with unit pairs, leaving G0=R3 and G1 of size at most
211, without changing the upper sizes936 and2184 for G2,G3.

## A complete finite presentation of K without Syz(P)

Write the resulting complete M resolution as
G3 --d3-> G2 --d2-> G1 --d1-> R3 -> M. Its physical presentation is
V:R40->R3. Retain source-bound comparison lifts

    d1*m=V,       V*n=d1,       d2*h=I-m*n.

They come from the actual DVR normal-form/Koszul contraction and source
relation ancestry, not from assigning values to unknown germs. Then

    S=[n*d2, I-n*m] : G2 direct_sum R40 -> ker V

is surjective. Its COMPLETE relation matrix consists of
(d3,0) and(-h,n). Indeed, if S(w,x)=0, set y=m*x-d2*w. Then x=n*y and
d2(w+h*y)=0, proving the assertion with exactness of the G resolution.

For each of the752 known comparison cycles v in ker V, construct w by
the same cycle contraction with d2*w=m*v. Its lift to this presentation
is(w,v). Adding all these lifts gives an actual cokernel presentation
of K. No source Syz(P) calculation is required. The coarse bounds are
at most976 generators and3147 relations over R. These bounds establish
finiteness of the operator recipe; they do not authorize its dense
materialization.

The contraction may be stored through regular Hadamard and unit-source
operators, retaining old r/g as their original regular functions. If a
needed regular-function interface or a cycle-lift identity is not
available, that lift stays an explicit missing gate. This preflight does
not invent those functions or compute their large derivative matrices.

## Finite torsion criterion and resource stop

Only AFTER this complete K cokernel presentation exists, tensor by
R/(I_A I_B). Its right-exact cokernel is still K because the product
annihilates K; no flatness or kernel base-change is assumed. The actual
product-base DVR basis has20 generators and14 torsion relations. The
coarse unminimized bound would be19,520 DVR generators and76,604
relations, too large for an automatic dense job. Unit cancellations and
actual size/dependency counts must come first. Actual DVR actions have
not been serialized by the previous26/326 preflight.

For the resulting exact finite O-presentation one needs only

    dim_QQ(K/ell K) > dim_Frac(O)(K[1/ell])

to detect nonzero torsion. The difference counts torsion cyclic summands,
not their lengths. Since R/(I_A I_B) is finite over O, O-torsion is
exactly H0_m(K)=H0_m(D0). This could decide the remaining depth gate
without a full Smith computation. K/m_R K is a different quotient and
cannot be used in place of K/ell K.

For a genuine comparison complex with only H1=K,H0=M, the DVR derived
special fiber has H2=K[ell] and an exact sequence
0->K/ellK->H1->M[ell]->0. The incomplete752->40->3 complex does not meet
that hypothesis. No K torsion, socle, or full D depth is computed here.

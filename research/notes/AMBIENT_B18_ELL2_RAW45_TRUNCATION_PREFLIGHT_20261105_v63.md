# Exact cubic source truncation for the ell-squared raw cokernel

This is a pure source theorem and finite-input preflight. No quotient rank,
source value, coefficient solve, or ell-squared comparison kernel is evaluated.
Let R=k[[z0,...,z6,r,ell]], with the actual ell=5(sA-rA) and r=rA. After
the accepted unit generator changes write

    IA=(z7,gA(r,ell)), IB=(z7+beta(r,ell),gB(r,ell)),
    ord(gA),ord(gB)>=3, ord(beta)>=2,
    J=IA+IB=(z7,(r,ell)^2).

Both actual conductor curves are D=k[[ell]]-free3. The common r is a
Weierstrass coordinate on each: its cubic leading coefficient at ell=0
is a nonzero rational number. This follows from the saved source cubic
and the actual parameter change, without declaring any nonunit a unit.

## The union and its finite Tor attachment

Put U=R/(IA intersect IB). It is the fiber product of the two conductor
algebras over A3=k[r,ell]/(r^2,r*ell,ell^2). Write each branch using
D-basis1,r,r^2. Matching in A3 means that the constant coefficients
differ by ell^2 and the r coefficients differ by ell. Thus U has the
following explicit D-basis, where the components are the two branches:

    e0=(1,1), e1=(r,r), u=(r^2,0), v=(0,r^2),
    p=(ell^2,0), q=(ell*r,0).

This retains the Tor attachment. U/(ell^2) does not in general inject
into the direct sum of the two branch quotients; that shortcut is false.

On either branch Weierstrass preparation gives

    r^3=A2(ell)*r^2+A1(ell)*r+A0(ell),
    ord_ell(A2)>=1, ord_ell(A1)>=2, ord_ell(A0)>=3.

Let W be the D-span of u,v,p,q. Modulo ell^2 one has r*W contained in
ell*W and r^2*W=0. In particular r*p=ell*q, r*q=ell*u; the A branch
formula for r*u uses A2*u+(A1/ell)*q+(A0/ell^2)*p. The B branch formula
for r*v is analogous, with the q,p terms negative. The extra terms
A1*e1 and A0*e0 vanish modulo ell^2.

Each zi is zero on the A branch and equals minus beta_i of order at
least2 on the B branch. This proves zi*u=zi*p=zi*q=0. For zi*v, reduce
beta_i*r^2 in the B Weierstrass basis: its r^2,r,constant coefficients
have orders at least2,3,4 in ell. In the six-element union basis this
has an overall ell^2, so zi*v=0 as well, including the Tor attachment.
Further zi*e1 belongs to ell*W and zi*e0 belongs to W+ell*k*e1.
Consequently all zi*zj, r^2*zi, ell*r^3, and r^4 vanish on U/(ell^2):

    m^4(U/(ell^2))=0.

## The product and the exact truncation theorem

Let B=(IA intersect IB)/(IA*IB). The accepted source theorem gives JB=0,
and J contains m^2. Since U is ell-regular, tensoring
0->B->R/(IA*IB)->U->0 by D/(ell^2) preserves its left injection. The
previous bound now yields, for S=R/(IA*IB,ell^2),

    m^4 S is contained in B/(ell^2 B),
    J*m^4 S=0, m^6 S=0.

Replace every generator of each conductor by its Taylor polynomial of
degree at most3. The omitted terms lie in m^4. A product changes by
J*m^4 S+m^8 S=0. The same argument applies to the truncated ideals:
their quadratic beta terms span all three plane quadrics, their cubic
generators are ell-transverse, and they have the same A3 union quotient.
Applying the argument in both rings proves the exact ideal equality

    (IA*IB,ell^2)=(IA_le3*IB_le3,ell^2).

The saved beta3_before_ell has z-dependent terms. To use pure-plane
graph generators, factor these degree3 z terms into a unit matrix
I+M2. Removing them changes the pure-plane generators only in degree4,
which the theorem permits. Thus the saved plane projection gives the
correct beta2 and beta3 graph coefficients for this equality; no unknown
fourth source coefficient is set to zero.

## A small exact module presentation for S

Use sA=r+ell/5 and ell^2=0. Eliminate the49 monic products
zi*(zj+beta_j). Over k[r,ell]/ell^2 the remaining module is generated
by1,z0,...,z6, with relations

    beta_j*zi-beta_i*zj (i<j),
    gA*(zi+beta_i), gB*zi (i=0,...,6).

The gA*gB product has degree6 and is zero. The zi actions are
zi*1=zi and zi*zj=-beta_j*zi. The21 pair relations make these commute.
Acting by zi on a pair relation or a gA relation gives zero; acting on
a gB relation gives a plane multiple of another gB relation. Thus the
displayed plane module presentation is closed under all nine R actions,
and no unrecorded overlap equation is discarded.

The m^6 bound gives74 ambient monomials:11 plane monomials of degree at
most5 and63 z-sector monomials of additional plane degree at most4.
The21 pair relations need at most5 plane multipliers; each of the14
remaining rows needs at most3. Hence at most147 rows and74 columns,
with a conservative1470 nonzero input bound. All required saved cubic
source data already exist. No such rank computation has been run here.

The accepted DVR types predict length(S)=33 and raw cokernel dimension
45=6+6+33 after ell^2. Those are independently justified target sizes,
not results of a new matrix calculation. Constructing physical actions
on raw45 would still require a separately gated finite source evaluation.
Constructing an honest T2 or its comparison resolution is a further
independent gate; this theorem alone does not determine H0(D0).

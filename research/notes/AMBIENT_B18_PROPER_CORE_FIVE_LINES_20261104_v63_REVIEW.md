# Independent review: actual five-line core and exact nilpotence three

2026-11-04. Accepted scope: the formal proper-source core and its WHOLE
normalized E4 germ. No D0, full D, ambient component identity or global
extension is inferred. This review has no dependency on the new physical
interface/Res4 operator under concurrent review; both use the already
independently accepted actual core CI as a separate common input.

## Independent coefficient checks

The checker imports no producer code or SymPy/CAS. It reads Q0,Q1,sr²
directly from the accepted preparation, with variables (a,d,x=v2,y=v6).
An independent sparse Fraction Laurent-polynomial interpreter verifies every
r/s restriction, the complete elementary prime case split and all five line
substitutions. It directly replays the N/P two-term source expansions in
QQ[t,t^-1,u], every exact C y-elimination coefficient in QQ[t,t^-1,u,v],
and all entries and signed 3x3 determinants of the R/L normal Jacobians.
No numerical value of t is selected and no implicit equation is freshly solved.

The formal unit P=d²+5 is retained. At r=0, the identities Q0=8h,
Q1=-8aj and dj-5h=aP(a+2d) give exactly N,C,R by separating d=0 from d!=0
in a prime. At s=0, Q0=4d(-a+2d+2y); the two prime cases give exactly N,P,L.
Each converse is a literal substitution, so there is no missing branch.

The N minor is -80t², with d=4u²/t and x=-3u+2u²/t through order two.
The remaining function starts -32u³; its actual local algebra is length
three with HF(1,1,1). The explicit IFT minor proves the higher corrections
cannot alter that leading coefficient. At C the exact y-elimination leaves
v² and a quadratic with u² coefficient -16(t²+5), a regular initial pair.
It has HF(1,2,1), length four. The R normal minor is
-256t²(t²+5), giving r². At P the unit minor is -240t² and s starts
-u²/(3t), giving length two with d=u. At L the normal minor is
16t²(t²+5), giving a reduced point. These computations embed into the
geometric generic fields C((t)); all quoted rational units remain units.
The five generic lengths are (3,4,2,2,1), summing to the CI multiplicity 12.

## Core nilpotence, with the injectivity hypothesis made explicit

The accepted one-dimensional complete intersection is CM and has exactly
its five minimal associated primes. If a finite submodule vanishes at all
these primes, it is zero: any nonzero submodule has an associated prime,
which is an associated prime of the containing module. Apply this to the
cube of the nilradical. The five Artin localizations have maximal-ideal
nilpotence at most three; at N and C its square is nonzero. Thus the exact
global FORMAL core nilpotence index is three. Reduced support and generic
lengths alone would not suffice without this no-embedded-prime fact.

## Embedded localizations of the whole E4

Let S=QQ[[u5,a,d,x,y]], J=(u5,Q0,Q1,sr²). The exact injection
S/J --d--> S/(dJ), with quotient S/(d), puts all five core primes in Ass(E4)
and bounds Ass(E4) by those five and (d). The prime (d) is minimal and hence
associated. Exactly N and P contain d, making those two embedded; C,R,L and
(d) are minimal. This proves equality of the six associated primes.

At either embedded prime q=N or P, the only minimal prime of E4_q is (d),
so its nilradical is (d)/(dJ). Multiplication by d identifies this ideal AS
A MODULE with (S/J)_q. Its multiplication is twisted: the kth nil-ideal
power corresponds to d^(k-1)*(S/J)_q, not the kth ordinary maximal-ideal
power of that core algebra. The exact local models prove d!=0 and d²=0
in both cases: at N d=4u²/t mod u³, at P d=u mod u². Hence nil² is nonzero
and nil³ zero at both embedded localizations. At C,R,L the element d is a
unit, giving the core local models; at (d) the ring is reduced. A finite
submodule of E4 vanishing at all its associated primes is zero. Applying
this to Nil(E4)³ proves its global vanishing, while the embedded localizations
give Nil(E4)²!=0. Localization of a nilradical commutes with localization,
so this last conclusion does not require the global element d itself to
be nilpotent on every branch.

This accepts exactly the claimed nil index three. It does not identify any
five coordinate lines with old physical Poisson strata and does not transport
the nil statement through the unresolved D0-to-E4 comparison kernel.

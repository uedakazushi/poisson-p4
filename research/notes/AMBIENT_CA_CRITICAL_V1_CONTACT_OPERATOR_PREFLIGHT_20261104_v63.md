# Existing V1: a single relative contact operator

2026-11-04. This packet computes only constant source directions from frozen
matrices. It does not evaluate the new source coefficients or decide kappa.
The proposed coefficient job uses the SAME old V1 point. It must include all
150 physical outputs and the nonlinear CB q*(-f) term. No new source IFT
coefficients, generic point or field computation are involved.

## The exact N fold and actual source coordinates

In the old CA G20-plus-fiber9 source the third binary coefficient is literally
source coordinate22. The fixed fiber is linear, its center coefficient is zero,
and group action transports the marked flag. Thus N is h_N=0 exactly, and
dh_N(kappa)=1. The actual u8 are the selected eight even odd-coefficient
functions b_i(R) in A; they are not arbitrary coordinates of its normalization.
Together with the actual g20 outputs and tau they have a unit source Jacobian.

The previously proved N-adapted fold tau-s_N(g,u) differs from h_N by a unit.
Changing p to pU(g,u,p), for any unit U, preserves the Ac8 subring

    A=k[[g20,u8,p^2,u1*p,...,u8*p,p^3]] inside B=k[[g20,u8,p]].

Indeed its conductor (u8,p^2)B is an ideal of B lying in A, so all new
generators lie in A. The inverse fold change has the same form and proves
the reverse inclusion. This keeps g20 and the ACTUAL u8 fixed. Consequently
one may use p=h_N, with no higher trace/root readback.

An odd term in U can change the cotangent x by a multiple of Z. It does
not change the cotangents u,Y,Z because U(0)=1. The checker verifies that
Zlinear annihilates the ENTIRE actual CB29 tangent. Therefore the saved
adapted CB even29 projection is also the correct differential for x=h_N^2;
no unknown higher fold coefficient is silently set to zero.

## A conditional finite contact lemma

CB is regular29 here and projects with unit Jacobian to (g20,u8,x=p^2).
The actual N is its smooth divisor x=0. Thus its remaining canonical CA
coordinates have unique expressions

    Y_i=x*b_i(u,x),       Z=x*c(u,x),

with b_1(0)=20/9 a unit and c(0)=0 (index1 means selected physical row60).
Every additional CA ambient graph pulls back to x*G_l; G_l(0)=0 because
the entire CB tangent lies in the complete CA tangent38.

Write v=u1/b1, l_i=u_i-v*b_i for i!=1. The pulled-back complete81
intrinsic cusp relations, after factoring x, generate EXACTLY

    (l_7, X=x-v^2, v*(v-c)).

For proof, the Koszul relation with index1 gives -b1*l_i; the square
Y1^2-x*u1^2 gives b1^2*X; and x*Y1-u1*Z gives b1*(x-v*c).
All other Koszul/product relations reduce to zero modulo these three
groups. This uses only b1 being a unit. The eight equations (l7,X) are
unit IFT equations. After them, put ctilde(v)=c(v*b,v^2), and similarly
Gtilde_l(v). The full divided ideal is

    (l7,X,v*(v-ctilde(v)), all Gtilde_l(v)).

If ctilde'(0)!=1 the third intrinsic equation is a unit times v^2.
Therefore a nonzero linear coefficient of ANY Gtilde_l cuts the contact
algebra down to k; if all these coefficients vanish, every remaining
Gtilde_l lies in (v^2), so the contact algebra is EXACTLY k[v]/v^2.
No higher graph coefficient can change this conclusion. If ctilde'(0)=1
the test stops inconclusive. G20 invariance permits the simultaneous formal
group product and reduction to this transverse computation.

In the first case the actual transverse intersection is S/(x*m_S),
with nilradical length1; in the second it is
S/(x*(u7,x,v^2)), with nilradical length2. The map to the actual normalizer
tensor sends x to t^2, where t^2=x and v*t=x. The intrinsic final equation
and the unit 1-ctilde'(0) imply v*t^2=t^3=0. The normalizer tensor has its
two independent nil classes t,t^2, with v*t=t^2; extra graphs x*G_l vanish
there. The actual image of the nil class x is nonzero t^2. Hence the
normalizer comparison kernel has length0 in the first case and length1
in the second. These are conditional kappa=0/1 statements, not current
coefficient results. They concern the old dense CA-critical V orbit only.

## The single source vector and why it suffices

Let beta=b(0)=nY/nx=(20/9)*nY. The checker chooses a_CA by the fixed29
Jacobian with dg=0, du=beta, dh_N=0. Choose a_CB by the same target data
and k_CB by dg=du=0, dx=1. It checks CB_J*a_CB=CA_J*a_CA on ALL150 rows.
Using the saved kappa mixed Hessian and W=CA_coeff(p^2), it verifies exactly

    V=CB_J*k_CB=CA_H(kappa,a_CA)+W+CA_J*r.

All vectors and constant inverses are saved. Proposed new coefficient work:

* CA on p*kappa+v*a_CA through total degree3: ten monomials.
* CA on v*a_CA+t*r through degree2: six monomials.
* CB on v*a_CB+x*k_CB through degree2: six monomials, retaining q*delta_f.

Define the ONE vector

    Delta=CB_coeff(v*x)-CA_coeff(p*v^2)-CA_coeff(p^2*v)-CA_coeff(v*t).

The constant complement projection Pi to physical150/CA_tangent38 gives
ALL the ambient graph linear contact coefficients Pi*Delta. Taylor's
formula proves this: CA source q1=v*J*a and q2 has p^2*W+p*v*H;
the p*v^2 and p^2*v cubic coefficients give respectively the quadratic
graph pairings with H and W, and the final mixed2 coefficient supplies
the pairing with J*r. Every intrinsic cusp quadratic vanishes on this
pair, so this graph class is independent of the chosen graph representatives.

The saved canonical linear Z functional annihilates CA_J, W and H, and is
one on the p^3 column. Since the exact source Z=h_N^3=p^3, the same Taylor
argument gives ctilde'(0)=Zlinear*Delta. No second CB reparametrization is
needed: Z and every CA graph have zero differential on the full CB tangent.
All source accelerations beyond these coefficients are absent from the
identity, rather than set to zero.

## Bounded execution target

The entire new output budget is 22*150=3300 rational slots. The CA ordered
20 exponentials and CB nine flag exponentials are evaluated in truncated
two-variable rings; each physical coefficient has at most ten slots. Store
all three complete150 outputs, verify cached source1 and CA fold/mixed2,
and record the CB bilinear ratio contributions explicitly. Enforce one
60-second/2GB/32MiB/150000-expanded-term cap, checkpoint before evaluation,
and stop if c_v=1. No automatic extra source order follows an inconclusive
result. Neither D0 origin depth nor a finite-origin socle is claimed here.

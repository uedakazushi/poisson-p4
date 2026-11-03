# Raw QQ surjectivity on the entire finite gcd-exception cover

Let p=32003 and V=Z_(p). Use the literal raw24-by56 matrix of the four actual
cleared quadrics times six light linears, with core coordinates(a,e,P). This is
not a later Schur residual. Its source SHA is
f9590655c56c4c7790d43bc780950ff9e7bf4a3339a4a0879e2901369a3b78aa.
All4824 rational coefficient denominators are p-units.

The rank-two gcd-exception proof defines DB=20P²+3Pe+20P+180a and a quartic
N(P,e), with e^4 coefficient72442188. Set a=−P*(20P+3e+20)/180. Its necessary
Q3-matrix minor, after dividing ONLY known units P, L=20P+3e+20 and J=5P−6e−4,
and a nonzero rational constant, becomes a polynomial M(P,e) of total degree9.
On the original core, L=−180a/P is a known unit. Thus every rank-two common-factor
exception belongs to the closed finite QQ cover

Sigma: DB=N=M=0.

The result here is that the raw matrix is surjective over the coordinate algebra
of ALL this cover, and hence at every complex point of Sigma. It does not prove
raw surjectivity on the whole three-parameter core. The pair and five-quadric
rank gates remain separate.

## A finite proper V-model, including every localizer

Take S=V[P,e]/(N,M), with a replaced by the above polynomial. An explicitly saved
QQ identity is

R(P)=A(P,e)*N(P,e)+B(P,e)*M(P,e).

R has degree28 in P. Every coefficient of N,M,R,A,B is p-integral; both the
leading e coefficient of N and the actual leading P coefficient of R are
p-units. Normalizing gives monic equations of degrees4 in e and28 in P inside
the defining ideal. Consequently S is a FINITE V-module, spanned for example by
e^i P^j with i<4,j<28. The original computation used a resultant, but the proof
uses the explicit identity and p-unit coefficients, not a black-box resultant
membership assertion. No flatness or freeness of S is assumed.

Define the FULL repository localizer in S by

F=500a+90e²+3e−24,
B0=(50a+18e²−9e−2)P+(60e−40)a,
C=100a+18e²−9e−2,
U=a*P*F*B0*C*J.

The finite special algebra S/pS has a verified monic5-generator GB and a28-monomial
basis. Exact independent polynomial arithmetic checks both ideal inclusions,
all10 S-pairs, reduction, every literal model coefficient and an explicit identity
Abar*Nbar+Bbar*Mbar+Cbar*Ubar=1. Thus U is a unit on the ENTIRE special algebra,
not just on a selected open part. The cokernel of multiplication by U on S is a
finite V-module with zero reduction modulo p; Nakayama makes it zero. Therefore
U is already a unit in S. This closes the escape issue for all localizers and
avoids any unproved clopen peel or localization-specialization inference.

## Direct finite-algebra raw-matrix check

The saved28-by28 multiplication matrices for P,e are built by complete normal
form reduction in the verified special GB. They commute. Native FLINT arithmetic
checks N=M=0 as operators and explicitly inverts the FULLU operator. It then
substitutes a=−P*(20P+3e+20)/180 into EVERY literal raw-matrix coefficient.
The resulting Fp-linear map has shape672-by1568 and rank672. A selected672-by672
inverse is saved, and its exact product with the selected source block is the
identity. This checks all24 source components over the full28-dimensional
algebra, including nilpotents. It does not assume the earlier global field
annihilator/saturation calculation.

Let K be the raw cokernel over S. Since S is finite over V, K is finite over V.
The block-rank certificate proves K/pK=0. Nakayama therefore gives K=0; tensoring
with QQ proves the asserted raw QQ surjectivity on the ENTIRE closed Sigma.
This is a proper finite-spread argument. A localized finite-field unit assertion
on a nonproper parameter open would not imply the same QQ statement.

## Saved evidence and scope

The verification-only public driver checks the QQ elimination identity and the
known-unit removal defining M, then independently verifies the field GB, source
lifts, S-pairs, U inverse, multiplication operators and complete raw block rank.
Its author replay passes in3.13seconds. The native block stage takes0.516seconds,
RSS37508KiB. All source data, explicit certificates and selected inverse are
retained, without ELF files. The archived generator also records the corrected
Cloud saturation API, but no saturation result is a premise of this proof.

Restore and replay from the repository root:

```sh
tar -xzf reports/S61_NONSQUARE_FINITE_GCD_EXCEPTION_NAKAYAMA_20261103_v63/finite_nakayama.tar.gz
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_finite_gcd_exception_nakayama_20261103_v63.py
```

This proves the finite raw-matrix statement. Translation to the original product
rank41 uses the repository's separately certified seventeen kernel directions
and actual-coordinate bridge. No whole-core rank24/product41, quartic coverage
or ambient removal is claimed here. Independent review is pending.

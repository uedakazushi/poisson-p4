# Actual proper-source core: one bounded exact membership computation

This is a new packet. Earlier frozen unresolved-core records remain unchanged.
The hypotheses are the independently accepted proper polynomial model, first
five rational pivots, and the subsequent five unit-linear core pivots. In
particular their source proofs retain all fourteen source variables, all
seventeen equations, and the whole normalized model's five multiplied pivot
relations. This packet concerns normalized E4, not the unnormalized D0.

## Result

The single authorized computation succeeded. The default M2 basis has 16
elements; all four graph remainders are exactly zero. All twelve quotient
coefficients and all four unit-cleared source identities were replayed over
QQ. The combined evidence stores 18,136 polynomial terms. M2 took 2.733
seconds; the complete preparation, computation and replay took 21.465 seconds,
with peak child RSS 227,824 KiB. No retry was made. The fixed initial-form
regular-sequence gates also passed, so the full actual seven-function core is
a one-dimensional complete intersection of multiplicity 12, CM and Gorenstein.
The whole normalized E4 has depth one and no finite-length submodule.

## Exact finite target

Use the seven literal polynomials K0,...,K6 in the four-variable core input.
Only D=1-d² and P=d²+5 (origins 1 and 5), and nonzero rational constants,
may be stripped. Write K0=U0 Q0, K1=U1 Q1 and K2=U2 r Q2, where each Ui is
a signed product of D and P, r=v2-a-d, and s=a+v2. The verified identity

    Q0 - 2 Q2 = -P s r

identifies the localized first-three ideal with (Q0,Q1,s r²). Macaulay2
1.24.11 is invoked once over QQ[a,d,v2,v6], GRevLex, using default gb and a
change matrix. Its four quotient/remainder identities are replayed in exact
QQ polynomial arithmetic. Global nonzero remainders or timeout have no
negative local implication. No saturation retry is part of this packet.

When all four remainders vanish, write Kj/Uj=A Q0+B Q1+C s r². Multiplication
by the common unit M=P U0 U1 U2 gives the literal source identity

    M Kj = Uj [U1 U2(P A-r C)K0 + P U0 U2 B K1 + 2 U0 U1 C K2].

All four such identities, including every coefficient, are independently
replayed. These are positive global polynomial certificates with a multiplier
whose value at the origin is nonzero; thus they prove the actual local
seven-generator ideal equals its first-three ideal.

## Conditional geometric consequence and the fixed arithmetic gates

The report records the successful positive membership target.
For the CI consequence, the script separately verifies that the two initial
quadrics f2,g2 have gcd 1, and that their restrictions to r=0 and s=0 each
have gcd 1. Thus QQ[a,d,v2,v6]/(f2,g2) is a two-dimensional homogeneous
complete intersection. Adding either r or s has height three, so both are
nonzerodivisors; s r² is a nonzerodivisor as well. The initial forms are a
regular sequence of degrees (2,2,3). The complete filtered ring argument
(successively compare lowest nonzero homogeneous parts) then gives a
one-dimensional complete intersection, Cohen--Macaulay and Gorenstein, of
multiplicity 12 for the actual first-three local ring. Positive membership
therefore proves the same statements for the full seven-function core.

The whole normalized model still retains

    E4 = S / d*(u1,...,u5,K0,...,K6),
    S = QQ[[u1,...,u5,a,d,v2,v6]].

Writing J=(u1,...,u5,K0,...,K6), multiplication by d gives the exact sequence
0 -> S/J -> S/dJ -> S/(d) -> 0. If the preceding positive gates pass, the
first term is the one-dimensional CM core; S/(d) is regular of dimension 8.
Consequently H^0(E4)=0 and depth(E4)=1. This does not make the whole E4
one-dimensional or CM, and it does not determine H^0(D0), T0, or full D depth.

## Reproduction and bounds

Run the matching Python source with a new --output directory. The report
records the phase file, exact M2 script/log and polynomial output, preparation,
source/cofactor evidence, hashes and actual resources. The single M2 process
has a 45-second cap within the 60-second whole budget, a 2 GB address limit,
32 MiB per-file limit and 150,000 expanded evidence-term guard. No source jet,
new point, parameter substitution, large Poisson GB, or QQ LinearAlgebra
backend is used.

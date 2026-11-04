# Rational coordinates for the FULLU Theta branch

Date: 2026-11-04. Status: exact QQ primary certificate; independent review
pending. This packet verifies a scheme coordinate change only. It computes
no heavy-block rank, minor, B-coordinate defect or additional source pivot.

Use `R=QQ[a,e,d,1/(adFBCJ)]` and the source polynomials of the frozen
p-branch profiles. Let

```
g1=(3e-2)d+20,  g2=(6e+1)d+40=2g1+5d,
Theta=B^2-C*g1*g2,  N=dB-g1*g2.
```

There is a literal square-completion identity in QQ[a,e,d]:

```
N^2-20*g1*g2^2 = d^2*Theta.
```

On R/(Theta), `C*g1*g2=B^2` is a unit, so both g1 and g2 are units. Define
`z=N/g2`. This is a regular coordinate on this scheme, with a representative
whose denominator is already a FULLU unit:

```
z = N*C*g1/B^2 modulo Theta,   z^2=20*g1.
```

The difference between these two representatives, after clearing g2 B^2,
is `-N*Theta`. In particular no arbitrary factor is being inverted here.

Put

```
L=z^2+50d,
K=(z+20)(z-10)+50d.
```

The inverse coordinate formulas are

```
e = 2/3 + (z^2-400)/(60d),
a = (z+20)*K/(500d^2),
g1=z^2/20,  g2=L/10,
B = L*z*(z+20)/(200d),
C = L*(z+20)^2/(200d^2).
```

All forward substitutions into the original polynomials are checked by full
QQ coefficient comparison. Conversely each of the six coordinates
`e,a,B,C,g1,g2`, with `z=N/g2` substituted, is checked after clearing its
displayed denominator. The resulting difference is saved with its explicit
quotient by the original Theta. The other composition is the exact identity
`image(N)=z*image(g2)`. These are scheme identities and make no pointwise,
reducedness or irreducibility assumption.

The remaining FULLU localizers have small translated numerators:

```
Fnum=40d^2*F
    =z^2(z+20)^2+d*(82z^2+2000z+7200)+720d^2,
Jnum=100d*J
    =z^2(z+20)+50d*(z+4).
```

They have seven and four terms, respectively. Thus the exact target ring is

```
S=QQ[z,d,d^-1,1/(z*(z+20)*L*K*Fnum*Jnum)].
```

Every source FULLU localizer maps to a unit in S by the displayed formulas.
For the reverse implication all target localizers are units in R/(Theta).
Explicit witnesses, in the notation above, are

```
g1^-1=C*g2/B^2,             g2^-1=C*g1/B^2,
z^-1=z*C*g2/(20B^2),       (z+20)^-1=g2*(z+20)/(20d^2*C),
L^-1=C*g1/(10B^2),         K^-1=(z+20)/(500d^2*a),
Fnum^-1=1/(40d^2*F),       Jnum^-1=1/(100d*J).
```

The coordinate maps therefore give the actual ring isomorphism

```
R/(Theta) ≅ S.
```

This identifies this auxiliary coefficient branch with a rational open
surface. It does not identify a new Poisson component or extend the result
to the physical a=0 locus. The old monic Theta profiles and one zero
extra-column defect remain intact; their remaining rank questions are not
solved by rationality.

Producer: `code/s61_nonsquare_Theta_rational_coordinates_20261104_v63.py`.
Evidence: `reports/S61_NONSQUARE_THETA_RATIONAL_COORDINATES_20261104_v63`.
The fixed identity run took 0.013 seconds under 57-second, 2-GB and 32-MiB
guards, with no solver, new point, prime or GB. The report pins the producer
and p-branch source hashes. Certificate SHA-256:
`a173a6f33bf18a8ba3292600e534c38cc6dc6f79a1885bae29aa758ba058ad4e`.

# P100 dual: the missing raw-H direction does not descend through NF117

Date: 2026-11-04. Exact QQ(z) primary certificate; independent review
pending. This uses the entire rational P100 chart, without evaluating a
parameter point.

The accepted source chart is

```
a=-10(z+20)/(9z^2), e=(z+100)/(9z), d=-3z(z+20)/100,
```

over the open ring
`C=QQ[z,1/(z(z+20)(z+60)(203z+2700)(z+12))]`.

## Actual raw-H dual

The seven-component dual is constructed from the fixed unit operators, not
by an unknown minor search. After the Theta row8 pivot, put
`n=37500C_core^2*g1` and form the B^4-cleared last column. Use the selected
three-by-three B adjugate divided by the already proved unit Delta to clear
the B coordinates. Lift through the saved second raw-H pivot, then through
the original Cdelta pair change and constant GL7 frame. The original
physical v3 component is the known unit `-a*g1*g2*B^4`; dividing by it
normalizes that component to one.

The resulting values, in original physical light order
`[v0,v1,v2,v3,v4,v6,v7]`, are

```
eta0 = 95/(2(z+20)),
eta1 = -(179z+11460)/(4(z+20)(z+60)),
eta2 = -425/(6z(z+20)),
eta3 = 1,
eta4 = 40/(z+20),
eta6 = -(429z+26460)/(4(z+20)(z+60)),
eta7 = 2125/(12z(z+20)).
```

Native exact univariate rational arithmetic checks that every one of the
sixteen actual raw H-linear source rows annihilates this vector. All
coefficients are read from the previously source-bound raw16-by-7 cache.
No sample, interpolation or field-rank computation is used.

## Six source quadrics and the descent gate

The original six quadrics are the exact source matrix also used by the
accepted joint L,ell frame. Their 174 coefficient entries have 2280 core
terms. All six are evaluated at eta using the full original v coordinates;
the H-linear coefficient of each quadratic is checked literally zero.
Their values are

```
q0(eta)=308125/(12z^2(z+20)^2),
q1(eta)=137500/(243z^2),
q2(eta)=(68750/243-25025z/1458)/z^2,
q3(eta)=(5500000/2187+5500z/243)/z^3,
q4(eta)=984500/(729z(z+20)),
q5(eta)=(308125000/729+38288125z/729+611125z^2/162)
        /(z^5(z+20)).
```

In particular q0, q1 and q4 have known-unit values on the whole curve.
The five eliminated H-quadratic directions span the original q0 through q4:
the first four prime generators span q1 through q4 by the accepted unit
source change, and the extra H pivot is q0H. Thus all five fail the
annihilation test. The functional obtained by squaring eta **does not
descend** to the NF117 H23 quotient.

As authorized, the job stops this quartic-dual route at that gate. No H23
dual is constructed, and neither cached q5H row112 nor quartic row118 is
evaluated. Both were merely input-bound for the conditional branch. The
full set of 21 quartics and the pure98 part remain unchanged.

There is a useful subsequent source-bound route: q1H is killed by NF117,
while q1(eta) is a unit with inverse `243z^2/137500`. The actual raw-H row
module on this curve is `ker(eta)` by the accepted exact unit-cokernel
reduction. Its products with all seven lights span the kernel of quadratic
evaluation at eta. Therefore the quadratic relation can potentially restore
the apparently missing H23 direction. A finite actual cofactor certificate
for that positive application is a separate next packet; it is not inferred
here solely from a numerical rank or a nonzero evaluation.

## Preserved implementation repair

The first attempt stopped at the raw-H annihilation assertion before any
quadratic evaluation. The saved fraction-free checkpoint retains pivot row1
unnormalized. The first lift omitted division by its known unit entry
`-BC_core/14062500`. The repair explicitly verifies that entry and divides
the deleted coordinate by it. Source, native input, executable and failure
receipt remain in
`S61_NONSQUARE_P100_H1_DUAL_FIRST_ATTEMPT_20261104_v63`; no result is claimed
from that attempt. The repaired attempt uses exactly the same fixed source
data and scope, with no new search or enlarged system.

Producers: `code/s61_nonsquare_P100_H1_dual_20261104_v63.py` and `.cpp`.
Evidence: `reports/S61_NONSQUARE_P100_H1_DUAL_20261104_v63`.
The repaired job took 3.532 seconds, peak univariate polynomial length 31,
and 86 output terms under the 57-second/2-GB/32-MiB/150,000-term guards.
There was no full134 Gaussian, POWER expansion, new point, prime or GB,
and no full cubic57/quartic239 assertion.
Certificate SHA-256:
`16ceea3aa59e441c98a2e855437b150ec8d3c26fc6ab435558e4b05af51baf6f`.

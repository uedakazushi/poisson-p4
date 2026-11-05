# Root source review of the T1-only C Hasse job

The7595 integer coefficient lines and common denominator are byte-identical
to the accepted T1 output. Put U=T1/a³, E=−18e²+9e+2 and C=100(a−E/100).
The planned layer is
`D*100^20*H_k = sum N_ied*binom(i,k)*100^(20-i+k)*E^(i-k)*e^e*d^d`.
The root reviewer independently derived this formula, all input bounds and
support/ownership. Since `||E^h||1<=29^h<=100^h` and `binom(i,k)<=2^20`,
every term and every partial sum fits the518-bit positive envelope. The
input has actual l1 bits328 and maximal coefficient bits322. The fixed
ownership35,188 counts all7595 records,441 power coefficients,1558 possible
layer terms and16,000 metadata/scratch slots.

The original root report's product counter denotes native convolution-inner
products. Its separate unchanged-report addendum includes scalar-prefactor
and root-power overhead:884,054 total integer products remain below the
producer's conservative3,972,185 bound. No layer was evaluated in this review.

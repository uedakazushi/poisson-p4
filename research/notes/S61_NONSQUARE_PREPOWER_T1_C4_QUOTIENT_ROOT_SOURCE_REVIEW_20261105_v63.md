# Root review of the finite C4 quotient operator

For U1=T1/a³ and R=−18e²+9e+2, expansion of(100a−R)^−4 at infinity gives
`D*100^20*Q_i=sum_(j>=i+4) I_j*binom(j-i-1,3)*R^(j-i-4)*100^(20-j+i)`.
The independently established C valuation4 ensures a polynomial quotient.
The root review recomputed the positive l1 bounds:472 bits for the full
quotient and500 for every multiply-back remainder, within the frozen518.
It checked all7595 copied integer lines, support10013, partial slab1406,
remainder14991 and full live ownership81310<81776. High-e terms cannot be
removed until an entire slab finishes. The code releases the root cache
before building the fixed divisor cache, then checks every coefficient of
`C^4*Qint−100^20*U1int`, which also detects any missing Laurent remainder.
No quotient value was evaluated in this pre-execution review.

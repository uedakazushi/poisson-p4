# Independent Horner proof that the exact C valuation of T1 is four

The root reviewer recomputed the first five Taylor layers using descending-a
Horner on five truncated jets. The recurrence multiplies only by
`E(e)+100*y`; it neither imports nor reproduces the native binomial/power
expansion engine. All four lower layers are literally zero and every one of
the665 fourth-layer coefficients agrees. Horner uses denominatorD*100^21,
exactly100 times the native scale. The518-bit native limits are unchanged;
this different peer representation has a proved526-bit stored/1580-bit
temporary envelope,44,112 owned slots and less than1.09 MB RAW.

As C=100(a−E/100) and a is coprime to C, this proves v_C(T1)=4 in
QQ[a,e,d]. C is already a named unit on FULLU. No explicit quotient,
other factor or tail, general gcd, bad ideal or strict chart enlargement was
computed by this certificate.

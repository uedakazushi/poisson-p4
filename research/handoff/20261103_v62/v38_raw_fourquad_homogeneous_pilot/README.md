# Raw four-quadric homogeneous field pilot

Extract after v36. The literal QQ 24-by-56 raw source has 4824 terms and
20/17-bit maximum numerator/denominator heights. Each entire source column
is homogenized to its fixed QQ maximum degree BEFORE reduction to GF(32003).
Every one of the fifty-six column degrees is seven. Target row shifts are zero,
variables are (a,e,P,w), and the order is `(dp,C)`.

The single 60-second / 2200000-KiB pilot completes in 10.75 seconds, maximum RSS
42004 KiB. The full monic homogeneous reduced basis has shape 24-by-925,
291802 terms and maximum degree twenty-four. Every original source column
reduces to zero. A separate native checker audits all leading coefficients,
column homogeneity and tail divisibility. The complete field basis is retained.
The generator uses a literal complete entry-by-entry native export, not a
matrix write that could silently preserve only the first row.

This field basis does not prove a QQ standard basis, homogeneous module equality,
localized surjectivity or original rank41. A separate direct QQ attempt reached
its sixty-second cap before producing a basis and is preserved in v37.
Any reconstructed candidate still requires exact QQ GB/S-pair checks, exact
original homogeneous-module containment, p-integral monic coefficients and exact
same-good-prime homogeneous equality before the graded Hilbert-function bridge
can prove the QQ statement. Exact QQ localized membership remains a further gate.

This immutable packet contains one prime only. A separately authorized two-prime
sequential continuation will preserve its own support/leading and reconstruction
comparison, with no automatic extension or QQ promotion.

# Auxiliary five-quadric rank gate: third safe pivot and bounded frontier

This continues the5-by21 coefficient matrix [Q0,Q1,Q2',Q3,D] of the separate
four-quadric structural packet. These pivot counts refer ONLY to this auxiliary
matrix, not the main cubic/product matrix. Its first two safe pivots use only
a,P,J and leave the pinned3-by19 native residual.

A further exact scalar-row identity is now independently verified overQQ. For
residual row1 (original five-quadric row3, namely Q3), with zero-based columns,

(55000/3)*M[1,5]−(11000/3)*M[1,7]=a³*P*J².

The two cofactors are nonzero rational constants, the right side is a known
unit, and the replacement coefficient55000/3 is a rational unit. An invertible
column replacement and unit pivot therefore give a third safe pivot and a2-by18
residual. The exact original-row identity, native output and replay logs are
saved. No residual unknown factor is divided. The log's COMMON_DEGREE4 is the
denominator-clearing factor degree; the actual target above has total degree6.

The3-by19 source is32KB/1844terms, core degree8, coefficient heights21/30bits.
The fraction-free third pivot produces89KB/4278terms, degree12, heights27/41bits.
To avoid this coefficient growth, the one authorized localized QQ module-GB
attempt used the earlier3-by19 source and only z*a*P*J−1. It timed out at59.81s
before a GB-done marker (RSS178700KiB). No QQ unit or nonunit certificate resulted,
and no longer retry is included. An initial0.033s initializer error is separately
retained; cumulative time remains below60s. No scalar or module GB extension is
performed.

Primitive-entry constant-span tests on the two remaining rows refused. A distinct
bounded degree1 cofactor search then used all four monomials(1,a,e,P), giving72
unknown coefficient directions per row, with known-unit monomials a^i P^j J^k
of factor-count at most8 as targets. Both searches finished in0.033s, each with
rank72,165 targets tested and zero modular target passes, hence no exact QQ solve.
These are negative searches in their stated primitive/known-denominator class,
not a proof that a QQ unit row or rank5 fails. A fresh fixture(e,P*e−1) verifies
the newly included P multiplier by the exact identity P*e−(P*e−1)=1.

The full five-quadric rank gate remains open, as do the uniform gcd and pair
gates. No all-core rank24, original product41, quartic or ambient result follows.
An exact selected residual3-by3 minor with a large nonunit factor is retained
as exploration only and is not inverted. Independent review of this frontier
packet is pending.

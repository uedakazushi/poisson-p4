# The three actual remaining quadratic L-minors

This packet evaluates exactly the three two-by-two determinants of
the actual remaining quadratic rows labelled 5 and 9, in original
target columns 21, 22 and 25. The preceding GL6 completion touches
only no-L pure targets and uses a no-L unit row. Therefore all six
entries of this L block are literally unchanged checkpoint references
in the new `16x23` presentation. The code checks that binding in its
saved scalar DAG, rather than assuming a projected surrogate matrix.

The input has 15 and 43 numerator terms in its two rows, rational
constant denominators, and degrees at most two and five. The three
raw multiplication support bounds sum to 414. Exact calculation
took 0.134 seconds. The three determinants are nonzero, with 18, 18
and 22 numerator terms and degree six. Only `a` was removed from
the first two numerators; the last numerator has no named unit
factor. None of the three remaining polynomials is a constant.

The certificate saves every input coefficient, all three actual
determinants, all exponents of the only permitted factors
`a,d,F,B,C,J`, and exact reconstruction identities. No remaining
factor was factored, inverted or declared a unit. The absence of an
individual named-unit determinant is not a rank-drop result for
this block or for the full presentation. No point, GB, Gaussian,
additional minor, or determinant expansion of Delta34 was used.

The subsequent fixed three-term Bezout row identity is a separate
packet. It adds a unit pivot even though none of these three
determinants is an individual named unit. This determinant frontier
remains immutable; it does not claim that its ideal is proper.
Whole-FULLU cubic57, full quartic239, and kernel38/gap42 are unchanged.

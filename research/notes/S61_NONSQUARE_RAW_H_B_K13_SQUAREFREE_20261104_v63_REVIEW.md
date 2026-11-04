# Independent review: K13 squarefreeness and geometric reducedness

2026-11-04. This reviewer performs no gcd, factorization, source solve,
new point/prime, CAS or producer run. The fixed source g13 is read directly
from the independently accepted D_EUCLID checkpoint. Its 14 coefficients,
degree 13, derivative and derivative content are recomputed with integers.

An independent dense integer-convolution engine verifies every coefficient
of the supplied full identity S*g13'+T*g13=R, with R a nonzero integer.
The earlier primitive-derivative identity is separately replayed, including
the exact content multiplication returning to the full derivative. Thus
g13 and g13' generate the unit ideal over QQ, without assuming that the
degree-13 quotient is a field or searching for factors.

Since g13 has degree 13 and nonzero rational leading coefficient, QQ[a]/g13
has dimension 13 over QQ. The verified derivative unit proves this algebra
is finite etale over QQ. Over an algebraic closure it is the product of 13
distinct residue fields equal to that closure. Combining ONLY with the
separately independently accepted ACTUAL_SOURCE_SCHEME isomorphism proves
the specified actual B-rankbad scheme consists of 13 geometrically reduced
points. It does not say these are 13 QQ-rational points, does not establish
QQ irreducibility, and does not transfer the result to the full cubic matrix
rankbad locus, full Poisson scheme or its nilradical.

The report binds the immutable source/checkpoint and actual-scheme review,
and saves the exact nonzero integer remainder with term/bit counts.

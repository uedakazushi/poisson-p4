# Forty-four exact homogeneous degree-seven seed identities

Let M^h be the literal 24-by-56 raw four-quadric product matrix, with every
whole column homogenized to its fixed QQ degree seven and target row shifts
zero. The 112 vectors reconstructed from three finite fields remain candidates.
Exactly forty-four of them have homogeneous degree seven.

All forty-four degree-seven vectors belong to the original homogeneous QQ
module. The saved constant rational matrix X of shape 56-by-44 satisfies
M^h X = S_7 exactly. This is a literal polynomial identity, not a modular
membership inference. The other sixty-eight candidate vectors have not been
proved to belong to the QQ module.

Expand each source and candidate vector in all 120 monomials of degree seven
in (a,e,P,w), together with its target component. This gives the complete
coefficient matrices C of shape 2880-by-56 and B of shape 2880-by-44. A modular
pivot search selects fifty-six rows of C. Exact rational Dixon solving on that
56-by-56 square produces X, after which full exact rational multiplication
verifies ALL 2880-by-44 entries of C X = B. Solving plus verification takes
0.011 seconds. The complete constant witness is 7.4 KB.

The prime is used only to choose an invertible square. The final exact QQ
polynomial identities provide the membership proof. These forty-four columns
are therefore certified redundant generators and may precondition a later QQ
GB. They do not prove a complete QQ GB, the correctness of the other candidates,
localized surjectivity, original product rank41, full cubic rank57, or an ambient
result. An augmented GB using all112 candidates still needs the complete QQ
GB, original containment, p-integral monic and good-prime homogeneous equality
checks and the graded Hilbert-function bridge. This packet awaits independent
review.

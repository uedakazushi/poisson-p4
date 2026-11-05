# Actual C2 second Betti number: independent source/minor review

The three accepted zero action matrices are exactly zA2,zA3,zA6. An
independent differential constructor uses the other five actions and
the ell action to reconstruct the active degree1,2,3 constant matrices.
It also reconstructs the full matrices and their literal signed
zero-wedge decomposition. The active224-by364 matrix has447 nonzeros.

The checker replays the provided126 independent-column operations,
checks all364 original column factorizations, and reconstructs the
selected original sparse126 minor entry by entry. Each normalized
pivot column has leading coefficient one. The source-change matrix is
triangular in the selected-column order. This proves the saved signed
product of nonzero rational normalizers without expanding it into one
large coefficient. No new elimination or rank search is performed.

The full degree-three rank is126+3*56+3*11=327. Consequently the actual
second Betti number over the original ell² base is518−89−327=102.
No T2 length, DVR torsion exponent or H0 calculation follows from this
rank alone. Transport of the active source units is a separate packet.

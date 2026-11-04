# Independent raw-H generic-rank and minor review

The independently accepted RAW_H_LINEAR_UNITS review binds every original
raw57 H coefficient, its physical frame, both known-unit directions and the
full saved 15-by6 checkpoint. This audit reuses that exact source gate and
reconstructs the new small rank certificate without a unit scan or solver.

For the saved polynomial p=M[2,0]=M[4,0], the literal rows 2,4 and columns
0,5 have determinant -(6/5) B p². The other three columns vanish in these
rows. Their complete support is exactly rows 7,9,10,11,12,13,14, whose
unchanged 7-by3 block has 908 terms. Clearing the two pair columns therefore
does not change this three-column block over the generic field.

The three already selected fraction-free steps are replayed. Each saved
nonzero pivot is independently checked as a determinant of the original
small matrix, and every quotient is checked by exact remultiplication. The
entire final checkpoint is compared. The selected 3-by3 determinant Delta
is independently computed by its six-term formula and has 407 terms. Its
adjugate reconstructs all seven row coordinates, and all 21 cleared upper
identities are checked. This proves rank3 for the small block. The actual
14-by5 minor is -(6/5) B p² Delta, with 1752 terms; its block-triangular
source binding is verified. With the two accepted known-unit directions,
the original actual raw16-by7 H-linear matrix has generic rank7.

No p or Delta is called a unit of all FULLU. Their inverses are used only
in the generic field and, separately, on the explicit principal open
described in the minor-open companion. No whole-FULLU rank7, fixed additional
global unit, full239, point, prime, GB or new kernel direction is claimed.

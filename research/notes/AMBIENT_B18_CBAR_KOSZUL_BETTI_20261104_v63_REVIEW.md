# Actual Cbar8 Koszul homology: independent review

The finite Koszul ranks and resulting minimal Betti numbers are accepted. The
source is the independently verified actual Cbar=C/(ell), of dimension8 over
QQ, with its eight coordinate actions in the SAME CA-F9 frame. They are
zA_0,...,zA_6,r_CA; after ell=0 the ninth s_CA action equals r_CA. Neither
raw physical rho coordinates nor the own CB parameter r_CB is substituted.

The independent checker reconstructs every differential from those actual
action matrices. Wedges are ordered lexicographically and coefficients by
the fixed8 module basis. For each wedge, deleting its p-th parameter applies
that coordinate action with sign(−1)^p. This reconstructs every sparse source
column in chain dimensions

    8,64,224,448,560,448,224,64,8.

There are1024 nonzero differential coefficients in total. Every action
commutator and every column of d² is checked in exact Fraction arithmetic.
No producer routine, new source value, jet or solver is executed.

Every saved elimination receipt is replayed against the independently
reconstructed source column. Each prescribed subtraction clears the claimed
pivot. Each declared zero residual is literally zero; each new pivot has its
declared nonzero coefficient, distinct leading row and normalized vector.
Thus the prescribed pivots are independent and span every source column,
proving both rank bounds. No new pivot or rank is searched. The ranks are

    5,35,105,175,175,105,35,5.

Subtracting adjacent ranks from chain dimensions gives the exact homology
dimensions

    3,24,84,168,210,168,84,24,3.

The coefficient ring Rbar=R0/(ell) is regular of dimension8 because the actual
CA-F9 functions are coordinates and ell=5(s_CA−r_CA) is regular. The Koszul
complex on these parameters resolves the residue field. Its tensor with
actual Cbar therefore computes Tor over Rbar. These numbers are precisely
the ranks of a minimal free Rbar-resolution; no explicit minimal resolution
matrices are claimed here.

The checker also reconstructs the inputs to every supplied maximal-ideal
filtration receipt and the stacked-action socle receipt. Their exact ranks
give filtration dimensions(8,5,2,0), minimal generator count3, Loewy Hilbert
function(3,3,2), and socle dimension3. Cbar is treated as a module, without
an asserted unital algebra structure.

No conclusion about T or K torsion, generic DVR rank, comparison kernels,
socle of the original intersection or D0 depth follows from this calculation.
The independent run is sparse and bounded by60s/2GB/32MiB/150000 coefficients.
All frozen primary files remain unchanged.

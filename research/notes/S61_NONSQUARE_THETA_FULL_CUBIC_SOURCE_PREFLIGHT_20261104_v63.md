# Full cubic57 on the two Theta curves: source and cost preflight

This packet only reads and counts immutable actual coefficients. It does
not substitute a curve, compute a rank, choose a new minor, evaluate a
POWER cofactor, or use a point or prime.

The exact source is `ELL_SIX_BORDERS/frozen57-source.txt`, with the accepted
row-basis hash. It contains the full actual 57-by-91 weighted cubic matrix:
41 pure product rows and16 genuine rows, pure84 plus heavy7 columns.
Every pure tail is retained. The frozen product41 theorem applies to all
FULLU; the reviewed Theta result gives heavy rank exactly6 on each P80 or
P100 curve. Neither uses the old extra Delta34 chart on these curves.

Over QQ(zeta) the exact remaining full-source test is therefore10-by-43.
Choose a nonzero pure41 product minor and six genuine rows with independent
heavy image. They give47 independent rows. Heavy rank6 makes the remaining
ten row classes purely light, in pure84/product41 of dimension43. Thus

    rank(full cubic matrix)=47+rank(ten-by-43 pure residual).

Full57 is equivalent to residual rank10. This is a field criterion on
each symbolic curve. A selected product minor is a tracked chart factor,
not automatically a whole-curve unit. The residual entries can be literal
48-by-48 bordered minors divided by the47 pivot determinant. For execution
one can reduce the original polynomial57-by-91 matrix directly; it avoids
constructing large multivariate completion cofactors.

A smaller fully serialized sufficient test is the accepted eight-unit
checkpoint:17-by-24 rational functions of a,e,d. Its rank17 implies rank25
for the actual projected ZW+psi matrix and hence full57. A smaller rank
does not refute full57: omitted pure directions may restore it. The formal
10-by-16 projected remainder after the later three quadratic and four
heavy pivots is not an evaluated small input; its old completion depends
on the degree207 power-cofactor circuit and is not proposed for expansion.

The coefficient file records exact input nonzero/term counts, degrees and
rational coefficient bit sizes for the full matrix, its product/genuine
pieces, the augmented25-by-32 source, and the17-by-24 checkpoint. It also
gives conservative bounds for a later one-variable substitution, without
performing it. For both curves, a and e have linear numerators divided by
zeta^2 and zeta respectively, and d is quadratic. Thus a^i e^j d^k has
Laurent support between k-2i-j and2k-i. Only known FULLU units are involved
in the existing source transformations and denominators.

An approved next job should specify one of these two exact tests and a
separate hard time/memory/output cap. A positive symbolic rank certificate
must retain its nonzero determinant polynomial and the actual source and
target indices. Uniform full57 on a whole curve would additionally require
removing all possible roots of maximal-minor factors in the permitted
open, rather than extrapolating a generic result. Full quartic239 remains
an independent gate.

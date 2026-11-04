# The degree-two Q10 coefficient-map kernel is zero

The exact right-kernel test for the fixed total-degree-two polynomial
class is complete. The QQ coefficient matrix has 347 rows, 260
columns and rank260. Its right kernel is therefore zero. The fixed
job took 2.536 seconds. It does not include a genuine source row,
a unit target, a new point, a GB or a degree-three test.

The polynomial input is the frozen `Q10x26` matrix with 317 core
coefficient terms. The ordered functional basis is exactly

`1,a,e,d,a²,ae,ad,e²,ed,d²`.

Each of the 26 target-functional components has ten rational unknowns.
Equating every coefficient of all ten equations `Q10*w=0` gives
347 rows. Their per-equation counts are

`[19,19,19,19,45,54,54,54,19,45]`.

The expanded coefficient input has 3,170 terms and maximum output
parameter degree six. All arithmetic is rational coefficient
linear algebra; no parameter-field rank or specialization is used.

The first 104 affine Q pivots and their exact source witnesses
were projected from the preceding saved functional test to its
quadratic equation rows. Their reconstruction was verified against
the current degree-two columns before use. Only the 156 new
quadratic-monomial columns were subsequently reduced. Every new
independent basis vector has its full exact witness in the input
columns. No column was dependent. The certificate saves all 260
normalized basis vectors, their pivots, the complete source witness
matrix and the coefficient row labels. These give an exact full
column-rank certificate.

The right-kernel basis is explicitly empty. The routine also contains
all-row and all-polynomial checks for any returned kernel vector;
those vector checks are vacuous here because the dimension is zero.
The rank certificate itself is not vacuous: all 260 independent
columns and their witnesses are saved and verified.

Together with the previous degree-zero and degree-one results, the
dimensions of the spaces of polynomial functional tuples of total
degree at most 0,1,2 annihilating Q10 are `0,0,0`. These are filtered
degree statements. The parameter coefficients of Q10 are not
homogeneous, so this is not presented as a homogeneous Hilbert-series
computation or a claim about module generators in all degrees.

In particular a next nonzero polynomial functional, if one is to be
constructed this way, must have degree at least three. This does not
exclude the known parameter-field annihilator or functionals of
higher polynomial degree. A common unit denominator with a degree-two
numerator cannot change the homogeneous annihilation equation after
clearing that denominator, but no broader denominator/degree class
is tested. No higher-degree solve is automatically launched.

The accepted global `14x21` presentation, full source and target
frames, heavy block, cubic57 and quartic239 frontiers, and all
uniform kernel counts are unchanged. The result is a finite-class
negative certificate, not a global rank or exhaustion conclusion.

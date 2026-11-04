# One symbolic P80 cubic test with no power-cofactor expansion

The exact job is positive: the remaining matrix has rank17 over QQ(zeta),
so the actual full cubic source has generic rank57 on P80. The saved
integer matrix has5,182 terms, degree34 and169-bit coefficients. Its
selected17-by-17 determinant is nonzero, degree447, with152 terms and at
most1,714-bit coefficients. The combined worker took5.393 seconds and
about51MB. The original input, all408 denominator classifications, all row
constant clearings, substituted matrix, pivot indices and determinant
coefficients are retained. Independent review is pending.

This is the one authorized combined substitution and generic-rank job for
the actual17-by-24 eight-unit checkpoint, on P80 only. It uses the frozen
source coefficients and accepted unit-operation/source-binding receipt.
No extra parameter point, prime, P100 run or full57 expansion is used.

Every input denominator is checked literally to be a nonzero QQ constant,
that constant times a, or that constant times B. On P80 use

    a=-8(7zeta+80)/(45zeta^2),
    e=2(zeta+40)/(9zeta), d=-3zeta(zeta+20)/80.

The single common scalar zeta^24*a*B is a certified unit on the whole
permitted P80 open. Multiplying by it removes all parameter denominators.
The exact coefficient substitution uses integer polynomial powers and
ordinary coefficients. Further row clearings are only nonzero rational
constants; their integer LCMs and removed integer contents are saved.
The resulting polynomial matrix is written before any rank operation.

The input support gives degree at most34 and at most14,280 ordinary
coefficient slots, below the150,000-term/32-MiB guard. The worker uses one
exact fraction-free17-by-24 elimination over QQ(zeta). It saves its pivot
rows/columns, the exact selected nonzero minor and its fraction-free
denominator. The selected minor is directly remultiplied by a separate
determinant evaluation on those same rows and columns; no adjugate is
expanded. A hard combined60s/2GB cap applies, with the source matrix
checkpoint retained if the rank stage fails.

The accepted checkpoint has removed eight genuine known-unit directions.
Thus its remaining rank17 implies rank25 for the actual projected ZW+psi
matrix and generic full cubic57 on P80. This implication is valid even
though the full raw-H projection has rank6 on the curve: the projected
test retains additional pure directions. A lower projected rank remains
inconclusive for the full cubic source. No uniform whole-curve result
follows merely from one nonzero polynomial minor, and quartic239 is not
decided by this test.

The possible full-cubic defect on the permitted P80 curve is contained in
the zero set of this saved nonzero determinant, hence is finite. The job
does not remove all roots, prove uniform rank57, or test P100. No further
rank job or factor stripping is triggered automatically.

# One symbolic P100 cubic test without power cofactors

The exact result is positive: remaining rank17, actual augmented rank25,
and hence generic full cubic57 on P100. The substituted integer matrix
has5,182 terms, degree34 and153-bit coefficients. The saved selected
17-by-17 minor is nonzero, degree447, with152 terms and at most1,471-bit
coefficients. The combined worker took5.659 seconds and about51MB;
independent review is pending. All coefficient/clearing evidence is saved.

This is the separately authorized one P100 substitution and generic-rank
job for the same frozen17-by-24 eight-unit checkpoint. The original source
rows, target columns and all408 denominators are retained. Every denominator
is checked as a nonzero rational constant times1, a or B.

The actual P100 chart is

    a=-10(zeta+20)/(9zeta^2),
    e=(zeta+100)/(9zeta), d=-3zeta(zeta+20)/100,
    B=zeta(zeta+60)/12.

The common scalar zeta^24*a*B is a known unit on the permitted curve. It
clears all parameter denominators; remaining row clearings remove only
nonzero integer/rational constants. The exact polynomial matrix and all
source and row clearings are saved before fraction-free elimination.
The support guard is degree34 and at most14,280 ordinary coefficient slots,
under the150,000-term and32-MiB guards. The combined cap is60s/2GB.

One exact fraction-free rank computation over QQ(zeta) supplies the chosen
rows/columns and their exact nonzero determinant, directly checked on that
selected submatrix. Rank17 implies actual augmented rank25 and generic
full cubic57 on P100 through the accepted source projection theorem. A
smaller rank would leave full57 open. A nonzero determinant alone gives
only generic rank, not uniform rank at all permitted zeta values.

There is no new source or beta point, prime, P80 repetition, degree207
cofactor expansion, full57 coefficient expansion, or automatic follow-up
minor search. Full quartic239 remains a separate problem.

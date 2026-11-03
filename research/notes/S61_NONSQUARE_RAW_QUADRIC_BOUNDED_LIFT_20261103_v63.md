# Raw four-quadric module: exact obstruction to degree-at-most-three cofactors

Let M be the literal QQ 24-by-56 matrix frozen in
`research/handoff/20261103_v62/v36_raw_fourquad_GF_pilot`: all six light-linear
products of the four cleared actual quadrics, after d=P/a and multiplication by
row powers a^[0,2,2,1]. Its coefficients lie in QQ[a,e,P], have total degree at
most seven, and all denominators are units at p=32003. Write J=5P−6e−4 and
h=a^5 J^3, of degree eight.

For every nonzero constant vector c in QQ^24 there is no polynomial vector
x in QQ[a,e,P]^56, with each entry of total degree at most three, satisfying
M x = h c. In particular no one of the twenty-four coordinate targets h e_i
has such a cofactor vector. This excludes only this finite coefficient-degree
class. It does not exclude higher-degree QQ inverses, and proves no pointwise
rank defect or failure of the desired QQ localized surjectivity.

For D=1,2,3, expand every cofactor in all monomials of total degree at most D.
The coefficient matrix C_D has 56*binomial(D+3,3) columns and
24*binomial(D+10,3) rows, indexed by output component and all core monomials of
degree at most D+7. Augment it by the twenty-four coefficient vectors of h e_i.
All source coefficients and targets are p-integral. Exact FLINT field RREF gives:

| D | rows | columns of C_D | rank C_D mod p | augmented rank mod p |
| --- | ---: | ---: | ---: | ---: |
| 1 | 3960 | 224 | 224 | 248 |
| 2 | 5280 | 560 | 560 | 584 |
| 3 | 6864 | 1120 | 1120 | 1144 |

Each entire augmented matrix has full COLUMN rank modulo p. Its columns are
therefore linearly independent over QQ: a nonzero full-column minor modulo p
is a p-integral minor nonzero over QQ. This proves the assertion, including
independence of all twenty-four targets modulo the bounded-degree image. No
QQ kernel reconstruction, GB, or specialization of a localized unit assertion
is used in this rank argument.

The complete coefficient matrices are regenerated directly from every literal
source polynomial, not from random evaluations or selected monomial samples.
The largest test takes 1.79 seconds with an RSS upper bound of 103660 KiB.
Fresh positive and negative 24-by-56 fixtures use the diagonal polynomial
`a^5 J^2`: the positive fixture has all twenty-four solutions with cofactor J;
removing its first diagonal entry leaves exactly twenty-three solvable targets.
Both fixtures pass. The source, fixtures, native arithmetic program and receipts
are retained. This packet has not yet received independent mathematical review.

# Independent source binding of all 670 five-functional entries

The selected-five formal-curve review proved rank5 by independently binding
25 entries. It did not bind the other 645 stored entries. This packet closes
that precise source-binding gap over QQ[H,L,s0,q0,t0], rather than extrapolating
from the one-variable curve.

The existing complete 134-by5 polynomial matrix has 692953 terms, at most
1685 per final entry and 7965 per row; the largest polynomial intermediate
has 2083 terms. Its original constructor subsequently stopped in an unrelated
field-elimination guard. That failure packet remains unchanged; no field
elimination or rank attempt is retried here. The larger total term count is
explicitly authorized for an independent streaming polynomial comparison.

The new checker is not the author's builder and imports no producer code.
For each original equation F it reads every original q monomial and its
five-variable coefficient, and directly constructs

    [P F(C), J_F(C) W0, J_F(C) W1,
     polar_QF(W0,d), polar_QF(W1,d)].

Here C, P, original W0/W1 and d=Q52-C are the independently accepted sparse5
source inputs. Linear terms, constant terms and both contributions of a
diagonal square to each polarization are retained. The derivative entries
are assembled term by term from polarized products, rather than reusing the
producer's collected J(C) row. Each of the five actual QQ polynomials is
compared exactly with the corresponding stored polynomial; the row is then
discarded. There is no rational-function arithmetic, GCD, row scaling, pivot,
sample point or finite prime in this check.

The existing matrix is streamed directly to the checker. No duplicate full
matrix or expanded output artifact is saved. Reports contain only preflight,
counts, hashes and the exact-comparison receipt; the small original-source
input and executable stay in ignored work/. The report's input hashes
explicitly bind the existing five-column-matrix.json.gz and the original134,
center and direction data.

All 670 QQ polynomial equalities pass, binding every stored functional to
its actual original equation. Combining this source gate with the accepted
rank5/functional-injectivity result makes the five-functional row space a
faithful representation of all 134 plane restrictions. Plane-zero of the
129 proposed residual equations also requires their separate exact
adjugate/row-combination circuit proof; that proof is not substituted by the
present coefficient check. No full affine-fiber exhaustion is claimed.

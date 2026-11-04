# Independent raw H-linear2 and pointwise rank126 review

This review checks both frozen primary packets
`S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63` and
`S61_NONSQUARE_H_LINEAR2_RANK126_20261104_v63` over
`FULLU=QQ[a,e,d,1/(adFBCJ)]`. It does not run the author's solver,
scan for a pivot, choose a point or run a GB.

The reviewer directly parses all 112 H-times-light coefficients of
the last sixteen rows of the immutable original `57x91` source.
It independently applies `v6=Z-3*x0+x1`, `v3=W` to the seven
coefficients and checks both constant-frame inverses. It rebuilds
`Cdelta=-64000*C^2*(C6-C1)`, its literal two-coefficient H form,
the Bezout identity, eta and gamma. All entries of the GL2 inverse
products are checked over QQ(a,e,d). The full-source subtraction
coefficients are saved; their pure parts are retained as part of
the source operation, not discarded in the H projection.

The reviewer reconstructs every coefficient of the polynomial
`15x6` remainder, compares all 90 entries to both saved formats,
and verifies the sole fixed completed step. Initial row divisors,
each saved pre-operation coefficient, every exact output divisor
and all 90 final checkpoint entries are checked. The cleaned pivot
is exactly `-B*C/14062500`; its original selected entry is exactly
`F*B^2*C/12500000`. Every factor divided in these operations is a
named FULLU unit. No new unit scan or further Gaussian step occurs.

The two actual H-linear source directions are therefore independent
at every geometric point of FULLU. This is a uniform rank lower
bound two, not full H7 generation. The native audit's field named
`unfinished_operation` contains the last completed-step label:
the author variable is not reset after a completed step. The exact
fixed transcript and final matrix verify that there is no partial
arithmetic update in the saved checkpoint. This metadata detail
does not add or remove a completed pivot.

For the rank126 corollary, work over any residue field of a FULLU
point. Two independent light linear forms can be made the first
two coordinates. Their seven light multiples span the 13 distinct
quadratic monomials containing at least one of these coordinates;
the other five coordinates have 15 quadratic monomials, so the
dimension is `28-15=13`. The accepted exact NF117 projection
preserves H degree and removes exactly five of the 28 H-quadratic
directions. Its image on this 13-dimensional subspace consequently
has dimension at least eight. This does not require a computed
basis of that image or of a projective complement.

All genuine cubic light multiples have H exponent at most one,
by their actual weighted source degrees. The separate accepted
quartic H-squared unit direction is therefore independent of this
cubic image after the NF117 elimination. The full actual source
matrix has pointwise rank at least `117+8+1=126` everywhere on
FULLU. The proof uses every residue field, not specialization of
a single sample.

The accepted result is a pointwise lower bound 126. It does not
supply a particular global unit 126-minor, an eighth constructive
H pivot, a free target complement or an automatically reduced
source basis. Whole cubic57, quartic239 and the known kernel38
with gap42 remain open/unchanged.

# Five fixed9-minors: exact coefficients and an inconclusive modular cover

Date: 2026-11-04. Exact coefficient primary; independent review pending.
The whole-P100149 lower bound and typed pure90 presentation remain valid.
This packet does **not** prove a whole-curve150 bound or a rank failure.

## Source and exact determinants

Append the actual pure source K2*v3 (row11 of the saved15-by98 matrix)
to each of the five weighted8-source selectors from the accepted149 cover.
Keep their eight targets and append old local98 target25. Each9-by9 matrix
is made polynomial by its known five-factor curve-unit row LCDs. All405
source clearing identities are checked; no extra old-row unit scaling is
needed. Source coefficients are cached whole-curve RF values, with no
new source or normal-form query.

The five determinants are computed exactly over QQ[z] using row-pivot
Bareiss. Every quotient is verified by polynomial multiplication,1020 exact
division identities in total. Each case saves its row swaps, pivots,
previous pivots, final-node sign, and a204-node exact-quotient replay DAG.
The resulting determinant degrees are `[214,209,198,212,215]`.

Strip only the five known curve units and a nonzero QQ content. The remaining
positive-leading primitive ZZ polynomials have degrees `[81,78,69,88,88]`
and82,79,70,89,89 terms. The full factor/content reconstruction is verified.
No unknown factor is removed or inverted.

## Existing-prime polynomial certificate

Every QQ coefficient denominator encountered in the input, exact quotient
and normalization checks is prime to32003;183428 coefficient checks pass.
Each nonzero primitive polynomial has leading coefficient nonzero modulo
32003. Thus all five reductions preserve their degrees.

The exact gcd of the five **coefficient polynomials** over F_32003 has
degree15, with successive degrees `[81,74,24,15,15]`. The saved modular
Bezout identity is

```
sum_i b_i(z)*Z_i_bar(z)=G_bar15(z),
```

not1. It is checked literally. No numerical zeta, beta or source point is
chosen, and no alternate prime is tried. The degree15 modular gcd alone
does not imply a nonconstant QQ gcd, a common complex root, an actual rank
drop, or a geometric component.

Because all primitive leading coefficients are32003-units, Gauss gives a
useful upper bound: a common monic QQ factor is integral over Z_(32003),
remains monic of the same positive degree on reduction, and divides G_bar15.
Hence the actual QQ gcd has degree **at most15**, possibly0. The common
zero set of this particular five-minor family in the P100 chart is therefore
an unresolved set of at most15 geometric points, possibly empty. This
concerns the chosen minor cover only; other minors or sources can rescue it.
No QQ polynomial gcd or QQ Bezout has been computed.

## Scope and saved frontier

The proposed whole150 gate is inconclusive. Away from the common zeros of
the saved five QQ primitive polynomials, one actual9-minor is nonzero, but
the earlier geometric generic result is not being reclassified as new.
Whole149 remains the accepted explicit uniform bound. No pure89/global9
right inverse, full239, extra target/source, two-row replacement, Laurent
series, new point, prime or GB is included.

Code: `code/s61_nonsquare_P100_next9_cover_20261104_v63.py` and `.cpp`.
Reports: `S61_NONSQUARE_P100_NEXT9_COVER_20261104_v63`. `native.json` saves
all five determinants, unit/content factorizations, Bareiss transcripts,
primitive polynomials, their degree-preserving reductions and the modular
gcd/Bezout. The frozen source preflight is a separate pinned dependency.

The one job took5.510 seconds including compilation, native0.697 seconds.
The peak was75226 polynomial terms including a51034-term literal reserve;
10629 terms were serialized. The combined source and result packet is
2157836 bytes, below32MiB. No cap was raised or job repeated.

Report SHA-256:
`5670647ac97cb82816852b21b50e70783d17c114e478a7cfdfd95b2553cba63c`.
Native output SHA-256:
`39fb05e69f71f3270fb0627e5d00ca5230104fd77c7605d2d2fcac75b12846b7`.
All primary bytes are immutable. The finite chosen-cover frontier is saved;
there is no automatic next calculation.

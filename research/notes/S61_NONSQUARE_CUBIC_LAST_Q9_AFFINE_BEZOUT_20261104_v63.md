# Fixed constant and affine Bezout coefficient classes

This packet tests two finite coefficient classes for a Bezout
identity on the nine actual last-row coefficients saved in
`S61_NONSQUARE_CUBIC_LAST_Q9_PREFLIGHT_20261104_v63`.
It performs exact rational coefficient linear algebra, not a
parameter-field rank test or a search for new unit entries.

Let the frozen nine polynomials be `f0,...,f8`. The first eight have
the literal recorded unit factor `a`; the last does not. The source
basis is therefore `g_i=f_i/a` for `i<8`, `g8=f8`. Every factor is
verified coefficient by coefficient; no other polynomial factor is
removed or inverted.

The first stage allows nine constant rational cofactors. The second
stage extends precisely this basis by multipliers `a,e,d`, giving
36 rational unknowns. The eight fixed target units are

`1, a, d, F, B, C, J, a*d`.

The union of coefficient monomials has 117 rows. Source and target
support together has 1,513 serialized terms, below the 10,000-term
preflight bound. The exact QQ span ranks are respectively 9 and 36.
None of the eight targets belongs to either span. The calculation
took 0.113 seconds.

The certificate saves the common monomial basis, all normalized
source polynomials, both exact span bases with witnesses in the
fixed source columns, and every reduced target remainder. It
verifies that each saved basis vector equals its source witness
and that every tested source column reduces to zero. A nonzero
target remainder is therefore an exact obstruction in this
particular finite coefficient class. If a candidate had existed,
the code would also have verified the identity against all original
QQ coefficients after multiplying through the known `a` content;
there were no successful candidates in this run.

The finite classes correspond to actual source cofactors
`p_i/a` on the first eight source entries and `p8` on the last,
where each `p_i` is constant or affine in `a,e,d`. Failure does not
exclude higher-degree cofactors, other regular FULLU coefficients,
different target unit powers, or another explicit completion. It
does not prove that the ideal is proper or that the remaining row
vanishes anywhere. In particular there is no conflict with the
accepted pointwise rank-one result for the remaining full row.

All pure and heavy tails stay in that full actual row. No higher
degree, new target, GB, saturation, point, unknown denominator or
automatic target/source transport was tried. Full cubic57,
quartic239 and kernel38/gap42 are unchanged.

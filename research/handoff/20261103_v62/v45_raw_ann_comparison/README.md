# Raw annihilator candidates at the same three fields; degree-four test

No new prime is used. The two existing large fields2147483647 and2147483629
recompute the raw24-by56 standard basis and its canonical annihilator sequentially
in12.46/13.91seconds with RSS about68MiB, under60-second /2200000-KiB caps. Their
annihilator leading monomials match the32003 field saved in v36. All polynomials
are aligned by leading monomial and normalized to be monic before CRT.

At the same77-bit three-field modulus, all coefficients of the first six of seven
annihilator generators agree under adjacent-prefix rational reconstruction.
The seventh has only16/28 stable coefficients. In particular the two small
QQ CANDIDATE polynomials are

A2=a^5 J²(8P+15)/288,
A3=a^5 J²(80P+48e+257)/1728,
J=5P−6e−4.

They satisfy the literal QQ polynomial identity540*A2−216*A3=a^5 J³.
Coefficient reconstruction does NOT prove that either polynomial annihilates
the QQ cokernel. All seven native reconstructed polynomials remain explicitly
marked candidates, including the unstable seventh one.

A single bounded degree4 coefficient test appends all48 coordinate targets
A2*e_i,A3*e_i to the same complete literal source coefficient matrix, using1960
unknown coefficients per target and8736 output rows. Exact field RREF gives
full source column rank1960 and full augmented column rank2008. Thus no nonzero
constant QQ combination of those48 explicit target polynomials is in the
cofactor-degree-at-most4 image. This is a rank-minor argument over the p-integral
literal QQ matrix; it excludes only that finite degree class, not higher lifts
or QQ annihilator membership. It takes6.64seconds /225024KiB. No degree5 test,
QQ solve or additional prime is included.

The exploratory stdout's key `rank_augmented_all24` is inherited from the earlier
program; in THIS run it describes all48 RHS columns, as explicitly implemented
and indicated by `all48_solvable`, the48 individual statuses and total2008.
The original code/stdout are retained without overwriting that provenance.

Full field basis/annihilator/saturation data, source pins, runtime receipts,
CRT/candidate scripts and the bounded test are saved. No QQ rank41, complete QQ
GB, localized unit or ambient result is inferred from these candidates.

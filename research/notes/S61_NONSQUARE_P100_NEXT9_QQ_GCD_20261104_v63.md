# Exact finite chosen-minor frontier on P100

Date: 2026-11-04. Exact gcd primary; independent review pending.
Only the five already saved primitive ZZ coefficient polynomials are used:
degrees81,78,69,88,88, with409 integer coefficients. No determinant, source
or normal-form coefficient is recomputed.

Their exact monic QQ[z] gcd H has **degree15**. Its positive-leading primitive
integer form H_Z satisfies `H=H_Z/5585717947953`. All five literal divisions
`Z_i=H*Q_i` are saved and checked. H and every Q_i have32003-integral QQ
coefficients; the quotient leading coefficients are32003-units and their
reductions preserve degree. The saved old-prime modular identity is

```
sum_i b_i*Q_i_bar=1 over F_32003[z].
```

It is checked literally. Gauss then proves that the Q_i have no common
nonconstant QQ factor. Together with the five exact divisions this certifies
H as the actual gcd, rather than merely the modular degree upper bound.
No QQ Bezout coefficients are expanded. The integer primitive-PRS fallback
was prepared but was not needed or executed.

The actual chosen9-minors differ from the Z_i only by nonzero QQ contents
and products of the five known curve units. Thus their ideal in T100 is
exactly `(H)`. This is an explicit finite **chosen-cover** scheme, not a
Poisson branch or a certificate of actual quartic rank failure. The whole
149 lower bound remains accepted; whole150 is not established.

The algebra `QQ[z]/(H)` has dimension15 and is not assumed to be a field,
irreducible or reduced. H has no known excluded linear factor because it
divides an input already stripped completely of those factors. Since the
monic coefficients are32003-integral, `Z_(32003)[z]/(H)` is a finite free
rank15 model. A future finite-fiber rescue must use its entire special fiber,
including any nilpotents, and cannot infer a global QQ result from one open
or sample. No extra minor, inverse or finite-fiber computation occurs here.

Code: `code/s61_nonsquare_P100_next9_QQ_gcd_20261104_v63.py` and `.cpp`.
Reports: `S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63`, containing the monic
and primitive gcd, all five quotient identities, their modular reductions
and cofactor Bezout. The early `gcd-phase.json` is a preserved checkpoint;
the final `native.json` supplies the completed coprimality certificate.

The one job took4.371 seconds including compilation, native0.003 seconds.
Peak11902 polynomial terms including the11038-term literal reserve;
716 terms were serialized. No budget was raised or job repeated. There was
no QQ Bezout, factor search, new minor, target/source, point or prime.
All primary bytes are immutable.

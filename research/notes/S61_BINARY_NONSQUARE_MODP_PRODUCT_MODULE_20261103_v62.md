# Exact finite-field right inverse of the pinned10-row residual

2026-11-03. This is a complete GF32003 polynomial certificate for a **reduced
QQ residual matrix**. It is not a characteristic0 rank theorem and does not
claim original product rank41 in characteristic32003 without the extra
prime-compatibility audit below.

The source is the first ten product rows of the acceptedV9 full26×60 QQ
residual, SHA256
`6d21ebfbeade3c5a0a7bd2dc091bde5ec274622ddf899add905242164d616ea3`,
after31 QQ product-unit pivots. Omit only its seven identically zero product
columns, leaving10×53. Reduce every coefficient modulo32003; all source
coefficient denominators are prime units. Here the third serialized variable
d means P=ad, B0=aB, C=100a+18e²−9e−2 and J=5P−6e−4.

The saved polynomial53×10 matrix L satisfies

    M L = a P³ C J² identity10

in ALL100 entries over GF32003[a,e,P]. L has183182 terms and total degree at
most18. A separate FLINT `nmod_mpoly` implementation reads the complete native
matrices and verifies every polynomial coefficient, independently of Singular's
GB and lift algorithms. It completes the100 checks in about0.27 seconds in the
current environment. The replay driver also compares all530 saved source
entries to direct reduction of the pinned QQ matrix. This is global polynomial
algebra, not evaluation at a set of finite-field points.

Consequently this reduced residual has full10-coordinate image over the named
open aPFB0CJ≠0, including all geometric points of that open in characteristic32003.
Discovery used Singular slimgb, quotient by the full free module, and saturation:
GB106 columns, annihilator10 generators and saturated annihilator1. Its first
annihilator generator is exactly aP³CJ²/648. The final acceptance needs none of
these GB claims: the explicit matrix identity supplies the right inverse.

## Exact scope and the missing bridges

To deduce ORIGINAL product41 over GF32003, one must additionally verify that
ALL earlier31 QQ pivots and row/column normalizations specialize invertibly at
32003. Merely knowing that the final matrix has prime-integral entries does not
prove this; dividing an earlier matrix by the prime can destroy that inference.
This audit is not claimed here.

To deduce product41 over QQ, one must produce a QQ right inverse (or another
valid exact/proper finite certificate). A localized finite-field unit ideal
alone does not imply the corresponding QQ statement; components can specialize
into the inverted divisor. The ongoing multi-prime cofactor reconstruction uses
finite-field data only as a search tool. Its acceptance requires a complete QQ
polynomial identity in all100 entries, followed by the already accepted QQ
31-pivot bridge. Full cubic57, quartic65 and ambient removal remain separate.

## Serialization failure retained and repaired

The first GF lift passed its in-memory identity, but a statistics variable named
`maxdeg` collided with a Singular library procedure. Moreover `write(file,L)`
silently wrote only the first matrix row. That47KB file is explicitly retained
as INVALID, with the failure log; reloading it failed the matrix identity. It
was never accepted as a QQ candidate or a full saved right inverse.

A fresh2×3 ASCII-link fixture established explicit per-entry export. The lift
was then rerun with unique variable names, all53×10 entries were serialized,
and the independent FLINT100-entry check passed. The complete native data and
failed partial data are preserved with hashes; no environment setup was repeated.

Replay, using existing activated FLINT/compiler tools:

    .venv/bin/python code/s61_binary_nonsquare_modp_product_module_20261103_v62.py --output work/modp-product-module.json

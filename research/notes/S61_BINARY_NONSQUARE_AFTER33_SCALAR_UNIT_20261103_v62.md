# Exact next scalar unit row after the33-pivot primary reduction

2026-11-03. For row7 of the actual24×58 QQ cubic residual after33 constructive
product pivots, a saved58-entry polynomial vector c satisfies

    row7 c = a P² J,    P=ad, J=5P−6e−4.

This is a declared unit on aPFB0CJ≠0. The complete original QQ identity is checked
by a separate FLINT polynomial dot-product implementation, and every saved source
row entry is compared to the pinned24×58 matrix. No claim rests on modular rank.
The row is one of the eight remaining product rows, not a genuine cubic row.

## Exact coefficient search and certificate

The search strips only declared unit factors from each entry of this SINGLE row
and normalizes nonzero QQ scalars. Polynomial cofactors of coefficient degree≤4
are allowed on those primitive entries:58 entries expand to2030 scalar columns.
A prime-field rank calculation selects one1008×1008 exact QQ system. The solver
uses Dixon over QQ, then verifies the FULL polynomial identity on the primitive
row, clears only named unit denominators, and verifies the FULL ORIGINAL row.
The final vector has975 terms in33 active components, total degree7 and maximum
3299 coefficient bits. Its unit factor exponents in `(a,F,B0,C,J,P)` are
`(1,0,0,0,1,2)`, with scalar1. The independent checker accepts only nonzero QQ
constants and these explicitly named factors.

The degree≤3 search found no candidate in its modular target-selection class.
This is a bounded search refusal, not a QQ nonmembership or arbitrary-degree
obstruction: a modular rank loss can conceal rational solutions. The degree4
positive result is accepted solely because the complete QQ identity holds.

## Reduction count and remaining work

The certificate provides a concrete next unit row. This packet deliberately
freezes the exact row, cofactors and independent checker BEFORE a growing Schur
expansion. The executed primary count remains33, with24×58 as its current saved
matrix. A redundant image column can always be appended and pivoted, giving
23×58 while preserving the cokernel over the named open. A23×57 replacement
requires a known-unit cofactor or another verified invertible column change;
neither is claimed here. A large generated expanded matrix is unnecessary for
this row-identity result and is not saved in this packet.

All-point product rank41, full cubic rank57, quartic65 and ambient exhaustion
remain open. This QQ scalar result is independent of the stalled raw eight-prime
cofactor reconstruction and of the accepted finite-field reduced-residual theorem.

Replay with the existing activated compiler/FLINT environment:

    .venv/bin/python code/s61_binary_nonsquare_after33_scalar_unit_20261103_v62.py --output work/after33-scalar-unit.json

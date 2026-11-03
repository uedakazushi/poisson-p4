# Frozen eight-prime cofactor frontier — no QQ rank claim

This preserves a bounded reconstruction experiment for a QQ right inverse of
actual pinned31-prefix product residual M10×53. The target is exactly
`aP³CJ² identity10`; P=ad, C=100a+18e²−9e−2, J=5P−6e−4. A final successful
QQ identity would combine with the accepted QQ31 pivots to prove original
product41. No such QQ matrix is certified by this archive.

All eight runs use the SAME explicit constructor:

    option(redSB); B=liftstd(M,T); H=lift(B,g*freemodule(10)); L=T*H.

Every finite-field run checks BOTH M*T=B and M*L=gI10 before serializing ALL
53×10 entries with ASCII link. The entire coefficient-free GBprotocol,
all GBleading monomials and Tleading/degree/termcount signatures match exactly
across eight requested31-bit primes. All124942 coefficient positions match;
cofactor degree is17. The primes are prime (also checked by SymPy), and a
separate native char check confirms large-prime arithmetic is not truncated
to a small field. Each completed run takes27–33sec and about159MiBRSS.

These structural matches are consistency checks, not a proof that arbitrary
modular right inverses reconstruct to a particular QQ solution. The original
32003 lift used a DIFFERENT constructor and must NOT be mixed into this CRT.

At three primes (93bits) and eight primes (248bits), rational reconstruction
has only the random baseline~60.8%coverage and ZERO coefficients stable across
the adjacent modulus prefixes. This remains true after source-column leading
coefficient/content scalings, and for asymmetric denominator bounds8…96bits.
No complete QQ candidate is available. Neither a coefficient-height bound nor
nonexistence of a QQ right inverse is inferred. Automatic prime extension stopped
at eight; any further batch should record a new method or measured height rationale.

Extract at repository root and reuse activation. Every field matrix is archived
with source, run, protocol and leading-signature receipts. The large generated
input `.sing` files are reproducible from `prepare_prime.py PRIME` and the pinned
QQ residual (restore V6/V9). Replay a prime with `run_prime.py PRIME`; at most two
single-process Singular jobs were used. `assess_three.py` and `assess_eight.py`
reconstruct the archived CRT states; caches are omitted with hashes because they
are quickly regenerated. The C++ assessment helpers use exact FLINT rational
reconstruction, but no reconstruction coefficient is a theorem until the entire
QQ polynomial identity is verified.

The separate complete32003 polynomial certificate IS accepted and independently
reviewed in S61_BINARY_NONSQUARE_MODP_PRODUCT_MODULE_20261103_v62. Its scope is
only the p-reduced QQ residual. Earlier31-pivot p-unit compatibility, originalGF
product41, QQ product41, full cubic57 and quartic/ambient coverage remain outside
that finite-field certificate.

Research continues separately: a canonical cofactor normalization modulo the
full kernel of M is being explored to avoid arbitrary lift/backtrace choices;
the first8-prime raw experiment here is immutable. Other scalar-unit attempts
and exact bounded-degree searches are not folded into this archive as QQ results.

# No-new-prime coefficient reconstruction comparison

This packet reuses complete field bases saved in v29-v34. It launches no GB or
new prime computation. The same three LARGE primes 2147483647, 2147483629 and
2147483587 give a fair 93-bit comparison between the earlier safe21 and later
31-pivot homogeneous source matrices. Leading monomials identify columns and
full nonleading support identifies coefficients; missing terms are zero-padded.
All source and order/homogenization checks remain those in the source packets.

Balanced adjacent-prefix agreement is 16692/531441 (3.14%) for safe21 and
10992/506068 (2.17%) for the later source. The earlier input has maximum QQ
coefficient heights 36/47 bits, degree 14 and 30731 terms; the later input has
217/222 bits, degree 26 and 205165 terms. Both still lack a full QQ candidate.

The additional bounded assessment pools denominators of candidate fractions
that agree between adjacent prefixes WITHIN each monic GB column. For each
column it takes their exact integer LCM D and compares signed integer lifts of
D times every coefficient residue at the same two prefixes. No denominator
above 4096 bits is processed; none exceeded that guard. The results are:

| source / primes | independently stable | pooled stable | observed complete columns | maximum D bits |
| --- | ---: | ---: | ---: | ---: |
| earlier safe21 / 3 | 16692 | 17943 | 29 / 1050 | 75 |
| later31 / 3 LARGE | 10992 | 11737 | 55 / 739 | 77 |
| later31 / 16 | 71024 | 85260 | 88 / 739 | 449 |

“Stable” and “observed complete” describe finite residue reconstruction only.
They do not certify the QQ fractions, that D is the true common denominator,
or that reconstructed vectors belong to the original QQ module. Every final
candidate would still need the full exact QQ GB and containment checks, p-unit
and homogeneous good-prime equality checks, and graded Hilbert-function bridge.
These exploratory data do not establish any QQ rank or membership theorem.

The compact packet saves the exact pooling source, generation script, per-column
receipts and fair-comparison generator/reports. Large reproducible CRT caches
and pair-input files are omitted with hashes and replay instructions. No ELF is
included. All prime computations remain frozen at the previously authorized
limits of sixteen for the later source and three for safe21.

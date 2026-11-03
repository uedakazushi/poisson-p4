# Frozen homogeneous-GB pilot — finite-field data only

One GF32003 computation of the complete column-homogenized10-row residual is
preserved here. This is a reconstruction frontier, not a QQ GB/module/rank theorem.

Input is the SAME pinned actual QQ31-prefix residual M10×53. Every ENTIRE source
column is homogenized to its fixed maximum QQ total degree, before reduction at
the prime. These degrees range10–26 and are recorded for all53 columns; all ten
row shifts are0. Variables are(a,e,P,w), serialized P as d. The exact order is
(dp,C); the native vector[a,a] smoke verifies that the higher component wins a
monomial tie. No degree drop at the prime changes the chosen QQ homogenization.

The single-process slimgb(redSB) run finished39.72sec/182MiBRSS under120sec/2.2GB.
The complete GB is10×739,506793 terms, maximum degree29. All columns are homogeneous,
all739 leading coefficients are1. A separate C++ data audit checks38254597 possible
tail/leading divisibilities and confirms the basis data are fully reduced. The
source generators reduce to zero in the computed finite-field basis. Native
ASCII-link export preserves every entry and all four exponents.

Field coefficient representatives use at most14 bits. QQ coefficient heights are
UNKNOWN: finite-field representative sizes do not bound them. The later bounded
three-prime canonical-GB coefficient assessment continues separately and does not
modify this single-pilot archive.

For a QQ result, every following step remains necessary:

1. Reconstruct the same monic homogeneous GB, with complete leading/support checks.
2. Verify its COMPLETE QQ GB property (all required S-pairs/standard-basis check).
3. Verify exact QQ containment of the original homogeneous columns Ih in G.
4. Verify p-integrality/unit leading coefficients and exact equality with the SAME
   good-prime homogeneous original module, in the SAME order and row shifts.
5. Apply the graded Hilbert-function sandwich to deduce Ih=G; only then dehomogenize
   and verify the desired localized rank criterion.

An inhomogeneous localized finite-field unit ideal alone is not a substitute for
these conditions. No reverse cofactor trace or originalGF31-pivot compatibility is
asserted by this pilot. The independent QQ33/34-pivot results are separate.

Extract at repository root, reuse activation, and run the saved source.sing with
`run-singular-optimized-v6.sh --cpus=0 -q`. The prepare script reconstructs all input
coefficients from the pinned QQ residual (restore V6/V9 if needed). The native GB,
source, fixed-degree metadata, runtime, order smoke, homogeneous/monic/reduced data
checks and hashes are archived. No ELF binary or large QQ source cache is included.

# Frozen canonical homogeneous GB sixteen-prime frontier — no QQ theorem

Extract after the immutable v29, v30 and v31 archives. The two new chunks
preserve all eight additional complete finite-field bases, source pins, execution
receipts, reduced-basis audits and reconstruction scripts. Each chunk is below
32 MiB; no executable binary is included. The large derived CRT cache is omitted
with its SHA and replay command in the manifest.

All sixteen computations use the same actual QQ residual 10-by-53 matrix,
homogenizing each entire column to its fixed QQ degree before reduction. Row
shifts are zero and the order is `(dp,C)`. Every basis is monic, homogeneous and
fully reduced, with 739 columns, degree at most 29 and the same leading module.
Large-prime supports have 506807 terms. The fourteen missing nonleading
coefficients at 32003 are explicitly zero-padded in the common 506068-position
nonleading support.

At sixteen primes the modulus has 480 bits. Balanced reconstruction agrees
between adjacent prefixes at 71024 nonleading positions (14.03%), compared with
31081 at eight primes and 42986 at twelve. Observed stable balanced heights reach
224 numerator and denominator bits. With a denominator bound of 128 bits,
50207 positions agree and observed numerator heights reach 320 bits, exceeding
the actual input numerator maximum of 217 bits. These are finite reconstruction
observations. They do not prove a height bound or provide a complete QQ basis.
The approximately 60% Farey coverage baseline is not counted as convergence.

This prime extension stops at sixteen. The next distinct bounded pilot uses the
earlier safe 21-pivot actual product residual, whose coefficients have only
36/47 numerator/denominator bits. Its 20-by-63 matrix expresses the same QQ
product-rank-41 target; it is not the stronger four-quad rank-24 criterion.
The source indices and measured comparison are preserved here, and its outcome
will be saved separately.

QQ acceptance remains unfinished: verify the full QQ GB property, exact original
homogeneous-module containment in the candidate, homogeneous columns with the
same shifts/order, p-integral monic coefficients, and exact equality at a good
prime with the same original homogeneous module. The graded Hilbert-function
sandwich then proves equality over QQ. Only after these steps can exact QQ normal
forms establish the desired localized membership. The saved complete finite-field
cofactors refine the concrete target to `w^22*g^h*e_i` in degree 30 for all ten
rows, where `g=a*P^3*C*J^2`; all ten QQ remainders still require verification.
No QQ rank-41, full cubic rank-57, quartic, or ambient theorem is claimed here.

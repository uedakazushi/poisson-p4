# One raw-coordinate finite-field lift, with independent complete checks

The single authorized GF(32003) liftstd run uses the literal raw24-by56 source
and ONLY the target h*e0, h=a^5*(5P−6e−4)^3. It completes in24.96seconds with
RSS440488KiB, under60-second /2200000-KiB caps. The saved56-by1 polynomial
cofactor has17990 terms and maximum degree18. No full24-target lift is computed.
All existing canonical bases remain unchanged.

The original script's internal subtraction compared a24-by1 matrix with a
rank1 target module and raised a matrix-size error. That failed check and stdout
are retained. The saved cofactor is accepted ONLY because a separate direct
FLINT polynomial multiplication against the entire literal QQ source reduced
modulo32003 verifies ALL24 equations M*L=h*e0. No GB or original script routine
is used in that independent multiplication.

Every cofactor has the common factor a^4*J^3, J=5P−6e−4. Exact FLINT division
removes it, yielding a primitive56-by1 cofactor with8918 terms and maximum degree11
for the much smaller target a*e0. A second direct multiplication against ALL24
source rows verifies M*Lprimitive=a*e0. The exact common factors and profiles,
original/primitive native data, both verification sources/receipts and failed
internal log are all preserved. A corrected replay script changes only the
internal target comparison to an explicit24-by1 matrix; it is not a new GB run.

These are complete finite-field polynomial identities for the RAW matrix only.
They do not prove QQ membership, original characteristic-zero rank41, or a
full inverse. The observed primitive support would require8918 QQ coefficient
unknowns; no QQ support solve has been launched. Cofactor degree11 is an observed
upper bound for this field lift, not a minimum. The separate degree-at-most4
obstruction concerns h targets, not this primitive target a*e0.

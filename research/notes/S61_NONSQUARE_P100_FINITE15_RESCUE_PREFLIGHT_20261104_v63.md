# Source-only one-border preflight over the finite H15 frontier

The exact chosen9-minor gcd H has degree15. Its monic coefficients are
32003-integral, so `S=Z_(32003)[z]/H` is finite free of rank15; neither its
QQ fiber nor special fiber is assumed to be a field or reduced.

The cheapest structurally admissible additional border from the saved
source matrix uses the accepted8-cover case9 (K2*v1 replaces q5), appends
K2*v3, and changes the ninth old local98 target to28. Its source order is
`[9,2,3,5,7,6,1,4,11]`, and targets are `[0,1,2,3,4,5,6,9,28]`.
The coefficient/support inspection uses no determinant value or rank.

The input has1043 RF terms, with1500 terms as the cleared-entry bound.
The one determinant has degree at most214 and215 output terms, with156735
bytes as the conservative raw JSON bound. Source QQ coefficients and H
coefficients pass the old-prime denominator gate. No source query occurs.

The safe future operation is to compute this one determinant N globally
over QQ[z] after only the five known-unit row LCDs, then reduce it by the
monic H. No old8-minor or finite-ring matrix pivot is inverted. A positive
gcd/Bezout `(Hbar,Rbar)=1` on the **entire** special fiber would imply
`S/(R)=0` by finite Nakayama, hence remove the finite QQ exception. This
requires p-integral coefficients but no leading-coefficient condition on N
or its residue: H alone anchors finiteness. No special-fiber open may be
discarded, including zeros of curve localizers.

Code/report stem:
`S61_NONSQUARE_P100_FINITE15_RESCUE_PREFLIGHT_20261104_v63`.
The source-only pass took0.178 seconds; no determinant or residue was
computed in this preflight. All bytes are immutable.

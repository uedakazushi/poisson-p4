# The saved P100 pure row is globally unimodular

Date: 2026-11-04. Exact characteristic-zero primary receipt; independent
review pending. This uses only the already saved actual pure98 source row.
It does not reevaluate any source or oracle coefficient.

Write the entries of
`S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63/pure98.json` as
`r_i=N_i/D_i`, in their saved order. All denominators are units in

```
T100 = QQ[z,1/(z(z+20)(z+60)(203z+2700)(z+12))].
```

The exact QQ[z] gcd of **all98 numerators is 1**. All98 quotient identities
`N_i=G*(N_i/G)` are checked. The extended gcd needs only positions0 and1:
the saved polynomials A0 and A1, of degrees21 and22, satisfy literally

```
A0*N0 + A1*N1 = 1.
```

Thus the two nonzero coefficients `A0*D0` and `A1*D1` form an explicit
regular dual for the whole98-entry row. The native verifier checks the
literal polynomial Bezout identity, and separately strips every original
denominator using only the five displayed linear curve units. No unknown
factor is inverted; in fact the numerator gcd is already1 before stripping.

The saved row is an actual source combination: producer q5H row0 minus its
full112-cubic H1 lift, with all pure cubic tails retained, followed by the
accepted NF117 operator. Its H1 and H2 parts vanish. The already certified
117 quadratic directions,23 H1 directions and one H2 direction give the
141 lower bound. At every geometric point of T100 the new pure row is
nonzero by its regular dual, so the original full weighted quartic map has
**pointwise rank at least142 on the whole P100 curve**.

This supplies a row-unimodularity certificate, not a selected fixed142
minor or an executed column completion. No rank98 or full239 claim follows.
The prior one-row note's unexecuted-gcd statement describes that earlier
stage; this separate receipt records the newly authorized completed test.

Producers are `code/s61_nonsquare_P100_pure98_row_gcd_20261104_v63.py` and
`.cpp`. Evidence is in
`reports/S61_NONSQUARE_P100_PURE98_ROW_GCD_20261104_v63`:
`native.json` contains all98 quotients, the two literal Bezout coefficients,
and their regular row-dual coefficients; `certificate.json` binds the actual
source row and original122 target labels. The job took0.992 seconds including
compilation, with0.028 seconds of native arithmetic and3235 stored terms.
It used no new source value, point, prime, minor, GB or Gaussian. There is
no automatic follow-up arithmetic.

Report SHA-256:
`f828554c9fc9123dd53038e8b69a9a0cffdb9ea98bf588319ee7e32e6260eda9`.
Native output SHA-256:
`e3d4540145ed4a618936a362c7c06dc93ed00700244c9a8167f85334808c4423`.
Certificate SHA-256:
`4de196041a86748985ad9d6f4ef151252066b87b2dfb66ce02042375fbf007be`.
All primary bytes are immutable.

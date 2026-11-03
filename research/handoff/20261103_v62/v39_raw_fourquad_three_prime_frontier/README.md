# Raw four-quadric homogeneous three-prime frontier

Extract after v36 and v38. This packet adds the two complete large-prime bases
at 2147483647 and 2147483629. They ran sequentially for 11.65 and 11.05 seconds,
maximum RSS about42 MiB, with sixty-second /2200000-KiB caps and no workers.
Both use exactly the SAME whole-column QQ homogenization degree seven, target
shifts zero, source pin and `(dp,C)` order. Both monic homogeneous reduced bases
have shape24-by-925,291812terms and maximum degree24, with the same leading
module as the32003 pilot. All original homogeneous generators reduce to zero.
The ten tail coefficients vanishing at32003 are explicitly zero-padded.

At a77-bit three-prime modulus,18438 of290887 nonleading coefficients (6.34%)
have equal balanced reconstructions at adjacent prefixes. The earlier late-source
pilot at the SAME77-bit modulus had7911/506068 (1.56%). Common-denominator
pooling observes20580 stable signed numerators (7.08%) and112 candidate-complete
columns, with pooled-denominator LCMs at most49bits. An independent pure-Python
Euclidean reconstruction repeats all18438 balanced agreements and constructs
all112 candidate vectors,6248terms total. Their degree histogram is
{7:44,8:3,13:57,14:7,15:1}. They are candidate relations; residue convergence by
itself proves no QQ membership. The separate exact-degree7 packet proves the
membership of44 of them by literal constant-column QQ identities. The other68
remain guesses for a subsequent augmented-module/Hilbert-function certificate.

The compact exact field annihilator shortcut is also saved: with zero-based
annihilator generator indices,540*A[2]+31787*A[3]=a^5*(5P−6e−4)^3 overGF(32003).
This is a coefficient-by-coefficient field polynomial identity, with no QQ lift.
Separate bounded coefficient tests exclude cofactors of degree at most3 even
overQQ for this fixed target; higher-degree QQ cofactors remain open.

The prime extension stops at three. No fourth field or automatic extension is
part of this packet. The derived CRT cache is only4.61MB and is omitted with
SHA/replay metadata; complete field bases are saved. There is no full QQ GB,
coefficient-height bound, QQ localized-surjectivity/rank41 or ambient conclusion.
All exact QQ GB/containment, p-integral monic, same-good-prime homogeneous-equality,
graded Hilbert-function and final QQ localized-membership gates remain explicit.

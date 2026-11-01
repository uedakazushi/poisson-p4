# The remaining small nilpotent strata at infinity

Date: 2026-11-01. Exact rational computation, with the v05 fixed-trace classifications explicitly retained as dependencies.

For sufficiently general beta, no point of the unsaturated 134-leading locus has either a regular-nilpotent modular matrix, or belongs to the open six-dimensional fixed-trace stratum of the minimal-nilpotent classification. Consequently these strata are also absent from the 135-leading locus. This statement excludes neither of the two nine-dimensional fixed-trace planes M0,M1.

## Minimal open stratum

The archived minimal-nilpotent certificate gives the full reduced fixed-trace fiber for s=x²/2. It is the union of two affine-nine planes and the closure of the open stratum (a,b)!=(0,0), which is an affine-three bundle over (A² minus zero) times A¹. Its dimension is six. The open stratum therefore has no additional exceptional subcases to omit.

The SL2 acting on (z,w) preserves the fixed modular matrix and acts on (a,b) as a nonzero two-dimensional covector. It puts every point of the open stratum in a=1,b=0. The six archived equations then give the following brackets, with four remaining parameters c,p,q,r:

```
Qxy=(1/2-c)x²-xz,       Qxz=0,          Qxw=-x²/3,
Qyz=(-3c²+3c/4)x²+(-6c+1/2)xz-3z²,
Qyw=-xy/3+p x²+q xz+(4c-1/2)xw+r z²+5zw,
Qzw=c x²+xz.
```

The symplectic substitution
(x,y,z,w) -> (x,y-Bw,z+Bx,w)
has pullback parameter action

c -> c+4B/3, p -> p+Bq+B²r, q -> q+2Br, r -> r.

All six identities are checked over QQ[c,p,q,r,B]. Choosing B=-3c/4 puts every point into c=0. On that entire three-parameter family, a fixed 12 by 12 minor of the beta-linear leading matrix is the nonzero constant -130977/327680. Thus its beta-constraint rank is at least twelve everywhere, including every value of p,q,r. Sp4 equivariance transports this bound to the whole open stratum.

The modular orbit has dimension four, independently checked by an exact commutator minor. Hence the affine stratum has dimension 4+6=10. The open stratum is conical, so its projectivization has dimension nine. Its incidence with beta has dimension at most

9+(16-12)=13.

Its image in W16 therefore has proper closed closure. The boundary (a,b)=0 is already in M0 or M1; it is not discarded or claimed absent.

## Regular nilpotent modular matrix

The archived regular-nilpotent theorem proves that the fixed-trace scheme consists of two affine lines, each a single centralizer orbit. Thus, up to Sp4 and an irrelevant projective scalar, all such points have one of the two explicitly stored representatives Q_A,Q_B. Each has an exact beta-constraint rank twelve, with nonzero rational minors stored in the report. Even the coarse parametrization by Sp4 and the four-dimensional beta kernel gives image dimension at most 10+4=14. Thus both complete strata are absent over a sufficiently general beta. There is no unexamined rank-drop parameter in this calculation.

## Scope and reproduction

Run:

```
.venv/bin/python code/s61_other_nilpotent_infinity_exclusion_20261101.py --output work/other_nilpotent_check.json
```

The frozen report is `reports/S61_OTHER_NILPOTENT_INFINITY_EXCLUSION_20261101.json`. It stores the uniform minimal-stratum minor, complete translation identity, modular-orbit certificate, regular representatives and their minors, with input hashes.

Historical dependencies are `legacy/v05/data/minimal_nilpotent_certificate.json`, `legacy/v05/proofs/minimal_nilpotent_theorem.txt`, and their regular-nilpotent counterparts. These are reduced fixed-trace exhaustion statements, not a new claim that the original nonreduced schemes have been classified here.

Together with the separately proved Jordan-(2,2) generic exclusion, the only nonzero nilpotent modular strata still requiring a complete generic infinity analysis are M0 and M1. They do survive: the separate minimal-nilpotent family certificate constructs two elliptic curves of such infinity points. Semisimple and mixed modular strata remain a separate task. None of these statements gives an upper bound for the full affine S61 fiber.

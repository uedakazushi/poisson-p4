# Generic-beta removal: exact source budget and remaining curve gate

The current rank strata do not by themselves remove the whole nonsquare
stratum for generic beta. They do provide a route strictly weaker than
uniform full239 on FULLU. The sufficient missing condition concerns the
actual relative bad incidence on each remaining core curve, not merely
the rank of the raw H projection.

## The actual all-nuisance source

Use Q0=xw partial_y wedge partial_z. In the K=5c+f=0,cL!=0 presentation,
with L=5b-e, the verified symplectic diagonal normalization sets c=L=1.
The normalized source H11 has core(a,e,d) and all eight coordinates

    u5,u25,u26,u27,u34,u37,u48,u51.

It maps into the actual primitive W16 through

    u4=a, u14=(e+1)/5, u35=1, u36=d, u49=e, u50=-5,
    u12=u13=0,

with those eight nuisance values unrestricted. The sixteen primitive
indices are4,5,12,13,14,25,26,27,34,35,36,37,48,49,50,51. The spatial
point is already fixed by the S61 setup; no extra four evaluation
parameters are added. Sp4 has dimension10. There is no additional scalar
factor in this parameter count.

The accepted full135 source action has positive nuisance weights
(2,4,2,4,2,2,2,2) and fixes the six core coordinates. It preserves each
full equation and its degree-two homogenization, including G6. The
normalized homogenizer t has weight-6. Only the source action must extend
to nuisance zero; extension of the entire projective scheme to lambda=0
is not assumed.

The verified Ga^3 subgroup fixes Q0 and every core coordinate on K=0 and
preserves the normalized source. Its action on Sp4 times the source is
free on the group factor and lies in the fibers of the joint map to
(beta,[Q]). This supplies exactly three legitimate fiber dimensions,
even when beta itself has a stabilizer. The two diagonal directions have
already normalized c,L; they cannot also be subtracted as extra fibers
of this H11 presentation.

Consequently for any core set E of dimension c,

    dim closure joint image(Sp4 times(E times A8)) <=15+c.

Finite core sets therefore have beta image of dimension at most15 and
are absent for sufficiently general beta, without any quartic calculation
at those points. A core curve gives only16, and can dominate W16.

## A weaker sufficient curve criterion

Let Z be the Q0-section pullback of the closed t-saturation bad support
of the full finite-type relative135 scheme over H11. It is closed and
Ga^3-stable; the verified source contraction gives
Z subset E times A8 when E contains the bad zero-nuisance cores.
Over the accepted general-beta smooth gauge open, faithful-flat base
change and equation equivariance identify this relative test with the
ambient bad support. A nilpotence calculation after imposing a single
fixed beta or just the transverse cuts is insufficient for this step.

For a geometric irreducible core curve C, C times A8 has dimension9.
If even one point is outside the FULL relative bad section, its closed
bad subset Z_C is proper, hence has dimension at most8. The same free
Ga^3 fiber argument gives

    dim closure joint bad image <=10+8-3=15.

Thus one full-relative good point, or equivalently a nonempty good open,
on each remaining C times A8 suffices to remove that entire curve's bad
contribution for generic beta. The actual bad set may still project onto
the whole core curve; it need not be finite in the core. No uniform
full239 theorem is required by this criterion. Alternatively one may
prove directly that the all-eight-nuisance beta map on a curve has rank
at most15 everywhere. That requires a genuine symbolic differential
identity, not subtraction of guessed stabilizer dimensions.

A concrete full-source t^N=0 identity, retaining all H11 parameters, can
establish the good-point condition without determining the whole normal
Hilbert function. The established free65 method is another sufficient
certificate: an admissible three-dimensional binary cut with actual
length at most65, the same dominating binary component's generic length
at least65, a retraction to its regular base, and the full-H11 relative
extension imply free65 and t^2=0. The regular base here has dimension14,
not the old H8 base11. The relative freeness, generic component and
ambient descent hypotheses cannot be replaced by low source ranks alone.

## What the accepted strata already give

The actual B7-by-3 bad locus on FULLU is finite or empty. Its contribution
is therefore beta-proper by the15-dimensional bound, regardless of the
individual exceptional quartic ranks.

The new raw-H classification leaves exactly two explicit Theta curves
P80,P100 in addition to such finite core exceptions. Both now have full
cubic57 at every permitted point. This alone does not establish their
transverse length or t-nilpotence. However the accepted older heavy-curve
quartic theorem applies: its15-by-7 block is a constant row combination
of the same actual16-by-7 heavy coefficients, with only a,F,B,C unit
scalings. A raw-H defect lies in that older defect locus. Its source-bound
proper-model theorem gives full quartic and actual length at most65
generically on each such curve. The accepted full-H11 free65/contraction
argument then leaves only finite bad core points there, whose beta image
is at most15. This reuses that accepted arithmetic and its combined
proof-audit note; no old field or special-fiber job is rerun here.

The common-frame gap B=D0=0 is a different core curve, outside B!=0.
Its all-eight-nuisance16-by-19 differential has a proved polynomial left
kernel and beta-image dimension15. This excludes that gap without
nilpotence, but the identity does not transfer to arbitrary core curves.

## The actual remaining gate

The accepted fixed-cut theorem bounds both the full cubic and full
quartic core defects by dimension1 on aFBC!=0. It does not put all those
curves in the raw-H defect locus. Unknown cubic or quartic defect curves
can still lie in the raw-H rank7 open. They have an available joint
dimension16 and are not removed by counting alone.

It would suffice to show that the actual residual bad core set there is
finite. More generally, for every remaining curve component it suffices
to prove either beta-image dimension at most15 with all nuisance values,
or that the full-relative bad section on C times A8 is proper as above.
The present finite-B and two-Theta-curve data do not yet supply this.

The checker saves a monomial logical countermodel: six quadrics have
product rank41; sixteen cubic rows give full rank57 and heavy rank6;
the low weighted Hilbert function is(1,7,23,34). Even adding H^2=0,
a pure weight-four monomial survives. This is not a Jacobi example. It
only rules out an inference from the accepted low ranks to full239.

Finally, this is the normalized nonsquare boundary chart. Its other
normalization boundaries and rank strata need their separate source
coverage results. Removing this boundary contribution would not alone
classify the affine residual S61 algebra, prove known-five exhaustion,
settle square strata, or describe special-beta local rings.

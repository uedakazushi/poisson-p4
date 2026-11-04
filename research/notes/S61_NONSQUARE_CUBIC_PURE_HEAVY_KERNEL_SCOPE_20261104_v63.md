# Pure projection versus the kernel of the heavy map

The proposed implication

`B7x3 drops rank => final14x21 drops rank`

is not established by the existing pure-rank certificates. The premise
that all eleven global pivots are pure is false: the old genuine pivot
is source label 11, target label 26, and that target is `H*L`.
There are ten quadratic pure pivots and one genuine heavy pivot.
This proposed implication was only sent for audit; no accepted ledger
theorem needs to be withdrawn.

The stage audit checks the frozen file headers and source metadata.
The old `completed-afterquad16x21.txt` has header `16 21` after ten
quadratic pivots. Five prescribed genuine heavy pivots, relative rows
`[1,2,8,10,11]`, then produce `completed-R11x16.txt`, header `11 16`.
The accepted generic canonical rank ten applies to those eleven
remaining source rows **after the five heavy Schur steps**. It is
not a rank upper bound for the fixed pure projection of all sixteen
original genuine rows. The five removed heavy-pivot source rows
can have additional pure components.

In the current route the single old genuine `H*L` pivot is already
among the eight completed pivots. Its normalized row has nineteen
nonzero pure coefficients. Removing this mixed source direction
cannot be counted as removing a pure target direction. The initial
augmented target has pure/heavy counts `26+6`; the eleven actual
pivots leave `16+5`, exactly as the accepted `14x21` presentation
records. The later GL6, GL3 and GL10 changes are pure target changes
with zero heavy tails in their quadratic pivot sources. They do not
change the classification of the earlier genuine heavy pivot.

In particular the old canonical rank ten establishes a lower bound
of ten for an appropriate pure image inside the remaining source
span. It does not prove that image has no additional pure directions
from the eliminated heavy rows. The claimed initial augmented pure
rank twenty was not proved as an equality by that certificate.
Consequently there is no proof here that all ten-minors of the
fixed final `14x16` pure projection vanish identically.

The distinction can already be seen in the matrix

`[[1,0,x],[0,1,0]]`, with the first two columns pure.

On D(x) a heavy pivot leaves one pure row. Nevertheless the whole
pure projection has rank two, and the full matrix still has rank
two when x=0 and its heavy rank becomes zero. A generic Schur
remainder therefore cannot supply the asserted global pure bound.

## The correct pointwise criterion

At any residue field write the actual final matrix as the map

`M=(P,H): k^14 -> k^16 direct-sum k^5`.

Put `r=rank(H)` and `K=ker(H)`, so `dim K=14-r`. Then

`rank(M)=r+rank(P restricted to K)`.

Indeed the kernel of M is precisely the kernel of the restriction
of P to K. Full source injectivity is equivalent to that restricted
map having rank `14-r`. On a heavy rank-five chart this is a
nine-source map into the sixteen pure targets. Where heavy rank
is at most four, its kernel has dimension at least ten; those ten
or more directions can still inject into sixteen pure targets.
There is no automatic total-rank upper bound thirteen.

There is a concrete finite-DAG chart recipe. Choose r source rows
I and r heavy columns C with nonzero minor Delta. For each source
u outside I, set

`k_u = Delta*e_u - sum_i lambda_ui*e_i`,

where `lambda_u=H[u,C]*adj(H[I,C])`. On the stratum `rank(H)=r`,
these vectors kill all five heavy columns and form a kernel basis
after inverting Delta. Their actual pure images are

`Delta*P_u - sum_i lambda_ui*P_i`.

All pure entries may remain references to the accepted coefficient
circuits. An unknown Delta is only a denominator of this chart,
not a globally certified unit. On lower-rank strata a corresponding
lower-r chart and the extra heavy rank equations are required.
This construction keeps precisely the pure directions relevant
to full injectivity rather than the unrelated fixed pure projection.

The accepted B-block rankbad-locus work controls part of H. If that
work produces actual points where B has rank at most two, then H
has rank at most four there, but the kernel-pure map must still be
tested. No full-cubic failure locus or its existence is concluded
from this audit. The accepted global `14x21` presentation, heavy
equivalence and generic full-rank certificates are unchanged.

No coefficient of degree207, new point, determinant, rank, GB or
full pure matrix was computed in this stage audit.

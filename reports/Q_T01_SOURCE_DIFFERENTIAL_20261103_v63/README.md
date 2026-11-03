# Qa / actual T01 source differential

Full canonical homogeneous150 coefficients over QQ(s), no projective
quotient and no truncation of the universal matrix-power map.

`report.json` records exact RREF ranks and nonzero minor witnesses.
`matrices.json.gz` retains both150x31 Jacobians, the independently computed
strict GL5 action, the stabilizer and actual integrable fiber directions.
All indices in the report and sparse matrices are zero-based.

The key ranks are24/24 for the two source images,30 for their sum,18 for
their intersection, and7 for their complete source-fiber tangent kernel.
The intersection is the17-dimensional strict GL5 orbit plus the Qa curve.

The mathematical conclusion is singularity of actual T01 at ALL finite
physical Qa, including a=0 and the ramified cover values.  The affine T01
tangent dimension is at least30 everywhere on that curve; its dimension
is26.  Source-fiber regularity7 is asserted only at the generic displayed
source point.  This packet does not count local branches or normalization
points, prove a completed ideal, or exhaust all Poisson components.

Proof: `research/notes/Q_T01_SOURCE_DIFFERENTIAL_20261103_v63.md`.
Reproduce to a fresh ignored directory:

```sh
python code/q_t01_source_differential_20261103_v63.py \
  --output work/q-t01-source-differential-v63/replay
```

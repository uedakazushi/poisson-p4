# μ5/r1 marked quartics: four full-component candidates and one frontier

2026-11-07. PRIMARY E. Producer proof, independent review pending. No global count is changed.

The complete marked-family argument is in `PLANE_HYPERPLANE_QUARTIC_MU5_MARKED_LIST.md` and `PLANE_HYPERPLANE_QUARTIC_MU5_A5_A4A1_MARKED.md`. The complete quadratic Saito/source batch and all-six-magnetic orbit normalization are in `PLANE_HYPERPLANE_QUARTIC_MU5_SOURCE_NORMALFORMS.md`. Each of the five actual sources has one dense projective orbit of dimension 23; its entire six-magnetic image, including the zero-invariant boundary, is in that orbit closure. The known103 A3+A2 boundary control is not one of the five ambient tests.

## One complete five-source ambient batch

The batch uses all 150 homogeneous coefficients and all 350 projective Jacobi equations. It constructs Euler24, all25 linear-coordinate action directions, scale, and all six actual magnetic derivatives. In every case the symmetry and whole-source ranks are 48, Euler rank is 24, and the full derivative annihilates every actual direction. The primitive integer rescaling changes no projective point. There is exactly one structurally justified canonical point per whole orbit, with no point search.

| Whole source | Canonical R | rank D | Conclusion |
|---|---|---:|---|
| A5 smooth hyperflex | h² |100| two normals; stop, membership unknown |
| A4+A1 node inflexion | hy |101| nonzero one-normal quadratic obstruction; component23 candidate |
| A4+A1 smooth hyperflex | h² |101| same component23 proof |
| D5 smooth hyperflex | h² |101| same component23 proof |
| 2A2+A1 smooth hyperflex | h² |101| same component23 proof |

For each positive case the report gives a genuine 48-dimensional actual symmetry subspace, one normal, 101 implicit directions, and a nonzero 150×150 direct-sum determinant. The normal is in the full derivative kernel and lies in the actual linear complement to the 48 symmetry directions. The full 350-vector quadratic term is evaluated against an explicit functional annihilating all 150 derivative columns. Its nonzero value is, respectively,

- `6041932323756943536/1516155105099253729`;
- `-1456129679381483677863776/1070464342813498931805`;
- `632765198693858076/11192591940071`;
- `-21729812880435/58571233780238336`.

The implicit function theorem eliminates the 101 derivative directions. The resulting one-variable residual has nonzero quadratic leading term, so its reduced support is the origin. Thus the raw local dimension is at most 48, equal to the actual orbit dimension. The irreducible 23-dimensional projective orbit closure is a whole ambient component. This argument claims no smoothness: rank D is 101, so these component germs can carry scheme thickness. It uses no higher jet or Gröbner basis.

For A5 the tangent upper bound is 25 and the actual orbit dimension is 23. Neither a component statement nor known-parent membership follows. Its two-normal frontier is retained without a higher-order or new-anchor computation.

The exact report is `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU5_AMBIENT_20261107.json`, phase COMPLETE, 44.267 seconds and maximum RSS 81,600 KiB. The launch with an incompatible external CPU hard limit failed before any mathematical computation; the single actual batch used the declared 180-second limit and completed. D5 cusp-tangent and the known103 A3+A2 control were not recomputed.

## Actual prior111 comparison

Use the immutable public111 catalogue at commit `37d56e2afb3fd467401d395e39f9c34f32bae7b3`, blob `dfc7c037b2050012fbe65c6c34b0e1ad53d3ff59`, projected in `inputs/mu6_prior111_37d56e2a_20261107.json`. The old five same-dimensional H+quartic components require actual tensor/divisor witnesses, rather than names. Those actual witnesses and their source bindings are detailed in `PLANE_HYPERPLANE_QUARTIC_MU6_COMPONENT_FRONTIER.md`; the argument here applies them to the new four sources as follows.

1. ESSENTIAL5 has an actual quintic directional-derivative minor of determinant −5, hence essential rank five. Every new orbit closure satisfies the closed essential-at-most-three condition.
2. MINIMAL_M has rank-one canonical cofactor on its dense23 orbit. The complete new sources have rank two.
3. A0 has actual quartic `h4+h²xz−z4/4`, of D5 type, with H=h passing through the D5 singularity as its cusp tangent. The new D5 smooth-hyperflex H avoids that singularity; the other three positive cases have different quartic singularity profiles. These are distinct marked pairs.
4. DPM has actual tricuspidal quartic (3A2), not any of the new μ5 profiles A4+A1, D5 or 2A2+A1.
5. Ga9constant has distinguished H inside the dual image of its actual canonical cofactor, whereas every new source has H outside. Equal23 components would force the two dense constructible actual orbit/source images to intersect, contradicting this intrinsic marked-tensor distinction. No unproved closedness of image incidence is used.

All other prior111 entries are separated by generic rank, full component dimension, or an explicitly closed divisor multiplication/nonreduced locus. The same23 MPS29 each contain an accepted actual toric-SNC point with nonzero quintic κhxyzw, which has essential rank five. That actual point cannot lie in the new closed essential-at-most-three component. This is a component-equality exclusion, not an inference from one boundary membership. Old irreducible-quintic components have an actual irreducible point outside the closed linear-factor locus containing every new orbit closure.

The four candidates are pairwise distinct: quartic singularity profiles distinguish all except the two A4+A1 cases. In that pair H contains the node in the node-inflexion case and avoids both singularities in the smooth-hyperflex case. Equal23 orbit closures would force equality of their dense23 actual orbits, impossible for these intrinsic marked incidences.

They are also distinct from all nine μ6 component candidates: the underlying irreducible quartic has total Milnor number five here and six there. This comparison does not presuppose adoption or counting of the μ6 candidates.

## Review and remaining scope

The four componenthood proofs, whole-family geometry, Saito/full-magnetic normalization, all-prior comparison and the A3+A2 known103 boundary bridge still require independent acceptance before any count changes. A5 smooth hyperflex remains unknown. The separately stopped E6 smooth-flex graph diagnostic does not exclude all parents. Ga9 E6 cusp-tangent/linear-w and stable/nilshort residual branches are unaffected by these four proofs.

The current runner is frozen in an exact public byte snapshot. The dependency manifest contains all geometry/source reports and certificates needed to replay the batch without a work-directory file. Historical jobs/reports are immutable even when a later isolated runner mode preserves their function ASTs.

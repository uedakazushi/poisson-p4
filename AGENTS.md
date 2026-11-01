# Poisson structures on P4 — Codex instructions

Version: v37-s61 · 2026-11-01 (Asia/Tokyo).

Current v37 additions: the entire canonical nonbinary xzw orbit is removed
by t-saturation over sufficiently general beta. A universal43-minor and one
square obstruction cover its open; five16D artifact divisors plus the earlier
nodes/zero49/proper-image exclusions cover all factors. Other R0 types z^2*w,
z*w^2 and mixed tails are separate. Binary rank-drop points admit free rank65,
I^4=0, nonzero t with I*t=t^2=0 and quotient64; transverse Hilbert function is
(1,8,22,34), with NO graded-piece freeness claim. Three16D endpoint strata u4=0, u14=0,
u49=0 are removed over general beta; not the whole17D rank-drop curve.
All four isolated M0/M1 g6/g7 axis presentations are removed over general beta,
using original134 full Jacobian rank45. Eight minimal junctions have generic ranks4,4,4,8,3,4,5,6 and nilpotent t;
finite-jet comparison removes their16D strata over general beta. Positive-
dimensional tails and other endpoints are not thereby covered.
AB normalization has D3=B3 dimension6691 and H1(E3)=0; all-degree generation
and origin normality remain open. Full affine finiteness/five-point exhaustion
and global component upper bound remain unproved; lower bound57 unchanged.

Previous v36 results: on the specified binary open, the completed full135
ring is finite free rank65 over C[[v0,v1,v3]], with reduced ring that base,
nilradical I^4=0, nonzero t in I^3, I*t=0 and quotient(t) free rank64.
The rank65 result does NOT cover all generic binary points: a dominant17D
binary rank-drop family has leading/full ranks33/34, hence gives a curve over
generic beta. A separate Delta!=0 example even has ranks32/33.
The nonbinary u49=0 endpoint is nonflat: rank7 over C[[u]], transverse length12,
three generic branch lengths2,3,2, finite-length torsion with5 minimal generators.
t is nilpotent; full t^2=0 is not asserted. Finite-jet specialization extends
saturation removal to the entire generic-beta fiber of its16D stratum.
The four corrected minimal baseline lines have square-zero local opens;
an M1 resonant tail has C[[u,e]]/(e^3), t=0. Junctions/endpoints are separate.
AB normalization has no extra degree2 generators: D2=B2 dimension761.
Origin normality and all-degree generation remain open.

Previous v35 results: generic reduced UNSATURATED infinity of the full135
homogenization is exhausted by R0 and minimal M0/M1 normal forms, conditional
on the historical Q45 all-point classification. Non-nilpotent, regular nilpotent,
nil22 and the other minimal open stratum have beta images of dimension<=15.
The corrected minimal list includes FOUR baseline lines missed by cancellation
of generic solution denominators. All case splits have independent review.
This does NOT remove all infinity by saturation or prove affine finiteness/five points.
Two nonbinary resonance nodes have full135 completed local rings free rank4
over C[[u]], complete intersections with reduced ordinary node, square-zero
nilradical and t=0. Their16D strata and the M0 isolated-axis16D component have
their entire generic-beta fibers removed. Other curve endpoints remain open.
Binary infinity has universal full Jacobian rank35, quadratic obstruction rank5,
and cubic relation rank<=50 via six full-Bianchi relations, giving normal length>=65.
Use the SAME smooth3D binary locus and constant-rank hypotheses when combining
this with the specified transverse quotient of length65. Do not extrapolate to all points.
The actual reduced T01 zero-first-jet fiber is closed AB, irreducible dimension10.
AB is normal AWAY FROM THE ORIGIN. Its two nonzero one-variable orbits have
embedding dimension/multiplicity16/7 and18/10. Origin normality and scheme fiber
remain open. Reduced connected fibers alone do not prove normality; a counterexample
is recorded. Prior minimal-curve square-zero formal rings and full135 equivariance
remain valid on their stated opens. Known reduced T01 U T39 U T40 is finite etale
degree5 over ALL Delta!=0, but full S61 is not exhausted. Global lower bound57 unchanged.

Current continuation: T39=q(5,1) and T40=q(5,2) reduced S61 closures are
explicit homogeneous quadratic graphs over ALL W16, hence both have degree ONE.
The complete generic parameter redundancy is 3000 for each. The two graphs
intersect in the five-dimensional degree-24 rank-one flag cone, with a length-10
square-zero nilradical at the vertex. T01 has degree THREE and geometric
monodromy S3. A canonical rational affine THREE-PLANE contains all five known
generic lifts; its intersection with the full 134 equations is reduced length
five with algebra A_T x K x K. CONTAINMENT OF THE FULL FIBER IN THIS PLANE IS
UNPROVED. The next exhaustion problem is to derive its 42 affine linear
conditions in the radical of the generic full ideal. The toric three-cover is now proved isomorphic to Jac(C_beta)[2] minus zero.
The canonical discriminant is geometrically irreducible, reduced, degree24.
C9 is globally polynomial: a projective pencil avoiding rank drop proves
codimension>=2, then normality extends the rational solution. A polynomial
weight-four X in beta,q gives the inverse coordinate. Therefore the ACTUAL
reduced T01 over Delta!=0 equals the finite etale monic cubic cover, degree3.
Globally T01 is the explicit affine modification described in the v31 notes,
not the finite normalization of the base. Its beta=0 fiber now has dimension exactly10;
the finite algebra model misses every nonzero point of that fiber. Do not
continue treating the algebra-model boundary as potentially exhaustive.
These statements do not exhaust full S61. The plane-adapted135-equation F4
probe hit600seconds without a completed basis; do not repeat it unchanged.
Do NOT add curve-invariant transport equations to full S61: an exact T02
counterexample has smooth canonical curves with varying j, even on the usual
frame/invariant open. T02 first-jet image is an irreducible hypersurface.
The 135-equation augmented input and infinity obstruction remain valid; a
local infinity point has dimension three, but this proves no affine finiteness.
Global component lower bound remains 57. Read current STATE/NEXT_TASKS and
code/S61_README.md; do not reopen solved elliptic degrees as unknown.

## Read first

Before mathematical work, read **RESEARCH_POLICY.md**, then
`research/STATE.md`, `research/NEXT_TASKS.md`, and the relevant entries of
`research/CLAIMS.json`. Read `docs/CODEX_SETUP.md` for execution details.
The current policy and state supersede old HANDOFF files as a description of the
frontier; historical mathematical claims remain source material to assess.

## Objective and scientific standard

Continue research toward complete classification of complex Poisson structures
on projective four-space. The first completion target is exhaustive reduced
irreducible components of the projectivized Jacobi scheme, with their generic
geometry, **before** quotienting by PGL5. Retain orbit classification,
intersections, and nonreduced structure as subsequent targets; do not silently
declare the full original problem solved after only the first target.

- The current lower bound is 57 = 40 MPS + 6 generic-rank-2 + 11 further
  generic-rank-4 components. It is not a proved total. T37, SS/NN, CA and CB must
  not be counted again.
- Rank <= 2 exhaustion has proofs in v14 and reproduced certificates. Use it
  unless a concrete flaw requires reopening it. Full independent review is not
  complete. The normal-actual-divisor theorem in v22 has an outstanding
  dependency audit; `D normal` and `D_red normal` are different hypotheses.
- Prioritize the **full S61 model**, including all 16 weight-1 directions.
  The 45-variable section and the v24 endpoint Q are not global exhaustion
  results. Separate tangent spaces, parameter differentials, quadratic
  obstructions, formal local rings, and reduced components.
- Distinguish established literature, bundled proofs, independently verified
  identities, conditional statements, conjectures, and observations. A passing
  program or a finite-field sample does not by itself prove global exhaustion.
  Record hypotheses, dependencies, proof locations, and unresolved steps.
- Maximize the value of the whole research program. Change method or priority
  when evidence warrants it; record why and what remains of the original goal.

## Work autonomously and preserve continuity

After a useful intermediate result, select the next justified task and continue
within the actual runtime and resource budget. Do not stop merely to ask whether
to continue, to save a file, or to perform routine reversible work already
authorized. Use bounded subproblems with explicit mathematical completion
criteria. Delegate independent exploration and verification when available.

Stop substantive work only on completion, user interruption, actual resource or
access limits, or a serious obstacle that cannot be responsibly resolved with the
available information. On interruption or a limit, save a checkpoint containing
the exact state, commands/results, blocker, and next executable step. Never mark
a goal complete because a time or token budget expired. These instructions do
not create a background scheduler or override platform permissions.

## Commands

```sh
bash scripts/setup.sh
.venv/bin/python scripts/doctor.py
.venv/bin/python scripts/verify.py --suite quick
.venv/bin/python scripts/verify.py --suite rank2
.venv/bin/python scripts/verify.py --suite integrity
bash scripts/build_notes.sh
.venv/bin/python scripts/make_bundle.py
```

`sources/v02.zip` ... `sources/v24.zip` preserve every supplied file byte for byte.
`scripts/bootstrap_sources.py` expands them into ignored `legacy/v02` ...
`legacy/v24`, with relative previous-version links. Do not edit imported sources
or their stored certificates. Put new work in `research/`, `papers/`, or `code/`.
Old runners that generate outputs must run through the wrapper or in a temporary
copy. See `sources/README.md` for source provenance.

## Updating and publishing work

- Write manuscripts in Japanese LuaLaTeX unless instructed otherwise. Do not use
  writable boxes or forms. Keep useful failed attempts and hard-to-recover
  observations in research notes; omit redundant execution chatter.
- Keep policy, state, claim/component registers, notes, manuscript sources,
  references and reproducibility commands synchronized. Separate unproved draft
  sections visibly. Audit primary literature before making novelty claims in
  a publication manuscript. Split or combine papers when that improves value.
- Work on **main** and push directly to `origin/main` at substantive checkpoints.
  During long work, also checkpoint approximately every 30 minutes when there
  are meaningful changes, including honestly labeled work in progress. This is
  an agent responsibility, not a timer installed by this repository.
- Before committing, inspect `git status` and the diff, stage only intended
  research files, and run relevant checks. `scripts/checkpoint.py` commits only
  already staged files; use `--push` to push. Never force-push, erase history, or
  discard another session's changes. If the remote advanced, fetch and integrate
  after inspection; do not overwrite it. Keep one integration writer when using
  parallel agents. **Do not add GitHub Actions workflows.**
- Versioned continuation bundles are authorized without asking the user again.
  Include policy, notes, available publication drafts, sources, evidence,
  bibliography and reproduction instructions. If ChatGPT Library tools are
  available in the executing environment, the user authorizes saving those
  bundles there without confirmation. Record actual save receipts. If that
  capability is unavailable, keep the checkpoint in GitHub, record the pending
  Library save, and continue unaffected research. Never claim an upload occurred
  without a successful receipt.

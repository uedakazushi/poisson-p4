# Poisson structures on P4 — Codex instructions

Version: v33-s61 · 2026-11-01 (Asia/Tokyo).

Current v33 additions: the actual reduced T01 zero-first-jet fiber has
dimension EXACTLY10 with AB as its unique10D component; any additional
component has dimension<=9. This uses the historical Q45 classification and stabilizer
bounds. This does not prove irreducibility or exhaust lower-dimensional parts.
Over ALL Delta!=0 the known reduced union T01 U T39 U T40 is disjoint finite
etale of degree5. This is not an assertion about the full S61 fiber or its
ambient local rings at all five points.
For sufficiently general beta, all reduced R0 infinity normal forms have been
exhausted, including resonances and rank-one quadratics. Proper incidences and
uniform beta-rank bounds exclude C0/L0 outside R0. Conditional on historical
Q45 classification, this exhausts the zero-modular reduced infinity locus.
At one nonsquare binary-plane point, cubic obstructions prove local dimension3
for the135 homogenization. Its completed local ring and t-saturation remain
unknown. Earlier flag/nonbinary completed local rings and t=0 remain proved.
A full135-leading nil22 counterexample has Delta!=0 and nonzero modular matrix
not stabilizing beta: do not assume either condition forces zero modular part.
T02's invariant equation has degree>=26; do not repeat fitting through24.

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

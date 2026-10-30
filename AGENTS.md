# Poisson structures on P4 — Codex instructions

Version: v26-s61 · 2026-10-30 (Asia/Tokyo).

Current continuation: T01 has dominant, generically finite primitive first-jet
projection to the full 16-dimensional first-jet space. This is not global
exhaustion. Read the updated STATE/NEXT_TASKS and current code/S61_README.md.
Rank-one first jets reduce to a preserved 45-variable, 134-quadratic fiber.
The CA and CB slice components have the same closed nine-dimensional first-jet
image, so they do not dominate the 16-dimensional first-jet space. Among the
known components, only T01, T39, and T40 can dominate; the latter two and all
potential unknown components remain to be assessed.
Do not reopen generic first-jet realizability as unknown; all extensions and
exceptional fibers remain open.

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

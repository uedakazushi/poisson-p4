# Poisson P4 — classification first

Version: v64-classification-plan / 2026-11-06 (Asia/Tokyo).
Mathematical baseline: v63 CO, a00a541b835b2bff18bafcba545ecbd4bae7175b.

## Start here

Read `RESEARCH_POLICY.md` first, then `research/classification/CURRENT.md`,
`research/classification/PLAN.md`, and `research/NEXT_TASKS.md`.
Run `python3 scripts/classification.py check` and `python3 scripts/classification.py next`.
The active queue is `research/classification/PLAN.json`.
`research/STATE.md`, old handoffs and old reports are mathematical evidence, NOT an active queue.
Read only the evidence needed for the selected task. Do not load the whole historical corpus.

## Primary outcome

Write a paper whose main theorem EXHAUSTS the reduced irreducible components of
Poisson structures on complex P4, before PGL5 quotient, and describes their general members.
Prove coverage of all points by the resulting closed families, componenthood and distinctness.
The present lower bound is 99, including 93 of generic rank four; the total is UNKNOWN.
All99 are individually identified, including the new rank3 double-quadric components of dimensions22 and24.
The entire H+2Q locus with Q irreducible, three distinct hyperplanes plus an irreducible quadric, and irreducible quadric plus irreducible cubic are covered. Other degeneracy divisors remain open.
Do not fix the lower bound as the answer. Add genuinely new components when necessary.
Full orbit classification and full nonreduced scheme structure are not claimed by this theorem.
They are secondary, not prerequisites silently added to the primary target.

## Decision rule before computation

Every nontrivial run must name an active classification obligation and the precise implication
from its result to that obligation. Prefer the WEAKEST sufficient statement: reduced support,
generic containment, or non-dominance, not a full completed ring or exact nilpotence exponent.
Use `PROOF_ROUTE.md` for the generic radical, scalar, and closed-complement routes.
A new minor, a longer jet, a smaller matrix or a new wrapper is not progress by itself.
After two nonclosing attempts on the same route, diagnose the mathematics and change route;
do not continue an open-ended sequence of coefficient probes.

## Work allocation

The default mathematics worker handles generic reduced exhaustion; a second worker handles
exceptional-base coverage/component audit. Use an independent reviewer when available.
The integrator alone writes main. With one worker, alternate between these PRIMARY tasks
when blocked. Do not use capacity for suspended secondary tasks merely because they are easier.
B18 exact nil index/H0/T2, higher contact/Betti data, and separate-paper polishing are SUSPENDED.
Reactivate only after recording an explicit indispensable dependency on the classification
proof and why the existing weaker result is insufficient. No separate publication sprint.

## Proof discipline

Separate repository-reported proof, proof checked in this run, exact computation reproduced,
conditional result, conjecture and observation. The plan checker does NOT verify mathematics.
Keep all 134 original equations, inverse maps, compatibility conditions and excluded charts.
Never infer full support from a chosen slice, a few primes, rank of a submatrix, a formal jet,
normalization, a quotient, or an unproved flatness assumption. Generic in FULLU need not mean
all points above generic beta. A boundary rank theorem is not itself affine five-point exhaustion.
Record every localization and its complement; distinguish base parameters from unknown fiber
coordinates. A finite algebra over K[n_1,...,n_42] need not be finite over K.
Retain the exact input hashes and a separately checkable certificate of each main claim.
Do not alter immutable producer bytes or old certificates to make a check pass.
Do not re-open solved results without an identified defect or missing application hypothesis.

## Autonomy and resources

Continue useful mathematics after checkpoints without asking to continue. This instruction
applies only within a running Codex task; it creates no daemon, scheduler or unlimited quota.
Use finite jobs with a mathematical success/failure criterion, measured resource limits,
checkpointable state and a fallback. Failed runs cannot be called proofs.
Old point/prime/support limits describe old experiments, not a perpetual restriction on new
classification methods. New primes, coordinates, algorithms and sensible resource budgets
are allowed when justified; keep exact characteristic-zero validation separate.
Use one reusable runner, not a growing family of guard/custody wrappers. Preserve only useful
counterexamples, irreproducible inputs and certificates, not execution chatter.

## Environment and validation

Reuse the existing environment when it works. `docs/CODEX_SETUP.md` is an environment reference,
not a source of research priorities. Do not reinstall dependencies or rerun every old suite.
Relevant existing commands: `.venv/bin/python scripts/doctor.py`,
`.venv/bin/python scripts/verify.py --suite quick`, `--suite rank2`, `--suite integrity`.
Quick/integrity PASS does not certify the global geometric classification.
New control-plane tests: `python3 -m unittest discover -s tests -p test_classification.py`.
No sudo, global Codex configuration changes, unapproved service setup or GitHub Actions.

## Save and publish

Write mathematics and manuscript sources in Japanese LuaLaTeX unless instructed otherwise.
Update the CURRENT/PLAN files IN PLACE; do not prepend dated progress logs. Preserve valuable
results in current papers or mathematical notes. The new main manuscript is
`papers/poisson_p4_classification_ja.tex`; its classification theorem remains a TARGET until proved.
Old papers/certificates remain intact as evidence until their useful mathematics is integrated.

Inspect status/diff; stage only intended files; use ordinary commits and non-force pushes to
origin/main at substantive checkpoints. If the remote advanced, inspect and integrate without
erasing another session's work. Do not change the claim catalogue merely to match a plan.
The user authorizes current, self-contained research bundles and Library saves without further
confirmation. Generate/validate/upload are separate states; report only actual receipts.
Do not package old ZIPs, logs or versions as a substitute for a current dependency-complete bundle.
Do not delete existing Library records without an explicit instruction.

# Poisson P4 — work through the final classification deliverable

2026-11-09 / v72-normalization-residue-integration

## User mandate: one end-to-end assignment

The user now explicitly asks for a precise final goal and continuous execution
until that goal is reached. Read RESEARCH_POLICY.md first, then
research/conceptual/FINAL_TARGET.md (the terminal contract),
research/conceptual/GLOBAL_PLAN.md (the route), and
research/classification/CURRENT.md (the live evidence/status).
This replaces milestone-sized stopping rules, the old stop-after-147 rule,
and any instruction to stop after one conceptual lemma or the next manuscript section.
Historical done/PASS fields are evidence of their stated scope, not this task's completion.

CONTINUE working within the available execution until M1–M4, C1–C2, P1, R1,
and S1 in FINAL_TARGET are all satisfied. A plan, new theorem, completed case,
passing tests, PDF, commit, push, or bundle is an intermediate checkpoint.
After saving a checkpoint, proceed directly to the next unmet requirement.
Do not end with “next steps”, “shall I continue?”, or a request for another user
instruction for already-authorized work. Short progress updates are fine.

## Precise outcome

Complete a proof-bearing, publication-ready Japanese LuaLaTeX classification
of the reduced projectivized Poisson scheme on complex P4 BEFORE PGL5 quotient.
Construct explicit irreducible families with saturated closed images C_i,
prove all-point coverage and pairwise noncontainment, and give each component's
generic construction, dimension, rank, and divisor/foliation geometry.
The current expected result is 147 components (141 rank four and 6 rank two).
Derive the count; do not assume it. Correct actual mathematical defects and all
dependent statements instead of fitting the proofs to 147.

The exposition must derive the classification through recognition → marked
reconstruction → actual Poisson lifting/deformation → maximal distinct closed
images. Reusable theorems must replace repeated proofs in the current manuscript.
Keep necessary exceptional cases and finite computations explicit and reproducible.
No requirement to eliminate every computation or achieve an absolute minimum
number of theorems/pages. Do not stop at a framework that leaves its applications open.

## Continue from the existing mathematics; do not restart it

The four shared files already integrate the marked-hyperplane/modular primitive
argument, two lifting mechanisms, quadric cofactor and reduced free plane-quintic
magnetic residue quotient, and the plane conductor/conormal reconstruction:
- papers/marked_hyperplane_reconstruction.tex
- papers/poisson_lifting_examples.tex
- papers/quadric_cubic_reconstruction.tex
- papers/conductor_conormal_reconstruction.tex
Read their actual proofs and preserve their valid scope. Novelty is unconfirmed.
The component navigation table is not a proof. Reuse verified inputs without
rerunning irrelevant historical experiments or assuming old PASS is sufficient.

The fifth shared file, papers/normalization_residue_reconstruction.tex, now
reconstructs the entire reduced fixed-D fibre as a compatible linear Poisson
space minus a quadric. Normalization uses all bounded-conductor closed one-forms
and the full B/A linear descent condition, including higher principal parts.
Coefficient-one/node residues, full scroll intersection kernels, finite MCM
partial conductors and unique relative ambient lifts replace duplicate degree5/6
and both actual T40 family arguments. All10 plane coefficients have an explicit
rank2 determinant boundary. Distinct ambient components with generically reduced
divisors cannot have the same closed divisor image; image inclusion alone still
does not prove component inclusion at kernel-rank jumps.
The new37 checks include one full150 toric kernel, three full M/R sources and a
node-cusp matrix factorization. They do not certify all147 or final dependencies.
These results are intermediate: CONTINUE through the remaining terminal contract.

Proceed through nonplane conductor recognition and normalization/descent,
rank-two support, maximal closed images, final manuscript integration,
primary-source comparison, scoped mathematical/computational checks, the current
self-contained bundle, main update and Library save. GLOBAL_PLAN identifies
starting files and dependencies. Do not stop after the nonplane stage.
Change method or order when blocked; a failed approach is not the end of the task.

## Mathematical guardrails

- A/B/C/R2 are entrances for points/generic families, not disjoint component closures.
- Weak P#df=fZ needs f!=0. Strong df+i_Z(alpha)=0 with the other primitive
  conditions is valid also at f=0; do not confuse the two statements.
- Simple H has deg f=4; multiple H has deg f<=3. Pf(P_2)=0 is not a simple-H assumption.
- r=1 includes 5H as well as integral quintics. Canonical symmetry zero is not all symmetry zero.
- Preserve translation cokernels, denominator and zero/unit complements,
  special-base components, markings, and generic/geometric-generic fibre distinctions.
- A divisor deformation need not lift. Nonlinear auxiliary maps do not prove PGL equivalence.
- The magnetic nu−1 formula needs the reduced free plane-quintic hypotheses;
  no automatic extension to nonreduced factors or arbitrary conductor surfaces.
- Keep the hidden binary-quadratic three parameters and all higher terms in
  the conductor source. The second Jacobi equation is inhomogeneous.
- Do not assume conductor reducedness, normalization-flatness, or simple normalized residues.
- Rank-three H+2Q has excess dimensions 1 and 6. Its local Koszul project is not
  a prerequisite to finishing the global classification.
- MPS40 concerns components meeting SNC. Dimension, a normal form, a generic
  divisor type, or a control-plane check alone proves neither coverage nor componenthood.

## Iteration, integration and completion

Use the unmet FINAL_TARGET requirement and its shortest necessary dependency
to choose work. Diagnose a method that repeatedly produces no structural result,
change the route and CONTINUE; do not convert that diagnostic into a stop rule.
Update CURRENT in place with proof labels, verification scope and remaining gaps.
Do not accumulate worklogs or duplicate claim ledgers. Plans support research;
revising plans alone is not a research milestone.

Retain all valuable current results and sufficient proofs; replace duplicates
by shared proofs and valid references, not by deleting content or parking old
manuscripts in an archive. The final paper may be one coherent manuscript or
several logically self-contained papers. Read primary references for material
claims and distinguish established novelty from unconfirmed priority.
Do not expand into all PGL orbits, all intersections, all nilpotent thickness,
precise nil/H0/T2, general Pn classification, Lean formalization, or English
translation unless indispensable to the stated target. Do not publish externally
or submit to a journal without a separate user instruction.

## Safe checkpoints, reproducibility and genuine hard limits

Reuse the environment. No new GitHub Actions. Check relevant identities/tests,
rebuild changed TeX and inspect rendered pages. Test success is not proof certification.
Use one Git integrator, inspect main and preserve concurrent changes. Commit/push
non-force, verify the remote head, then CONTINUE. No automatic resets or force-push.
The user permits current research bundles and Library saves without confirmation.
Package the current necessary proof/input/certificate dependency closure, not old
ZIPs, histories, caches or redundant drafts. Never delete existing Library files.
Generation, validation, main persistence and Library persistence are separate states.

The only normal stop is the complete terminal contract. Genuine exceptions are
a user stop or an actually reached runtime/resource/access/input limit that
prevents useful authorized work. Try safe alternatives and continue unblocked
requirements first; never bypass controls. Before an unavoidable interruption,
save the exact remaining obligations, reproducible blocker, current commit and
next concrete operation in CURRENT. Say incomplete, not done. Do not imply
background continuation or automatic future resumption. These instructions do
not enlarge the platform's execution limits or permissions.

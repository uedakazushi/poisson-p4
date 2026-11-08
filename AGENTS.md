# Poisson P4 — global reconstruction, not a local-lemma sprint

2026-11-09 / v67-global-reconstruction

## Read first

RESEARCH_POLICY.md → research/conceptual/GLOBAL_PLAN.md →
papers/poisson_divisor_symmetries_ja.tex → research/conceptual/NOTES.md.
Then inspect the relevant original proof and the main classification manuscript.
The user explicitly requests deeper conceptual organization of the WHOLE
classification and main updates, including any new mathematics. This supersedes
the previous priority of proving an H+2Q Koszul lemma and the old stop-after-147
instruction. classification/PLAN.json records construction, not the active queue.

## Primary output

A proof architecture that replaces repeated recognition, reconstruction, and
closure arguments by reusable theorems. The axes are:
1. recognize all structures via marked geometric data;
2. reconstruct full sources, including all exceptional parameters;
3. lift divisor deformations to actual Poisson families;
4. count maximal distinct closed images, not normal-form labels or all orbits.

The current manuscript/catalogue asserts 147 reduced components (141 rank four,
6 rank two), before PGL5 quotient. It is not an axiom. Preserve valid results;
correct affected claims if an actual defect is found. Do not say the entire
classification was independently verified just because existing receipts say PASS.

## Global entrances and available proved tools

A: actual degeneracy divisor has a hyperplane, including 5H.
B: no hyperplane and reducible, hence an irreducible quadric plus cubic.
C: reduced integral quintic, normal versus nonnormal conductor geometry.
R2: saturated rank-two foliation plus its determinant section/zero divisor.
These are pointwise entrances, NOT disjoint sets of irreducible components.

The structural paper proves the marked-hyperplane affine quadratic equivalence,
mult_H(D)=5-deg(Pf P), and fixed-f,Z reconstruction:
    div P=Z, P#df=fZ, Pf(P)=f, f != 0.
The first two are affine-linear in P; the last is quadratic. Jacobi follows.
It also gives a common resonance decomposition, relative inverse-form model,
first/second-order lifting equations, graph-closure/Rees interpretation, finite
factor markings, and the rank-two foliation/zero-divisor reconstruction.
Novelty is unconfirmed. These are not already a full classification of f,Z or
all relative fibres. Check proofs and application hypotheses before using them.

## First concrete assignment

Follow GLOBAL_PLAN P0. Make one compact mapping from existing component IDs and
main theorems to A/B/C/R2, using actual generic divisors and proof sources, not
names guessed from IDs. Do not spend the whole task building another ledger.
Then integrate ONE common entrance/reconstruction proof into the main manuscript
that covers BOTH simple H+quartic AND multiple-H sources, checking a third
marked-H configuration for reuse. Read the full old source equations, inverse
maps, and zero/unit complements. Identify exactly which old lemmas or repeated
arguments are replaced and which finite calculations remain indispensable.
The first task must not end with another plan, audit receipt, or renamed table.

Next apply the moving-divisor lifting equations to the existing TRIPLE-to-Ga9
family and one distinct splitting/collision mechanism: explain the direction,
necessary correction terms and remaining obstruction, not just replay a PASS.
Then unify B via its canonical symmetry in aut(Q), C via conductor/normalization/
residue/conormal data, and connect the rank-two support theorem. See GLOBAL_PLAN.

## Non-negotiable mathematical distinctions

- f != 0 is essential; the paper contains a rank-two counterexample.
- Fixing f,Z is not classifying them. A linear presentation alone is not progress.
- Simple H: deg f=4. Multiple H: deg f<=3. Pf(P_2)=0 belongs only to the latter.
- r=1 includes 5H as well as an integral quintic. Canonical symmetry dimension
  zero does not imply absence of all symmetries.
- Preserve translation cokernels, resonance walls, and every zero complement.
- Generic fibres, geometric generic fibres, their Galois orbits, and special-base
  components differ. Never discard vertical rank-four components.
- Divisor deformation need not lift to a Poisson deformation. Nonlinear auxiliary
  coordinate changes do not establish projective equivalence.
- Rank-three H+2Q has excess dimensions 1 and 6; uniform affine-linear primitives
  are false. That local project is optional, not the global prerequisite.
- MPS40 classifies components meeting SNC, not all components.
- A normal form, source dimension, or generic divisor type does not establish
  ambient componenthood, all-point coverage, or closure noncontainment.

## Scope, progress, and integration

Do not pursue full orbit classification, all component intersections, complete
nonreduced local rings, precise nil/H0/T2, or unlimited new component searches.
Do not demand a huge Rees-ideal Groebner computation. Use the weakest necessary
support statement. After two nonclosing attempts, diagnose and change the route.
A mathematical milestone replaces actual repeated proofs and explains a family,
dimension, exceptional stratum, or boundary. More tests/pages/receipts alone do not.
Keep all valuable existing proofs in current manuscripts; do not replace the
368-page proof by summaries or a link to an old version. Split papers only after
clarifying dependency and exposition. Use Japanese LuaLaTeX, primary references,
explicit novelty and verification status. No writable boxes.

## Validation and safe persistence

Reuse the environment; no GitHub Actions or reinstall-everything workflows.
The scoped identity check is:
    python code/check_global_structure.py --output work/global-structure.json
It does NOT verify all 147 components. Rebuild changed TeX twice and inspect PDFs.
Use one integrator, inspect main, preserve others' changes, and commit/push
non-force with an expected-head check. Synchronize README/policy/current/next.
The user authorizes current self-contained research bundles and Library saves
without confirmation. Include necessary current proofs and inputs; exclude old
ZIPs, logs, caches, backups, and nested archives. Do not delete Library records.
Generation, checking, main update, and Library upload are separate states. Report
only actual outcomes. If Library access is unavailable, create/validate the bundle
and state that save was not performed. Never claim future/background execution.

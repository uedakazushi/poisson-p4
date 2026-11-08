# Poisson P4 — conceptual structure first

2026-11-09 / v66-conceptual-structure

Read RESEARCH_POLICY.md, research/conceptual/NOTES.md,
research/classification/CURRENT.md, then research/NEXT_TASKS.md.
The user's latest direction is to organize the classification conceptually,
not merely verify more computations. This supersedes verification-first as the
primary queue. Verification remains necessary and supports the structural work.

## Primary outcome

Reconstruct the paper around geometry, reconstruction, and deformation:
(1) the actual degeneracy divisor and its canonical commuting symmetries;
(2) closed rational forms with bounded poles on a fixed divisor;
(3) degeneration and noncontainment of the resulting families.
Do not replace this by new audit receipts or renamed case tables. A useful
structural theorem should actually replace repeated calculations or explain
families, dimensions, or boundaries.

## Current proved structure

papers/poisson_divisor_symmetries_ja.tex contains full written proofs:
- r reduced irreducible factors give a canonical (r-1)-dimensional abelian
  subalgebra of global Poisson symmetries on P^(2n), including nonreduced D;
- on P4, projective orbit dimension is at most 25-r;
- fixing the full anticanonical section F turns the rank4 problem into the
  linear kernel V_F of F d(alpha)-dF wedge alpha and Pfaffian quadrics;
- the two-factor symmetry recovers the affine modular field.
Novelty is NOT established. These results do NOT by themselves prove the 147
count, all componenthood statements, or all scheme-theoretic deformation claims.
Do not confuse global vector-space dimension of symmetries with pointwise rank.

## Next concrete work

Replace the five rank(Q)>=4 H+2Q linear-kernel computations by a
filtered de Rham/Koszul argument, but NEVER assert affine-linear primitive
exactness for all rank(Q)>=3. Rank-three transverse has an additional
1-dimensional class, and tangent has six filtered classes. See
research/conceptual/QUADRATIC_POLE_ORDER.md. Keep degree bounds, genuine
cohomology versus filtration obstructions, and degenerate d(beta) cases distinct. Then unify
the two-factor semisimple/nilpotent/resonance entrances using the canonical
symmetry. Keep irreducible quintics and rank2 foliations as separate structural
inputs. Track actual zero/unit complements and source coverage throughout.

The existing 147-component manuscript, catalogue, certificates and valid results
remain intact. Do not shorten the 368-page proof by deleting mathematical
content. Integrate structural results into a comprehensible order only after
their applications are justified. The old classification PLAN.json records
completed construction; publication/REVIEW.md is a supporting audit queue.
The existing TRIPLE-to-Ga9 verifier remains available; do not rerun old suites
without a specific mathematical need.

## Workflow

Use Japanese LuaLaTeX, existing environments, and finite relevant calculations.
No GitHub Actions. Do not reopen precise nil/H0/T2 or full orbit classification
without an indispensable dependency. Preserve immutable proof inputs.
Check current main, use non-force/expected-head updates and preserve other work.
The user authorizes current self-contained research bundles and Library saves
without asking. Distinguish generation, validation, GitHub update, and Library
save. Never delete existing Library records without explicit instruction.
A conceptual subproject package is not a replacement for the full classification
bundle; keep its scope and any missing global assets explicit.

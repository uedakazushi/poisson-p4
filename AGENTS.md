# Poisson P4 — independent referee pass before public release

2026-11-09 / v74-publication-candidate

## Read order and current assignment

Read RESEARCH_POLICY.md, research/publication/REFEREE_PLAN.md,
research/classification/CURRENT.md, research/conceptual/FINAL_TARGET.md,
then the current main manuscript and original proof sources.

The specified high-risk proof review and finite runtime-path review are now
closed at the scope stated in REFEREE_REPORT.md. CURRENT.md records the remaining
delivery state. Finish the current full-bundle validation and main/Library save;
do not restart completed derivations, add speculative lemmas, or rerun all old
experiments. Reopen a proof only for new evidence or a concrete defect.

The 147-component classification has been DECLARED COMPLETE, with a 418-page
main manuscript, three companion PDFs and five verified Library ZIP volumes.
This is the starting point, not a command to stop or to believe every claim.
The new phase is an ADVERSARIAL, theorem-level, independent referee review
and final publication editorial pass. It supersedes the previous instruction
to keep inventing conceptual lemmas or to stop when old M1–S1 flags say PASS.

## Proof criteria applied in this review — retain their scope

1. First challenge the reduced-fixed-divisor theorem:
   papers/normalization_residue_reconstruction.tex,
   lem:quintic-ambient-lift, thm:reduced-fibre-linear,
   cor:reduced-fibre-boundary, cor:reduced-divisor-component-images.
   Derive the theorem afresh and attempt counterexamples for singular/reducible
   D, q_F=0, kernel rank jump, generically irreducible fibres and divisor-image
   comparisons. Clearly distinguish all-Poisson compatibility from rank-four
   nondegeneracy. Do not accept the PDF paragraph or a PASS receipt as proof.

2. Independently audit papers/rank_two_support_reconstruction.tex:
   stable degree-two foliation -> full Rat113/Rat122 sources, five-section and
   spinor arguments, nonstable/divisorial cases, six ambient rank-two components
   versus rank-four boundaries. Check the actual assumptions of Correa–Muniz
   and Loray–Pereira–Touzet. Show the rank-four support is completed BEFORE the
   rank-two closure term is used; no circularity.

3. Independently audit papers/maximal_image_synthesis.tex and
   source:last-five-explicit in papers/component_family_catalogue.tex.
   Check all five families including all six magnetic parameters, zero/unit,
   proper flags, the NEW25, T36/T38 and T37 separating conditions, the actual
   noncontainment relative to all older 142, and full pointwise source coverage.
   A dimension or generic-divisor comparison is not enough. Try explicit
   degeneration/collision counterexamples. If found, repair all dependent claims.

4. Verify that all 51 serialized rational circuit families really define
   irreducible DENSE generic sources with UNIVERSAL Jacobi, not just one
   randomly chosen exact rational sample. Match the MPS40 theorem's precise
   normal-crossings scope and all source-to-S61 coordinates. Check denominators,
   nonempty opens, roots, branch limits and the eleven additional families.
   Choose proof-level cross-checks, not 51 expensive symbolic explosions.

5. Independently trace the minimal proof-dependency cut set from all-rank-four
   pointwise support through the six rank-two sources to the 147 maximal
   noncontained images. Use the claim graph as a locator, not mathematical proof.
   Re-run selected exact finite certificates with independent implementations
   when feasible; do not indiscriminately replay all historical tests.

6. Verify reproducibility using the complete current repository sources and
   locally generated or already available ZIP volumes. Check input hashes,
   references, essential proof assets, and rebuild the four current Japanese
   LuaLaTeX PDFs. **Do not fetch or download binary ZIP files from Library.**
   Remote Library file readback is not a condition of referee acceptance. Check the current proof layout, introductions and primary
   literature comparisons for a human referee, especially the 40 prior SNC
   components vs the alleged 147 total. Do not publish or submit externally.

## Deliverables and standards

Maintain ONE current concise research/publication/REFEREE_REPORT.md with
each audited theorem's exact hypothesis, proof location, independent reasoning,
failed counterexample attempts or actual counterexample, and downstream
dependencies. Do not write a large worklog, copy old receipts, or mark PASS
only because scripts succeeded. Distinguish: written proof / independently
read and checked / exactly recomputed / sample only / still open.

Correct actual mathematical gaps in the main theorem, companion proofs,
catalogue, research policy, Introduction, and necessary certificates; do not
preserve 147 by arbitrary assumptions. Keep all valuable existing results and
proofs, and avoid making 418 pages shorter by discarding needed mathematics.
A new general theorem is welcome ONLY if it removes a real proof dependency;
do not restart speculative component hunting or full orbit classification.
After hard proof gates close, finalize a reader-ready manuscript (split into
logical papers only if it improves proof clarity), update README/POLICY/STATE,
build a minimal current reproducible bundle, validate and save it.

The audit PASS criterion is proof-level closure of the high-risk claims,
plus clean-room reproducibility and accurate literature comparison. There is
no requirement to independently rerun every one of 147 original computations.
If material uncertainty remains, state the precise blocker and do not claim
publication readiness.

Use current existing environments, no GitHub Actions, no redundant installers,
no bulk new old-version archives. Commit/push verified milestones non-force to
main after checking head and preserving parallel changes. Library save is
already authorized; do not delete older Library files. Report separately
generation, mathematical checking, PDF validation, GitHub persistence and
Library persistence. External hard limits may interrupt work, but do not claim
asynchronous continuation or automatic future resumption.

# CURRENT OVERRIDING RULE: GITHUB MAIN IS THE ONLY SOURCE OF TRUTH (2026-11-10)

These rules override every older Library persistence/packaging instruction in this repository and apply to all modules, including the 147-component classification and Poisson subvariety incidence.

* **GitHub main is canonical.** Every new proof, current Japanese LuaLaTeX manuscript, active research policy/notes, exact input, certificate, script and build recipe must be committed to main before a milestone is called saved. Neither Library, a previous ZIP, an untracked work/ directory, nor a previous conversation may be required to resume or reproduce current work.
* **NO GITHUB ACTIONS.** Never create, configure, dispatch, or run GitHub Actions or hosted CI workflows. Use manual local/Cloud Codex shell scripts only. Do not install actions/workflows as a convenience.
* **Library is emergency-only.** Do not upload/save to ChatGPT Library on the normal path. Only when an actual GitHub push fails may the current small self-contained bundle be saved there as emergency fallback (if possible). Clearly mark GitHub push failure and restore main as canonical at the next opportunity. Do not delete old Library assets without explicit permission.
* **Pre-push and post-push requirements:** prove or scope-check mathematical changes; run the relevant code/TeX manually; check complete tracked-input closure; ensure no essential source lives only in work/ or Library; non-force commit/push with expected remote head; re-read main and verify its new SHA. GitHub write success, mathematical verification, PDF generation, and Library fallback are distinct outcomes.
* Read RESEARCH_POLICY.md -> REPRODUCE_GITHUB.md -> current module research policy -> notes -> paper. Keep current documents self-contained, preserve existing valid results and proof certificates, avoid histories and duplicate bundles. Record unverified external claims/newness as unverified.
* Manual reproduction: bash scripts/verify_github_source.sh incidence ; for classification preflight and current papers: bash scripts/verify_github_source.sh full . No Library retrieval and no Actions.

---
# Poisson P4 — independent referee pass before public release

2026-11-09 / v74-publication-candidate

## Read order and current assignment

Read RESEARCH_POLICY.md, research/publication/REFEREE_PLAN.md,
research/classification/CURRENT.md, research/conceptual/FINAL_TARGET.md,
then the current main manuscript and original proof sources.

The specified high-risk proof review and finite runtime-path review are now
closed at the scope stated in REFEREE_REPORT.md. The current full bundle has
been validated, saved to Library, and its source pushed to main. CURRENT.md and
reports/publication_candidate_delivery.json record completion. No mandatory
work remains for this goal. Do not restart completed derivations, add speculative
lemmas, or rerun all old experiments. Reopen a proof only for new evidence or a
concrete defect; a new research or external submission scope needs user direction.

The 147-component classification has been DECLARED COMPLETE, with a 418-page
main manuscript, three companion PDFs and five verified Library ZIP volumes.
The subsequent ADVERSARIAL, theorem-level, independent referee review
and final publication editorial pass are now complete at the specified scope.
Their retained criteria below supersede the previous instruction
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
main after checking head and preserving parallel changes. Library is only an emergency fallback after an actual failed GitHub push. Do not delete old Library files. Report mathematical, PDF and verified GitHub-main persistence separately. External hard limits may interrupt work, but do not claim
asynchronous continuation or automatic future resumption.

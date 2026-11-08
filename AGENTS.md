# Poisson P4 — publication review

2026-11-09 / v65-publication-review

## Entry and current authority

Read RESEARCH_POLICY.md, research/classification/CURRENT.md,
research/publication/REVIEW.md, then research/NEXT_TASKS.md.
The user explicitly authorized publication-focused independent review and main
updates on 2026-11-09. This is the new instruction that supersedes the previous
stop-after-classification instruction. Do not interpret completed construction
items in research/classification/PLAN.json as a prohibition on this review.

## Mathematical status

The manuscript/catalogue assert exactly 147 reduced irreducible components of
the projectivized Poisson scheme on complex P4 before PGL5 quotient: rank four
141, rank two 6. Do not change that count without a mathematical reason. The
new review has NOT independently verified the entire classification. Existing
repository review receipts are not external peer review or formal proofs.
The main manuscript and 147-component catalogue remain unchanged in this pass.

## Work now

Follow the priority queue in research/publication/REVIEW.md. Inspect original
proofs rather than repeatedly composing PASS receipts. Track each theorem's
hypotheses, actual source, inverse maps, denominator/zero complements, and
conclusion. Distinguish full source coverage from dense-source containment,
componenthood from family dimension, and global components from local branches.
Prioritize normal/nonnormal quintic coverage, stable rank-two recognition, and
T37 dense-orbit/residue/leaf arguments. No new component search, precise nil
indices, H0/T2, or full orbit classification without a specific dependency.

## Checked in the present review

code/verify_final_bridge.py independently rebuilds the TRIPLE-to-Ga9 family
from a product bracket and the Ga9 potential. It does not import old producers.
reports/final_bridge.json records 81 exact checks. This verifies the displayed
family, not the coverage of all original TRIPLE sources or Ga9 componenthood.
research/publication/verification_note_ja.tex supplies explicit proofs of this
bridge and the finite-cover/rank-closure/specialization mechanisms.
Their application-specific geometric inputs remain separately reviewable.

## Validation and editing

Use the existing Python/SymPy environment:

    python code/verify_final_bridge.py --output work/publication-review/final_bridge.json
    python -m unittest discover -s tests -p test_final_bridge.py

classification.py check/next and the old tests are control-plane checks only;
a completed queue is not a mathematical re-verification. Do not reinstall the
whole environment or rerun every legacy experiment. No GitHub Actions.
Write papers in Japanese LuaLaTeX. Preserve valid mathematical content. Do not
append endless historical addenda or retain contradictory current-status text.
Do not replace the 368-page proof by a short summary. Paper splitting comes
after the dependency review, not before it.

## Save safely

Inspect the current main head and changes. Use ordinary non-force commits to
main, or equivalent Git data API tree/commit/ref updates with an expected-head
lease. Rebase/reconcile if main advanced; do not overwrite another session.
Keep original proof/certificate inputs unchanged. Generate, validate, push, and
Library upload are distinct states. Only report actual successful operations.
The user permits current self-contained research bundles and Library saves
without asking. Exclude old ZIPs, logs, intermediate build products and caches;
never delete existing Library records without explicit instruction. A focused
review package missing the main proof is not a full research-continuation bundle.

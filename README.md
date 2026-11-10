# Poisson structures on projective four-space

2026-11-10 / v16-local-geometry + v74-classification.

**GitHub main is canonical.** Current Japanese LuaLaTeX sources, exact code, inputs and certificates suffice to restart the research. Generated PDFs are not tracked. No GitHub Actions. Library is used only after an actual failed GitHub push; no Library download or old ZIP is required.

## Current research

Start with [RESEARCH_POLICY.md](RESEARCH_POLICY.md), then the [module policy](research/subvariety_incidence/RESEARCH_POLICY.md) and [current notes](research/subvariety_incidence/NOTES.md).

The new standalone manuscript [isotropic_log_tangent_cone_ja.tex](papers/isotropic_log_tangent_cone_ja.tex) proves the complete seven quadratic initial equations at the isotropic Log point, a nonzero square-zero degree-two class in the **true tangent cone**, local dimension 28, Hilbert–Samuel multiplicity 1, **smoothness of the 28-dimensional logarithmic component at the intersection** and the 23–24 bound for its intersection with the 25-dimensional rank-four component, and cancellation of the mixed ac obstruction in the full ambient second-order deformation. Ambient reducedness and the full higher initial ideal remain open. The 26-/24-dimensional components of the quadratic approximation are NOT counted as new Poisson components.

The explicit common three-parameter Poisson family is independently verified by [check_isotropic_intersection_dim23.py](code/subvariety_incidence/check_isotropic_intersection_dim23.py). Its logarithmic epsilon specialization and affine-linear long-source rectification give an actual 23-dimensional family in the component intersection. The numerical lower bound 23 also follows from the ordinary dimension inequality (28+25-30); the explicit geometry is the added result. The exact intersection dimension (23 or 24) remains open.

The [universal invariant-quadratic pencil](papers/isotropic_nilpotent_pencil.tex) proves a strict rank-four smoothing for a five-parameter smooth quadratic slice preserved by a nilpotent symmetry. Its rank-two fibres belong to Log1112 and its rank-four fibres to the long nonconstant-g 25-dimensional component, via an actual affine-linear source map. An explicit PGL5 shear reduces the five quadratic coefficients to a single residual parameter eta/(1-delta*eta); this does not determine whether the entire component intersection has dimension 23 or 24. The [exact producer](code/subvariety_incidence/check_universal_nilpotent_pencil.py) and [tracked certificate](research/subvariety_incidence/certificates/universal_nilpotent_pencil.json) are included in the manual tangent-cone verification.

An explicit [attracting PGL5 cocharacter](papers/isotropic_nilpotent_pencil.tex), diag(1,1,rho,1,rho^-1), sends the normalized nilpotent pencil parameters (eta,lambda) to (rho² eta,rho lambda); its rho→0 limit is a rank-two fully isotropic Log point. Invariant quadratic forms have precisely two rank strata, rank(q|W)=0/1, under the subgroup preserving A and the base three-form line. This does not classify full PGL orbits or resolve the global 23/24 component intersection. The [exact cocharacter producer](code/subvariety_incidence/check_nilpotent_torus_attraction.py) and [certificate](research/subvariety_incidence/certificates/nilpotent_torus_attraction.json) are part of manual reproduction.

The [torus-weight decomposition](papers/isotropic_nilpotent_pencil.tex) gives exact 30-dimensional tangent multiplicities (1,5,12,10,2) in weights -2 through 2, and 95 normal-direction multiplicities (2,7,17,26,28,12,3) in weights -2 through 4. All seven quadratic initial relations are homogeneous; the nonzero square-zero tangent-cone class has weight 2. This enables weight-block elimination for the still-open cubic initial equations. The [symbolic producer](code/subvariety_incidence/check_nilpotent_tangent_weights.py) and [tracked certificate](research/subvariety_incidence/certificates/nilpotent_tangent_weights.json) are self-contained.

The [incidence companion](papers/poisson_subvariety_incidence_ja.tex), [local dimension proof](papers/quadric_log_component_intersection.tex) and [long rank-four membership proof](papers/rank_four_long_component_membership.tex) retain generic quadric extensions, cubic/del Pezzo rigidity, the explicit 25-dimensional rank-four component through the 28-dimensional logarithmic component, and the non-Cohen–Macaulay consequence at their stated dependencies.

## Manual reproduction

```bash
bash scripts/verify_github_source.sh tangent-cone  # new result and 6-page PDF
bash scripts/verify_github_source.sh incidence    # existing incidence plus new result
bash scripts/verify_github_source.sh full         # classification preflight plus incidence
latexmk papers/isotropic_log_tangent_cone_ja.tex
```

The new producer is [check_isotropic_quadratic_cone.py](code/subvariety_incidence/check_isotropic_quadratic_cone.py); its [exact certificate](research/subvariety_incidence/certificates/isotropic_quadratic_cone.json) is regenerated without historical inputs. Full instructions: [REPRODUCE_GITHUB.md](REPRODUCE_GITHUB.md).

## Existing 147-component classification

The current classification manuscript claims 147 reduced projective Poisson components before the PGL5 quotient: 141 generically rank four and six rank two. Original mathematical content and verification assets are retained, not re-audited by the new local calculation.

Entry points: [classification state](research/classification/CURRENT.md), [referee report](research/publication/REFEREE_REPORT.md), [referee plan](research/publication/REFEREE_PLAN.md), [target theorem](research/conceptual/FINAL_TARGET.md), [main manuscript](papers/poisson_p4_classification_ja.tex), [family catalogue](papers/component_family_catalogue.tex), and [classification reproduction](research/reproduction/README.md). Previous delivery receipts remain at their stated scope; they do not establish external acceptance or priority. Do not publish or submit externally without explicit authorization.

Codex instructions: [AGENTS.md](AGENTS.md).

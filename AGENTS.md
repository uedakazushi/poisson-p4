# Poisson P4 — current execution instructions

2026-11-10 / v14. Read RESEARCH_POLICY.md first, then research/subvariety_incidence/RESEARCH_POLICY.md, NOTES.md, and papers/isotropic_log_tangent_cone_ja.tex. Japanese LuaLaTeX is the manuscript format.

## Current active priority (v14)

The 28-dimensional logarithmic component and 25-dimensional long nonconstant-g rank-four component have an explicit common 23-dimensional family at the isotropic Log point. The universal strict four-parameter Poisson family has a logarithmic epsilon-source and an actual affine-linear long-source identification, with all six brackets exactly checked. The GL5+3 derivative has rank24 and nonzero minor 10871635968; Euler gauge plus family rank48 leaves projective rank23. The numerical lower bound 23 also follows from the standard dimension theorem; the new content is the explicit source family. Exact intersection dimension 23 versus 24 is open, as are the full higher initial ideal and ambient reducedness. Read papers/isotropic_log_tangent_cone_ja.tex and code/subvariety_incidence/check_isotropic_intersection_dim23.py. No Actions, routine Library copies, or generated PDF tracking.

## Mathematical status

The full quadratic initial ideal at the isotropic Log point is now rigorously computed: seven relations on a 30-dimensional tangent space. With s_i=t_(20+i), a=t27, b=t28, c=t29, it is
(s0*c-s2*a,s1*c,c*(s2+c),s3*c,s4*c+s6*a,s5*a,s5*c).
The nonzero homogeneous degree-two class h=a*(s2+c) survives in the TRUE tangent cone and squares to zero, since ALL quadratic initial relations were computed. Thus the tangent cone is nonreduced. Its unique top-dimensional support is the plane a=c=0, generically length one; local dimension is 28 and Hilbert–Samuel multiplicity is 1. The logarithmic source minor is independently rechecked. The mixed ac obstruction can be removed to second order in the full ambient deformation by an explicit correction; the pure c^2 obstruction remains.

Do NOT infer nonreducedness of the original ambient local ring. The radical of the quadratic obstruction scheme has dimensions 28/26/24, but the lower two are not yet actual tangent-cone or Poisson components. Higher initial equations and the true support of h remain open.

Previous valid results, including the explicit PGL5 membership in R4_MODULAR_RANK1_LONG_NONCONSTANT_25, the 28/25-dimensional component intersection, non-CM/depth<=25 with their component-theorem dependencies, generic quadric incidence and the original 147-component classification assets, remain in their current companion papers. Do not discard them or count quadratic-model components as new Poisson components.

## Next work

Compute cubic and higher initial equations by solving the 95 normal coordinates to order two and combining the residuals with the quadratic syzygies. Identify the lower-dimensional true tangent cone, the annihilator/support of h, and intersection with the known long component. Treat ambient reducedness and all-order lifting as separate unsolved problems. Change direction if mathematical value favors a more conceptual result.

## Verification and persistence

- GitHub main is the only source of truth. Non-force update with expected head and re-read the remote SHA. No GitHub Actions or hosted CI. Library only as emergency fallback after an actual failed push, never the normal path. Do not delete existing Library assets without permission.
- Track current TeX, code, exact inputs, small certificates and manual reproduction instructions. PDFs and compile products are generated, not tracked. No history bundles, old ZIPs, work logs or caches in current deliverables.
- New independent verification: bash scripts/verify_github_source.sh tangent-cone. It rebuilds all 150 linear columns, all 465 tangent pairs, the logarithmic minor, mixed correction and a six-page standalone manuscript. Old incidence checks remain available via incidence; full also runs the existing classification preflight.
- Separate proof in manuscript, independent reasoning, exact computation, PDF build, remote persistence and unconfirmed novelty. New code must run before its certificate is called PASS.
- Existing classification referee criteria and scope are in research/publication/REFEREE_PLAN.md, REFEREE_REPORT.md and research/classification/CURRENT.md. Do not reopen completed computations without new evidence; do not treat prior PASS reports as substitutes for checking a newly used hypothesis. No external publication/submission without a separate user instruction.

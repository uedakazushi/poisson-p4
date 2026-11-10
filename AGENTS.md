# Poisson P4 — current execution instructions

2026-11-10 / v16. Read RESEARCH_POLICY.md first, then research/subvariety_incidence/RESEARCH_POLICY.md, NOTES.md, and papers/isotropic_log_tangent_cone_ja.tex. Public manuscripts use Japanese LuaLaTeX.

## Mathematical status and priorities

At the totally isotropic Log point p, the projective Poisson tangent has dimension 30 and the complete degree-two initial ideal is
(s0*c-s2*a, s1*c, c*(s2+c), s3*c, s4*c+s6*a, s5*a, s5*c)
with s_i=t_(20+i), a=t27, b=t28, c=t29. h=a*(s2+c) survives as a nonzero square-zero element of the TRUE tangent cone. Local dimension 28, multiplicity 1; the unique 28-dimensional Log1112 component is smooth. The 25-dimensional long nonconstant-g rank-four component meets it in local dimension 23 or 24, with an explicit common 23-dimensional Poisson family; the exact dimension is open. The ambient local ring's reducedness is not proven. Keep the separate hypotheses supporting the full 147-component classification and long-component membership.

The universal nilpotent-invariant quadratic Poisson pencil has a genuine PGL5 cocharacter g_rho=diag(1,1,rho,1,rho^-1). In normalized coordinates its projective action is (eta,lambda)->(rho^2 eta,rho lambda), contracting all members to the rank-two fully isotropic model. Quadratic forms preserving the fixed nilpotent A and base three-form line have exactly two rank strata, rank(q|W)=0 or 1. The ratio eta/lambda^2 for lambda!=0 is only a cocharacter invariant, NOT a proved full-PGL invariant. Proof: papers/isotropic_nilpotent_pencil.tex; independent script code/subvariety_incidence/check_nilpotent_torus_attraction.py.

**Current priority: weight-graded cubic initial equations.** At p the tangent weights (-2,-1,0,1,2) have dimensions (1,5,12,10,2). The linearized Jacobi normal directions have weights (-2,-1,0,1,2,3,4) with dimensions (2,7,17,26,28,12,3); all seven quadratic initial equations have weights (2,1,2,3,1,1,1) and nilpotent h weight 2. Exact 30-eigenvector proof and script: papers/isotropic_nilpotent_pencil.tex, code/subvariety_incidence/check_nilpotent_tangent_weights.py. Solve the 95 normal variables to quadratic order by these seven weight blocks, then compute cubic residuals modulo syzygies of the seven quadratic equations. Identify the actual low-dimensional tangent-cone support, the annihilator of h, and whether the global Log/long intersection has dimension 23 or 24. Explore symmetry-breaking deformations; avoid restarting fixed-A quadratic casework.

Do NOT infer ambient nonreducedness from its nonreduced tangent cone. Do NOT label the 26-/24-dimensional quadratic-obstruction components as real Poisson components. Distinguish first- and second-order lifting from all-order lifting, and manuscript claims from independently checked computations. Newness and external referee review remain unconfirmed.

## Reproduction and persistence

- GitHub main is the only source of truth. Preserve existing current proofs and certificates. Update non-force with an expected remote head, then re-read main. GitHub Actions and hosted CI are forbidden.
- ChatGPT Library is emergency-only after a genuinely failed GitHub push; no routine Library upload, no historical ZIP bundles or history logs, no deletion of Library assets without explicit instruction.
- Track Japanese TeX source, current policy/notes, exact code, minimal input and small certificates. Generated PDFs, auxiliary/cache files and redundant previous versions stay out of Git.
- Manual current verification: bash scripts/verify_github_source.sh tangent-cone; existing incidence and full modes retain their recorded scope. New torus and tangent-weight producers must match their exact JSON certificates. Do not claim fresh-clone/full-classification verification without actually executing it.
- Proofs and theorem assumptions in papers/isotropic_log_tangent_cone_ja.tex and its shared TeX; broader incidence in papers/poisson_subvariety_incidence_ja.tex. Existing classification referee scope: research/publication/REFEREE_PLAN.md, REFEREE_REPORT.md, research/classification/CURRENT.md. Do not publish/submit externally without explicit user instruction.

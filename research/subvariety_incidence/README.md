# Poisson subvariety incidence — current v8

GitHub main is the only authoritative research source. Read [research policy](RESEARCH_POLICY.md), [notes](NOTES.md), and [the Japanese LuaLaTeX paper](../../papers/poisson_subvariety_incidence_ja.tex). Shared proofs are in papers/quadric_log_extension.tex and its nested sources: quadric_log_orthogonality.tex, quadric_log_isotropic_boundary.tex, and quadric_ambient_quadratic_obstruction.tex.

**New theorem:** At the special totally isotropic Log point with fixed-fibre extension Spec C[a,b,c]/(ac,c²), the embedded c tangent is obstructed already at order two in the FULL ambient P4 projective Poisson scheme. A universal five-coefficient functional annihilates all 150 ambient linearized directions but evaluates to -9 on Jac(K). The 350×150 projective linearization has rank95 (augmented rank96), so the ambient projective tangent dimension is30. This shows singularity, NOT ambient nonreducedness.

**Manual exact reproduction** (9 independent SymPy producer scripts and 9 checked JSON certificates; no generated PDFs tracked):

```bash
bash scripts/verify_github_source.sh incidence
```

The original general Q3 Rat/Log/Aff results, resonant orthogonality, smooth isotropic examples, cubic and del Pezzo geometry are retained. Full special Q3 fibre classification, ambient nonreducedness and component membership remain unresolved. GitHub Actions forbidden; Library is emergency-only after an actual failed push.

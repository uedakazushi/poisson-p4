# Poisson subvariety incidence — current v7

GitHub main is the only canonical source. Read [research policy](RESEARCH_POLICY.md), [notes](NOTES.md), and [Japanese LuaLaTeX manuscript](../../papers/poisson_subvariety_incidence_ja.tex). Shared proofs: [general Log](../../papers/quadric_log_extension.tex), [nondegenerate resonant blocks](../../papers/quadric_log_orthogonality.tex), [isotropic restriction boundary](../../papers/quadric_log_isotropic_boundary.tex).

New result: for three explicit smooth quadrics with degenerate q|W, full normalized ambient Poisson extension schemes are respectively A², two reduced transverse lines, and nonreduced Spec C[a,b,c]/(ac,c²). The latter has a rank-four dense open of its reduced A² and an embedded nilpotent direction on a line. rank(q|W) alone cannot classify extension fibres. This is NOT a classification of all special Log points or all q|W degenerations.

Eight exact SymPy producers and eight matched JSON certificates are tracked; run:

```bash
bash scripts/verify_github_source.sh incidence
```

No tracked PDFs, no GitHub Actions, and no routine Library use.

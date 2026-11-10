# Poisson subvariety incidence — current v5 research

Source of truth: GitHub main. Read [policy](RESEARCH_POLICY.md), [notes](NOTES.md) and [the main Japanese LuaLaTeX paper](../../papers/poisson_subvariety_incidence_ja.tex).

The full common Log extension proof is [quadric_log_extension.tex](../../papers/quadric_log_extension.tex), which includes [quadric_log_orthogonality.tex](../../papers/quadric_log_orthogonality.tex).

**New result:** For resonant Log structures with nondegenerate H (2-dimensional weight-one form) and nondegenerate qU, the full ambient extension fibre is A1 (rank 2) when the cross-term b is nonzero; precisely when b=0 the fibre is Spec C[t,v]/(v(t+2)) with a rank-four branch. The orthogonality locus is codimension two in this quadric family, not merely the previous one-parameter wall. The case c=0 is included if qU is nondegenerate. Degenerate H or qU remains open.

Six tracked SymPy scripts reconstruct six exact JSON certificates. The general classification theorem is proved geometrically, not inferred from the five sample tests.

```bash
bash scripts/verify_github_source.sh incidence
```

No PDF tracking, no GitHub Actions, no routine Library storage; all sources and proof inputs remain on GitHub main.

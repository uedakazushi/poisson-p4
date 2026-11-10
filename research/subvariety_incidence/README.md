# Poisson subvariety incidence — current v6

GitHub main is canonical. Start with [research policy](RESEARCH_POLICY.md), [notes](NOTES.md) and the [Japanese LuaLaTeX manuscript](../../papers/poisson_subvariety_incidence_ja.tex), which includes [Log fibre proofs](../../papers/quadric_log_extension.tex) and [the invariant orthogonality classification](../../papers/quadric_log_orthogonality.tex).

New theorem: **the entire family of resonant Log Poisson extensions with nondegenerate 3-variable quadratic qU is classified**. For rank(H)=2, a rank-four branch exists exactly when the cross-term b vanishes (fibre Spec C[t,v]/(v(t+2))). For rank(H)=1, a rank-four branch exists exactly when c vanishes in the adapted normal form qU=u0²+2u1z+cz² (fibre Spec C[t,v]/(v(t-1))). Other fibres are reduced A1 with rank-two Poisson extensions. This is NOT a classification of every special Log Poisson structure.

Seven exact SymPy producers and seven matched JSON certificates are tracked; run:

```bash
bash scripts/verify_github_source.sh incidence
```

No tracked PDF, no GitHub Actions, and no routine Library use.

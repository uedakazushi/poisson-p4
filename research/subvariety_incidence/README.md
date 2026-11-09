# Poisson subvariety incidence — current research module

2026-11-10 / v2. Read [current research policy](RESEARCH_POLICY.md), [notes](NOTES.md), and [the paper](../../papers/poisson_subvariety_incidence_ja.tex).

- Paper: 13-page Japanese LuaLaTeX, containing the smooth-quadric modular multiplicity criterion and generic Rat(1,2) ambient rank-two rigidity.
- Three exact SymPy scripts: [code/subvariety_incidence](../../code/subvariety_incidence).
- Three exact JSON certificates: [certificates](certificates).
- Reproduction manual: [REPRODUCE_GITHUB.md](../../REPRODUCE_GITHUB.md).

```bash
bash scripts/verify_github_source.sh incidence
```

The next problem is the Log(1,1,1) and Aff extension fibres on smooth Q3, including strict stabilizers, actual Jacobi solutions, rank 2/4, Q+C3 vs 2Q+H and the ambient 147-component relationships. Full Q3 classification and external novelty confirmation are not yet claimed. The generic del Pezzo extension rigidity and the earlier valuable results are retained.

GitHub main alone is canonical. No GitHub Actions. Library only after an actual failed GitHub push.

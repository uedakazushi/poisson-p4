# Poisson subvariety incidence — current GitHub module

2026-11-10. This module studies smooth hypersurfaces and (2,2) del Pezzo surfaces with Poisson structures, their ambient extensions, and the geometry of rank-2/rank-4 components. It does not replace the 147-component classification.

- Research policy: [RESEARCH_POLICY.md](RESEARCH_POLICY.md)
- Current notes: [NOTES.md](NOTES.md)
- Japanese LuaLaTeX paper: [poisson_subvariety_incidence_ja.tex](../../papers/poisson_subvariety_incidence_ja.tex)
- Exact producer scripts: [code/subvariety_incidence](../../code/subvariety_incidence)
- Exact, already checked v1 certificates: [certificates](certificates)
- Canonical reproduction instructions: [REPRODUCE_GITHUB.md](../../REPRODUCE_GITHUB.md)

From a complete GitHub checkout with pinned SymPy 1.14 and LuaLaTeX/Poppler:

```bash
bash scripts/verify_github_source.sh incidence
```

This regenerates both certificates and the 11-page PDF, compares JSON values, and checks the generated PDF's text SHA256. A PDF need not be committed as a binary to be fully reproducible from the canonical TeX. The v1 generic del Pezzo rigidity theorem holds for a nonempty Zariski-open set of smooth pairs, **not** for every pair. Novelty/priority are unconfirmed.

**GitHub main is the source of truth. No GitHub Actions. Library is used only if a real push fails.**

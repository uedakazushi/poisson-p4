# Poisson subvariety incidence — current research (v4 / 2026-11-10)

GitHub main alone is authoritative. Read [research policy](RESEARCH_POLICY.md), [notes](NOTES.md), the Japanese LuaLaTeX [main paper](../../papers/poisson_subvariety_incidence_ja.tex) and the included [Log quadric proof](../../papers/quadric_log_extension.tex).

Main results: smooth Q3 generic Poisson extensions for Rat(1,2), Log(1,1,1), Aff are respectively a point, A1 and two reduced points. A special Log boundary has a nodal two-line extension fibre with a rank-four branch. A nondegenerate one-parameter family q_c=Σx_i²+2c*x0*x2 has an orthogonality wall at c=0; the minor 13824c/(c²−1)^4 proves a strict stabilizer jump and the accompanying fibre jump. Full classification of all special Log and Aff fibres remains open. Generic (2,2) del Pezzo rigidity, cubic rigidity, and earlier classification work are retained.

Five independent exact SymPy producer scripts in [code/subvariety_incidence](../../code/subvariety_incidence) recreate the five JSON certificates in [certificates](certificates).

```bash
bash scripts/verify_github_source.sh incidence
```

The command uses LuaLaTeX twice, checks source-generated certificates and validates the resulting PDF. No GitHub Actions; Library only if an actual GitHub push fails. Do not claim that this paper independently re-verifies all 147 components.

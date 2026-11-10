# Poisson P4 — current Codex research instructions

2026-11-11 / v18. Read RESEARCH_POLICY.md first, then research/subvariety_incidence/RESEARCH_POLICY.md, NOTES.md, papers/poisson_subvariety_incidence_ja.tex and its shared papers/del_pezzo_pfaffian_incidence.tex.

## Highest-level objective

Seek a **geometric structural classification** of irreducible Poisson components on complex P4 through Poisson-invariant hypersurfaces and complete intersections, Pfaffian degeneracy quintics, modular vector fields, conormal isotropy Lie algebras, and incidence maps. Explain the existing 147 reduced components and their boundaries; do not blindly increase the number of cases or assume 147 completeness in an independent lemma. Novelty and external priority are not confirmed. Existing long 147-component proof and exact certificates must remain intact.

## Immediate research priority: smooth del Pezzo (2,2) special extensions

For S=V(q0,q1) smooth and 0!=sigma_h in H0(-K_S), all ambient bivector extensions are J(q0,q1,h)/4+q0*A+q1*B with 20 constant antisymmetric parameters. General pairs have a unique reduced rank-two Jacobian extension. Exceptional rank-four extensions for special smooth pairs are STILL OPEN. An arbitrary rank-four Poisson preserving S has Pfaffian D in I_S I_C, where C=Z(sigma_h) is an anticanonical (2,2,1) curve. Necessary degree-five quintic spaces have dimensions h0(I_S(5))=65 and h0(I_S I_C(5))=41.

New v18 proven structural input: the conormal isotropy bracket kappa:wedge²N*_S -> N*_S is a section of O_S(2)^2, and the first normal jet of D is sigma_h*kappa. Hence D is in I_S²(5) iff kappa=0, with 15-dimensional necessary quintic space, compared to 41 for arbitrary kappa. In the nonabelian case kappa has no zeros on the symplectic open S\C; for smooth elliptic quartic C its zero scheme is all C or is finite of length 16. These facts are necessary criteria, NOT a construction or existence proof of exceptional rank-four extensions, nor 26 additional equations beyond the 20-variable Poisson Jacobi scheme. Exact geometric proof in papers/del_pezzo_pfaffian_incidence.tex. Local first-jet and Hilbert-number check: code/subvariety_incidence/check_dp_conormal_lie_pfaffian.py. General first-order conormal theory is known from Fernandes–Mărcuț (JEMS 2026); distinguish its results from the concrete projective Pfaffian formula.

Next analyze the fixed smooth elliptic quartic C and its 10-dimensional family of smooth del Pezzo lifts; compare exact zero-schemes and isotropy brackets of actual rank-four component sources to the 15/41-dimensional Pfaffian ideals. Solve the full 20-parameter Jacobi scheme at special pairs to prove existence or universal nonexistence. Vanishing of a minor is never itself an exceptional rank-four example. Follow with special quadric Rat/Log/Aff incidence and cubic/hyperplane structural unification.

## Existing theorems and independent secondary program

Preserve smooth cubic rank-two unique extension, no smooth degree>=4 hypersurface Poisson-invariant for nonzero ambient Poisson, general del Pezzo rigidity, Q3 generic one/A1/two extensions, resonant isotropic quadric boundaries, hyperplane transverse symmetry criterion. Preserve the separate isotropic Log tangent-cone paper (dimension30 tangent, dimension28 local, nonreduced tangent cone, 23/24 component intersection). Its higher cubic initial equations and ambient local ring reducedness are unresolved but **not the current top research priority**.

## Reproduction and persistence

GitHub main is the only canonical source. Update with non-force expected-head commits and read back the remote SHA. Track full current mathematical proofs in Japanese LuaLaTeX, exact code, small certificates and input closure. Do not track generated PDF, cache or logs; no GitHub Actions. Library only after an actual failed GitHub push. No history archives or redundant ZIPs, and do not delete Library assets without instruction. Run bash scripts/verify_github_source.sh incidence to reproduce current 13 exact producers and paper, and record unrun verification honestly. No external publication without explicit authorization.

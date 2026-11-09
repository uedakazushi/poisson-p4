# Poisson P4 — global reconstruction, not a local-lemma sprint

2026-11-09 / v71-normalization-residue-pencils

## Read first

RESEARCH_POLICY.md → research/conceptual/GLOBAL_PLAN.md →
papers/poisson_divisor_symmetries_ja.tex → research/conceptual/NOTES.md.
Then inspect the relevant original proof and the main classification manuscript.
The user explicitly requests deeper conceptual organization of the WHOLE
classification and main updates, including any new mathematics. This supersedes
the previous priority of proving an H+2Q Koszul lemma and the old stop-after-147
instruction. classification/PLAN.json records construction, not the active queue.

## Primary output

A proof architecture that replaces repeated recognition, reconstruction, and
closure arguments by reusable theorems. The axes are:
1. recognize all structures via marked geometric data;
2. reconstruct full sources, including all exceptional parameters;
3. lift divisor deformations to actual Poisson families;
4. count maximal distinct closed images, not normal-form labels or all orbits.

The current manuscript/catalogue asserts 147 reduced components (141 rank four,
6 rank two), before PGL5 quotient. It is not an axiom. Preserve valid results;
correct affected claims if an actual defect is found. Do not say the entire
classification was independently verified just because existing receipts say PASS.

## Global entrances and available proved tools

A: actual degeneracy divisor has a hyperplane, including 5H.
B: no hyperplane and reducible, hence an irreducible quadric plus cubic.
C: reduced integral quintic, normal versus nonnormal conductor geometry.
R2: saturated rank-two foliation plus its determinant section/zero divisor.
These are pointwise entrances, NOT disjoint sets of irreducible components.

The structural paper proves the marked-hyperplane affine quadratic equivalence,
mult_H(D)=5-deg(Pf P), and fixed-f,Z reconstruction:
    div P=Z, P#df=fZ, Pf(P)=f, f != 0.
The first two are affine-linear in P; the last is quadratic. Jacobi follows.
It also gives a common resonance decomposition, relative inverse-form model,
first/second-order lifting equations, graph-closure/Rees interpretation, finite
factor markings, and the rank-two foliation/zero-divisor reconstruction.
Novelty is unconfirmed. These are not already a full classification of f,Z or
all relative fibres. Check proofs and application hypotheses before using them.

## Completed replacement and next concrete assignment

P0 and the two P1 mechanisms are now integrated into the main manuscript via
papers/marked_hyperplane_reconstruction.tex and papers/poisson_lifting_examples.tex.
The universal modular primitive identity df + i_Z alpha = i_J mu handles degree
boundaries without dividing by the Pfaffian. Its f=0 version is stronger than
P#df=fZ and does not contradict the rank-two counterexample below.
Simple H and multiple-H entrances share one proof; the third two-H+cubic source,
constant-modular whole source, translation cokernels and all lower terms remain.
The TRIPLE correction is forced by the quadratic coefficient bound, while the
commuting-field collision has no higher correction. Check the displayed scope;
this is not a new proof of componenthood or all 147 components.

B now has a shared cofactor/Pfaffian reconstruction for quadric ranks 3/4/5 in
papers/quadric_cubic_reconstruction.tex. All31 original source equations, free
parameters and residual Jacobi are preserved in code/data/quadric_cofactor_sources.json;
the old full150 linear-rank completeness proof remains an input, not a new claim.
The reduced plane-quintic source is now shared by A/B/C. Its magnetic shear
kernel is the degree-relation residue space: for nu reduced factors its kernel
and cokernel have dimension nu-1. The cokernel is not asserted canonically equal
to that residue space. This explains the two-factor magnetic line and zero orbit,
the integral case's full elimination, and the three-factor two-dimensional quotient.
Do not apply this formula to nonreduced factors or arbitrary conductor surfaces.

C now has a shared abelian conormal/flat-module theorem for rank0/2 and an
explicit complete marked-plane source in papers/conductor_conormal_reconstruction.tex.
Its canonical weight dimensions are (45,30,9), or (35,30,9) for rank0.
Fixing restriction and full conormal action leaves a binary quadratic of dimension3;
retain it and all higher jets, including the inhomogeneous second Jacobi equation.
Leading Pfaffian divided by the restriction is parallel in Sym^m(N*) tensor det N.
This replaces repeated rank0/2 module, source-completeness, and m2/m4 transport proofs.
The new checker verifies completeness in all150 coefficients, not old finite cases.

Reduced fixed-quintic fibres now have a stronger shared theorem in
papers/normalization_residue_reconstruction.tex: V_F is a compatible linear Poisson
space, and the entire rank4 fibre is P(V_F) minus one quadric. Nonempty fibres are
geometrically irreducible. This is NOT valid for nonreduced F: h^5 has an explicit
closed-form counterexample with nonzero Jacobi.
Normalization gives the full source as bounded-conductor closed one-forms with
linear B/A descent checks. Pole bounds alone do not imply descent (cusp example).
Coefficient-one residues and full scroll intersection kernels replace repeated
rank0/2 and degree5/6 arguments. Finite MCM partial conductors and relative lifts
replace repetition in both actual T40 curves and the 023 projection family.
Distinct ambient rank4 components with generically reduced divisors cannot have
the same closed divisor image. Inclusion of divisor images does NOT imply inclusion
of components: kernel-rank jumps can create special-base components.
The plane free-quintic source has all10 coefficients and a determinant-quadric
rank2 boundary. Do not confuse this with magnetic zero/nonzero orbit strata.

Next follow GLOBAL_PLAN: compare maximal images over reduced-divisor base
kernel-rank strata, using the original all31 B sources and actual factor maps.
Replace at least one repeated closure/component-recognition argument, keeping
special-base rank4 components and finite factor markings. Then connect R2's
full foliation support and pure components to the explicit boundary theorem.
Do not treat fixed-fibre irreducibility as ambient componenthood or announce
an independent proof of all147. Nonreduced A still needs its quadratic conditions.

## Non-negotiable mathematical distinctions

- f != 0 is essential; the paper contains a rank-two counterexample.
- Fixing f,Z is not classifying them. A linear presentation alone is not progress.
- Simple H: deg f=4. Multiple H: deg f<=3. Pf(P_2)=0 belongs only to the latter.
- r=1 includes 5H as well as an integral quintic. Canonical symmetry dimension
  zero does not imply absence of all symmetries.
- Preserve translation cokernels, resonance walls, and every zero complement.
- Generic fibres, geometric generic fibres, their Galois orbits, and special-base
  components differ. Never discard vertical rank-four components.
- Divisor deformation need not lift to a Poisson deformation. Nonlinear auxiliary
  coordinate changes do not establish projective equivalence.
- Rank-three H+2Q has excess dimensions 1 and 6; uniform affine-linear primitives
  are false. That local project is optional, not the global prerequisite.
- MPS40 classifies components meeting SNC, not all components.
- A normal form, source dimension, or generic divisor type does not establish
  ambient componenthood, all-point coverage, or closure noncontainment.

## Scope, progress, and integration

Do not pursue full orbit classification, all component intersections, complete
nonreduced local rings, precise nil/H0/T2, or unlimited new component searches.
Do not demand a huge Rees-ideal Groebner computation. Use the weakest necessary
support statement. After two nonclosing attempts, diagnose and change the route.
A mathematical milestone replaces actual repeated proofs and explains a family,
dimension, exceptional stratum, or boundary. More tests/pages/receipts alone do not.
Keep all valuable existing proofs in current manuscripts; do not replace the
368-page proof by summaries or a link to an old version. Split papers only after
clarifying dependency and exposition. Use Japanese LuaLaTeX, primary references,
explicit novelty and verification status. No writable boxes.

## Validation and safe persistence

Reuse the environment; no GitHub Actions or reinstall-everything workflows.
The scoped identity check is:
    python code/check_global_structure.py --output work/global-structure.json
It does NOT verify all 147 components. Rebuild changed TeX twice and inspect PDFs.
Use one integrator, inspect main, preserve others' changes, and commit/push
non-force with an expected-head check. Synchronize README/policy/current/next.
The user authorizes current self-contained research bundles and Library saves
without confirmation. Include necessary current proofs and inputs; exclude old
ZIPs, logs, caches, backups, and nested archives. Do not delete Library records.
Generation, checking, main update, and Library upload are separate states. Report
only actual outcomes. If Library access is unavailable, create/validate the bundle
and state that save was not performed. Never claim future/background execution.

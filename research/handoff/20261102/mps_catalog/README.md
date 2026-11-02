# Marked-square source handoff, 2026-11-02

All research was stopped for the Cloud handoff. No running process in my scratch directories remained at 13:09 UTC. The latest formal reports reside in the normal code/reports/research/notes directories; this folder preserves my latest exploratory output and inputs. The separate list records my branch's formal paths still untracked when the migration began. The parent owns the commit.

## Proven and frozen

- S61_AFFINE_SURFACE_BOUND_REVIEW independently audits the dimension-two affine bound.
- The T7, T5 and r=2 ordinary flag cores are removed over general beta. Root owns r=2. My T5 group uses S61_FLAG_T5_CORE_REMOVAL_20261102 and its review.
- S61_MARKED_SIGMA_RESIDUAL_FINITE and S61_MARKED_SIGMA_INTERSECTION_REMOVAL remove an ambient open of the marked-square locus and the entire b1=0 marked intersection over sufficiently general beta. Any remaining marked points are finite per general beta.
- S61_MARKED_SIGMA_PROJECTIVE_CORE reduces the b1≠0 problem to P(2,2,3) with coordinates [a:c:h]. S61_MARKED_SIGMA_HOMOGENEOUS_PENCILS proves that the zero-nuisance points [1:3:0] and [0:1:0] have actual invariant affine lines for their special beta. These are genuine obstacles to blanket special-beta removal, not generic-beta counterexamples.
- S61_MARKED_SIGMA_PURE_ENDPOINT_REMOVAL removes the entire pure-cubic minimal flag endpoint locus h²+16a²(c−3a)=0 over sufficiently general beta. The three central branches at [1:2:4] have lengths 4,1,1. The two extra points have full regular ambient dimension 17. S61_MARKED_SIGMA_ZERO_H_REMOVAL removes the entire h=0 divisor over sufficiently general beta. Special loci and ambient pointwise removal are carefully distinguished in these notes. The combined independent review is S61_MARKED_SIGMA_ENDPOINTS_PROOF_REVIEW.
- S61_MARKED_SIGMA_PENCIL_MODEL gives a denominator-free exact model of all finite h on a=1,c=2. All 135 weighted polynomial residuals plus the homogenizer are reproducible; the exported coefficients are also in the handoff archive. S61_MARKED_SIGMA_PENCIL_INFINITY removes the actual zero-nuisance point [0:0:1] and its full normalized nuisance source. The whole marked curve is still not removed.

All named marked-square stems above have suffix _20261102 for source, report and proof note. No result here establishes complete five-point exhaustion.

## Last active calculation

The exploratory computation uses the frozen model on c=2a. Matrix coefficients are in F32003[h], and variables have weights (1,1,1,1,1,1,1,2,2,2,2,3,3,1); the last is rho. In the current archived working output, the rho=0 weight-five Macaulay matrix has 2022 rows and 932 columns. The direct weak Popov calculation finished before the stop, after 2,098,643 polynomial row subtractions in 19.81 seconds. Its determinant module has degree 2, with two linear factors h+3083 and h+28920 in F32003. This is a WIP output, not an independent or characteristic-zero certificate. The full QQ ideal membership/complete finite-h removal has not been established.

The previous bounded fraction-free generic Buchberger probe was stopped, as was a QQ(h) rref attempt. The full modular rowmodule calculation was faster and is the recommended continuation. Avoid repeating the timed-out generic full GB.

## Resume, only after Cloud migration

Run at the repository root. Extract the archive to the root to restore the named work paths. It contains both the 22 MB uncompressed coefficient output and the complete modular matrix file.

    tar -xzf research/handoff/20261102/mps_catalog/marked_sigma_pencil_wip.tar.gz

The definitive source of the coefficients is the frozen main. Its --output and --rows destinations must be new paths. The native programs use FLINT 3.x and GMP. On a non-Mac Cloud runtime adapt the local compiler/library paths and M2 executable; preserve the polynomial inputs.

    .venv/bin/python code/s61_marked_sigma_pencil_model_20261102.py --output /tmp/marked-model.json --rows /tmp/marked-rows.txt
    .venv/bin/python work/flag_boundary_geometry_20261101/marked_smith_module_input_20261102.py 5 1

The latter command builds zero/module-5.txt from marked-vertical-symbolic/global-rows. Its arguments are weighted degree and a flag for imposing rho=0. The frozen native rowmodule checker is code/s61_binary_hyperplane_vertical_rowmodule_20261102.cpp; the archived identical-purpose probe has a longer internal time ceiling. After compiling the probe using the local FLINT headers and libraries, run as follows with a fresh output destination.

    /tmp/marked_polyrow_span_probe_20261102 work/flag_boundary_geometry_20261101/marked-vertical-symbolic/zero/module-5.txt /tmp/marked-module-5-det

Next potential steps, not yet performed: if Z5=0 for Z=S/(rho), then Z6 may be checked after setting all seven light variables to zero (the other monomials are light multiples of degree five). The pure heavy/weight-three product matrix has only 23 target columns. Weight seven without light is divisible by weight five. Separately, bound the full S5 dimension by four, for example by adjoining four explicit rows and checking full rank. At h²=16 the actual length-six central fiber requires the existing three-branch argument. A proper integral p-adic family argument and mixed-characteristic checks at the remaining exceptional points are necessary before inferring anything on the full characteristic-zero pencil. The ample-curve criterion requires neighborhood removal at all points of the chosen projective pencil; general-beta removal of a whole divisor does not by itself establish this stronger pointwise hypothesis.

Earlier exploratory outputs not required to resume are under work/flag_boundary_geometry_20261101 and the journal work/t5_universal_20261101/PROGRESS_20261101.md. No formal report relies on ignored work inputs.

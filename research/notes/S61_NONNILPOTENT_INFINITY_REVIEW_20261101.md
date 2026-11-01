# Independent review of nonnilpotent infinity exclusion

Date: 2026-11-01. Review of the corrected `S61_NONNILPOTENT_INFINITY_EXCLUSION_20261101` result.

The corrected theorem is accepted: the image in W16 of all nonnilpotent modular points satisfying the 134 leading equations has closure of dimension at most fifteen. Consequently it is absent over a nonempty base open; adding G6 preserves that exclusion. This does not exclude such points over every smooth-curve first jet, identify the saturated projective closure, or bound the affine fiber.

## Coverage correction and historical dependencies

The first audit found a genuine missing chart. The displayed G2 chart sets A=1 and therefore does not cover A=0, B nonzero. The one-zero-trace classification does not repair this omission by itself: its semisimple SL2 matrix T was referred back to the same log charts. In particular the trace-xy plane with d=1/2 has Qxy=Qzw=xy/2, whereas the d=0 specialization of the displayed charts has Qxy=xy, Qzw=0. These two planes cannot be conflated.

The corrected code includes the entire additional family N(0,1/2,0,1/2,e), with all four rank drops. This removes the gap. The review report hashes only the corrected files.

I read the archived v03 regular-semisimple proof, v04 semisimple-boundary proof, and v05 mixed-Jordan proof and compared their formulas with the representatives. They supply point coverage, including resonances and zero parameters, rather than a claim about dense open subsets of components. Specifically:

- The regular semisimple resonances are represented by T1,...,T4 and log points, with the Darboux exchanges and irrelevant nonzero scalar allowed.
- In the one-zero case, the two nonzero nilpotent SL2 orbits are both retained. The semisimple cases are now covered by all three log charts.
- In the equal-eigenvalue case, all ten nonzero exceptional curves in the thirteen-branch table are represented, including both K signs. The three baseline curves and U=V=0 orbit are log cases; the dual semisimple chart is covered by the stated symplectic duality. The archived unit-ideal computation excludes the remaining nonzero nilpotent U or V.
- The two nonnilpotent mixed Jordan types are accounted for: one has empty fixed fiber, the other has exactly the three complete affine lines used in the new calculation.

This review reads and uses those archived exhaustive calculations; it does not independently rerun their original primary decomposition computations.

## Independent coefficient and exceptional-set verification

The review script does not import the main exclusion script or its leading-matrix helper. It reconstructs each displayed tensor in S61, substitutes beta and h times its quadratic coordinates into every stored equation, and extracts the h-squared coefficient. It checks all entries of the reported block matrices, all sixteen fixed-representative matrices, and all three mixed-line matrices.

For every two-parameter rank-one block, an independent rational Groebner calculation determines the complete zero scheme of its entries. Its quotient length equals the number of distinct listed rational zeros, and each listed point satisfies the ideal. This proves there are no missing complex zeros or hidden zero-dimensional lengths. The unions contain exactly 27 and 23 points; full matrix ranks at those points reproduce the report. In particular this does not rely on testing a finite grid or merely reproducing the factor-line enumeration.

For the corrected one-parameter boundary, the independent entry ideals give precisely e=-1,-1/2,0,1/2. The generic rank is twelve and the rank at all four points is eleven. The nonzero fixed and mixed minors are recomputed exactly; the latter remain constants in the complete line parameter.

The independently verified calculations use the same archived S61 reconstruction and equations, which are an explicit shared dependency.

## Geometric dimension argument

The diagonal Cartan fixes each log tensor, so each Sp4 orbit has dimension at most eight. The two-dimensional open log charts have incidence bound 2+8+(16-12)=14, and their finitely many exceptional points have bound 8+(16-10)=14. The additional one-dimensional chart has bound thirteen, both generically and at its four exceptions.

The sixteen fixed representatives have beta rank at least eleven, so even the unquotiented Sp4 parametrization has dimension at most 10+5=15. For a mixed line the nonzero modular vector field is a strict infinitesimal symplectic stabilizer; its orbit dimension is at most nine, giving 1+9+5=15 uniformly. Equivalently, the stabilizer acts freely on the group coordinate of the incidence presentation and preserves its beta kernel. This justifies subtracting its dimension even when the stabilizer varies.

Scalar multiplication of Q contributes no parameter: Q is projective and its homogeneous leading equations have the same beta kernel after nonzero scaling. The finite union of constructible images therefore has closure of dimension at most fifteen. Properness of these parametrizations is unnecessary for that dimension assertion.

## Reproduction

`code/s61_nonnilpotent_infinity_review_20261101.py --output <new-path>` produces `reports/S61_NONNILPOTENT_INFINITY_REVIEW_20261101.json`. The report records the independent Groebner bases and zero-scheme lengths, exact exceptional ranks and minors, and input hashes. The review introduces no additional full-fiber or global-exhaustion claim.

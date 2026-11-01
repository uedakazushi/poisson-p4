# Independent review of the minimal-infinity completed local rings

Date: 2026-11-01. Accepted, with the source-to-canonical M1 normalization now made explicit in the main proof note.

The independent checker does not import either main certificate generator or its polynomial-specialization helper. It regroups the original encoded S61 equations by the forty-five quadratic coordinates, reconstructs their projective homogenized first and second derivatives, and checks the saved data against those derivatives.

At each of the four source points it verifies the full Jacobian, tangent basis, selected cokernel rows, quadratic obstruction coefficients, full obstruction rank and one-dimensional common Hessian radical. The resulting transverse length upper bounds are respectively 4,3,3,2 for M0-134, M0-135, M1-134, M1-135.

Over QQ(c,b0,...,b8), it independently rebuilds the canonical tensors and normalized coordinates, substitutes the complete beta kernel, and reconstructs the leading Jacobian and homogenizer column. It verifies the first-block kernel and pivot determinant, reconstructs the residual matrices and columns, and checks each saved rational residual kernel by exact multiplication. The free-coordinate submatrix is the identity, so these kernels have the asserted independent dimensions. A rational specialization also confirms the corresponding rank lower bounds. Thus all universal rank upper bounds used by the proof are identities, not sample-based extrapolations.

The finite-module argument is valid. The smooth curve gives a retraction S to R=C[[u]] and a section after lifting u. A finite transverse quotient makes S finite over R with at most r+1 generators. At the generic point of this same curve, the universal Jacobian upper bound gives transverse embedding dimension at least r, hence local length at least r+1. This forces S to be free of rank r+1 over R. Its augmentation ideal I has rank r; minimal generic length gives I_K²=0, and torsion-freeness gives I²=0. No additional generic factor or branch can remain within this completed local ring.

The homogenizer belongs to I. Its nonzero linear part in the first three cases proves t is nonzero there, and t²=0 follows from I²=0. In M1-135, the universal extra rank of the full Jacobian makes t zero in the generic cotangent space; since I_K²=0, t=0 generically and hence everywhere by torsion-freeness. The completed-local identities descend through faithful flatness. All four points are therefore absent from the t-saturated closure.

One normalization detail was clarified during review: the universal M1 generator uses Qxy=Qzw=x² and coefficients -(c+1)xz, (c-1)xw. The original source family becomes this tensor after projective scaling by two, with c_can=2c_source+1 and a_can=2a_source,b_can=2b_source. The main note now states this before its shear reduction. It changes no computed rank or conclusion.

Frozen independent files:

- `code/s61_minimal_formal_rings_review_20261101.py`
- `reports/S61_MINIMAL_FORMAL_RINGS_REVIEW_20261101.json`

Reproduction:

```
.venv/bin/python code/s61_minimal_formal_rings_review_20261101.py --output work/minimal_formal_review_check.json
```

The exact run passed in 141 seconds. The review concerns the stated points and the open part of the two curves justified in the main note. It does not extend the completed-ring description to resonance fibers, all infinity points, or the full affine S61 fiber.

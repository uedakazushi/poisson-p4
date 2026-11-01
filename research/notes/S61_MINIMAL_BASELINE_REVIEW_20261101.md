# Independent review of the four baseline thickenings

Date: 2026-11-01. Accepted on the specified points and justified opens. The M0 completed local ring is R⊕R³ with square-zero second summand and nonzero t satisfying t²=0; the M1 ring is R⊕R² with square-zero second summand and t=0, where R=C[[u]]. These statements do not cover all endpoints or intersections of the four baseline curves.

## Exact checks and the additional beta direction

The independent script uses the original134 equations and universal G6 polynomial directly. It does not import the main specialization, beta-matrix or local-data routines. At both canonical tensors it reconstructs the rank-six beta equations and verifies that the saved rank-ten matrix is their full kernel. Every leading135 equation vanishes identically on that kernel.

Both leading44-column Jacobians and full t columns are reconstructed over QQ(b0,...,b9). The saved rational kernels are checked against all135 rows. Lifting the residual kernels back through the first-block elimination gives full independent tangent kernels of dimension four for M0 and three for M1. Their t-zero subkernels have dimension three. A rational specialization attains the corresponding ranks, proving leading rank41 and full ranks41/42, together with their universal rank upper bounds.

The M1 top70 compatibility vector is explicitly nonzero and retained. The review confirms it is zero only in M0. Thus the extra tenth beta direction at resonance is included, rather than borrowing the earlier nine-variable compatibility calculation.

At each stated rational point, all135 full Jacobian entries are independently reconstructed. The reviewer computes a fresh left kernel, the complete tangent quadratic coefficient matrix and its cokernel projection. The obstruction ranks are six for M0 and three for M1. Their common Hessian kernels are one-dimensional and contain exactly the nonzero tangent of the same actual fixed-beta curve. The full curve satisfies every leading equation. The t coordinate on the full tangent is nonzero for M0 and zero for M1.

The parameter-dependent shears eliminating a,b and the finite (z,w) exchanges are verified as full symplectic tensor identities. They transport the same curve to the canonical tensors and identify the paired baseline strata. No parameter denominator is discarded in these affine shears.

## Why the transverse bound is exact

Let the full tangent dimension be n=4 or3. The quadratic forms have common radical equal to the curve tangent and span dimension n(n-1)/2. They therefore span all quadratic forms on the (n-1)-dimensional normal quotient. Cutting by a linear parameter u whose restriction to the curve has nonzero derivative leaves all normal quadratic monomials in the initial ideal. The transverse completed ring has length at most n. Its embedding dimension is n-1, so length is also at least n; the length is exactly n.

This proves finiteness over R=C[[u]] and an R-module generating bound n by complete Nakayama. It does not yet prove that the higher relations are square zero; that requires the generic lower bound along the same curve.

## Generic lower bounds and finite freedom

The canonical full-kernel rank identities are polynomial upper bounds, valid at every beta after specialization. The finite symplectic normalization carries an open of the actual fixed-beta curve to this canonical family. Full homogenized equivariance is an explicit prior dependency. Thus:

- For M0, the full relative tangent dimension at the generic curve point is at least45-41-1=3. Its finite local factor has length at least four.
- For M1, the leading scheme S/(t) has relative tangent dimension at least44-41-1=2, so its generic curve factor has length at least three.

The upper generating bound and these lower bounds force rank_R S=n. The kernel of the surjection R^n→S is rank zero inside a free module over the domain R, hence zero. Therefore S is free of rank n. The curve factor uses the entire generic rank, so it is the unique generic support.

The smooth embedded curve gives an R-algebra retraction S→R. Its kernel I is free of rank n-1. In the generic Artin factor, length n and embedding dimension n-1 imply I²=0. Since S is R-torsion-free, that relation extends to the whole ring. Choosing an R-basis of I identifies S with R⊕R^(n-1) with square-zero multiplication on the second summand. Freeness alone would not have implied this multiplication law; the generic length and tangent bound are essential.

For M0, t vanishes on the reduced curve, so t lies in I and t²=0. Its nonzero tangent class proves t itself is nonzero. For M1, the quotient S/(t) already has full generic rank three; hence (t) has rank zero inside the torsion-free module S and must be zero. Faithful flatness of completion descends the corresponding nilpotence/vanishing to the local algebraic rings, so both specified points disappear after t-saturation.

The finite symplectic exchanges transport these conclusions to nonempty opens of the paired M0 and M1 baseline curves. The theorem correctly does not assert that a generic-point saturation statement removes every point of a seventeen-dimensional family over W16.

## Reproduction and scope

`code/s61_minimal_baseline_review_20261101.py --output <new-path>` produces `reports/S61_MINIMAL_BASELINE_REVIEW_20261101.json`. It records the independently verified full beta kernels, universal rank bounds, sample Hessian dimensions, actual curve radicals and finite changes. The historical baseline classification and global homogenized-equivariance theorem remain explicit dependencies. No affine-finiteness or five-point-exhaustion conclusion is drawn.

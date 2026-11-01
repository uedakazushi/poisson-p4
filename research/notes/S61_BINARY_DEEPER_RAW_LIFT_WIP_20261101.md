# Deeper binary: a shorter lift and two cubic vectors (WIP)

This records a continuation after v50. It proves no new core-curve removal,
ambient nilpotence, or generic-fiber upper bound. The frozen archive
`research/archive/binary_deeper_pid_20261101/` is unchanged.

## Computations completed

On the nuisance-zero core, the nonzero weight coordinates of the primitive jet
are

\[
 b_1=d,\quad b_9=-3(2e-1)(15e+8)/500,\quad
 b_{10}=(3e-1)/5,\quad b_{13}=3e/5,\quad b_{14}=-2.
\]

Their torus weights satisfy \(a-b=-1\). Combining the corresponding diagonal
symplectic action with jet scaling fixes the core. Direct inspection of the
saved tangent matrix, after setting the five nuisance parameters to zero,
gives the nine normal weights
\((2,2,2,2,2,4,2,4,2)\). The two weight-four variables are exactly \(v_5,v_7\),
whose products with every normal variable already belong to the quadratic
ideal. Consequently the remaining seven variables all have weight two: this
action does not divide the 210-dimensional fourth-degree space into smaller
weight blocks.

It does divide the canonical quadratic basis into seven weight-four rows,
fourteen weight-six monomials \(v_5v_i,v_7v_i\), and three weight-eight
monomials. Here \(i\in\{0,1,2,3,4,6,8\}\). A new exact contraction over QQ
reconstructs the canonical cubic tails from the **original** quadratic
implicit lift. After projecting \(v_5=v_7=0\), all ten tails in the weight-four
and weight-eight groups are zero. The fourteen remaining tails form two
seven-component cubic vectors. Their rational coefficients are saved in the
scratch output `/tmp/poisson_deeper_raw_canonical_cubic_20261101.txt`; that file
stores an extra factor \(v_0\), which the extraction routine removes exactly.
This is a smaller description to investigate, not yet a certificate of full
fourth-degree rank on the entire core.

The earlier construction first reduced the quadratic lift modulo the
quadratic initial equations. That introduced large denominators. The new
construction retains its 29 nonzero coordinate rows, including the
homogenizing coordinate. Its Hessian contraction therefore includes the
mixed terms involving that coordinate. The established change-of-lift
argument shows that these cubic tails change by \(I_2S_1\), and their exact
linear-syzygy quartics change by \(I_2S_2\). It is essential to retain the mixed
terms; simply reusing the former 17-coordinate operator is incorrect.

All 86 former syzygies have also been contracted with the shorter lift. The
first nonzero numerator decreased from 99,940 to 56,092 terms; the complete
scratch output decreased from about 909 MB to 481 MB. This is useful but does
not by itself settle the remaining exceptional curves.

## Current bounded attempt and limitations

The seven projected quadrics have 49 linear multiples of rank 48 at the
stored finite-field witness. Their unique linear syzygy has coefficient-matrix
rank two there. An exact QQ(e,d) elimination of their cubic multiples is
running; the intended next operation is to reduce the fourteen cubic tails
modulo those multiples before forming fourth-degree relations. Any new pivot
factor must be recorded and covered, rather than silently inverted.

Other completed probes did not close the problem. Allowing individual
Gaussian pivots supported only on the seven treated horizontal factors gave
no further pivot after the original 55. The canonical quadratic presentation
gives 104 such pivots in degree four, but its remaining block is unresolved.
At \(e=2\) modulo 32003, every row and column content of the earlier residual
matrix is constant, so a simple common-content division cannot remove the
large annihilator powers. A degree-36 coefficient Ansatz contains all target
columns for the known-factor annihilator at this finite-field slice; degree
28 does not. No QQ(e)[d] left-inverse identity has been certified.

Reconstruction entry points in ignored `work/` are
`binary_deeper_raw_lift_prepare_20261101.py`,
`binary_deeper_raw_lift_quartics_20261101.cpp`, and
`binary_deeper_reduce_canonical_cubic_20261101.cpp`. Their mathematical inputs
are the frozen triangular-chart, uniform-quadratic and cubic certificates,
plus the v50 archived native reconstruction inputs. The new aggregate route
has not received independent review. The remaining unknown fourth-minor core
curves, and all subsequent global exhaustion steps, remain open.

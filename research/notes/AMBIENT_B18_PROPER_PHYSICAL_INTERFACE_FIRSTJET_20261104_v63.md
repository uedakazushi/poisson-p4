# Actual proper-source / physical coefficient interface: unit first derivative

2026-11-04. The exact full source coefficient map is an algebraic unit-implicit
operator. This packet pins its original physical inputs, source flags, root
signs, group gauge and target-relation row identity. It computes their exact
constant Jacobians by first differentiating the already frozen two48 systems;
no new scalar jet, Gröbner basis or source worker is run.

The actual nine-coordinate minor is **5/48**, and the actual twelve-relation
row minor is **687194767360000**. Both are nonzero rational numbers. The one
successful arithmetic run took less than one second and about 72 MiB RSS.
The serialized matrix-typed operator lists every input reference, all fifteen
gauge equations, all flag/frame operations, and the finite first-five pivot
telescopes. It has no unspecified coefficient or nonunit quotient.

## 1. Actual normalizer input and root signs

Use exactly ZERO_NORMALIZER_TENSOR_MODEL's same-R0 coefficient nodes and two
48 unit systems, with their recorded unknown/equation order. The own source
parameters are p_i=tau_i+s_i/2. Each outer gamma and core k parameter is its
frozen pair a+b*p_i, not a newly chosen trace basis. The underlying five by
five source frame is H_i=G_gamma_i*core_group_i(p_i), using the core groups
and parameter conventions in QUADRATIC_IMAGE_REDUCTION. The separate scalar
is 1+gamma0. This pins both source-root signs before any coordinate adapter.

The physical coefficient algebra C has two monics. Their unit minor in
(x27,x61) is -20/9, leaving actual regular coordinates

    (x108,x109,x112,x113,x117,x118,x119,tauA,tauB).

Both eliminated coordinates have zero first derivative in those nine. The
code independently differentiates each frozen48 system and checks every
linear residual and the original48 determinant. The complete first source
derivatives, including the means s_i/2, are retained in the evidence.
The literal paired-unknown equation order has the negative of the previously
stored even24/odd24 block determinant; the exact row and column permutations
are recorded and both determinant conventions are verified.

## 2. Literal pointed source gauge operator

For any frame H near identity, normalize row_x(H^-1) to have x coefficient
one; call it ell. Subtract row_h[x]*ell from row_h(H^-1), and normalize its h
coefficient to one; call it k. Their five remaining coefficients are the
flag chart (a,b,c,d,e). Every inverse just used has closed-point value one.

Let G(t)=I+sum t_j E_j use the frozen fifteen chosen GL4 positions. Apply it
to the CA source, set its five normalized flags to zero, and set the ten
own fixed-flag fiber readbacks to the frozen center. Its fifteen-equation
Jacobian in t is the frozen matrix of determinant -2/3. This operator lives
in the CA rank-two algebra; replacing each t_j by a+b*tauA gives a pointed
thirty-scalar unit system, determinant (-2/3)^2. It requires no conductor or
root inversion. Its free nine CA fiber readbacks and the normalized flags
of G(t)H_B are exactly the fourteen proper-source functions.

Retain the nine functions (a,d,v2,v3,v4,v5,v6,v7,v8). The evidence computes
their exact Jacobian in the nine physical regular coordinates and records a
nonzero determinant. Hence this is an actual pointed algebraic formal
coordinate map, not a comparison of abstract ring types. Its inverse is the
corresponding nine-unit Hensel operator; it is not claimed rational.

## 3. The whole target-relation row identity

Write Delta=P_A-P_B. In the actual same-R0 cut, only the twelve frozen target
differences F12 remain independent. The common group gauge and relative CB
flag frame give the exact identity

    source_H17 = Ann_CB * Action(F_CB G)(Delta).

The right side is linear in those twelve target differences, with actual
unit-implicit coefficient functions. The CB term is annihilated because its
frame is its own actual flag frame and its nineteen source coefficients
remain in that fixed linear source space.

Now use the five explicit rational pivot roots from PROPER_RATIONAL_ELIMINATION.
Finite telescoping at those roots expresses each of the remaining seventeen
source rows as its residual plus multiples of the five selected source rows.
More precisely, if p=(v0,v1,b,c,e) and p_star are the five explicit roots,
write H_sel=A5*(p-p_star), H_rest=H_rest(p_star)+B12*(p-p_star).
Each entry of A5,B12 is a finite geometric sum of monomials, evaluated at the
actual source readback and those explicit roots. The determinant of A5 has
origin 8/3, so its inverse is a known local unit inverse. Scaling by the
frozen row contents N converts the twelve residuals to dH12. Thus

    R12=diag(N)*(R17_rest-B12*A5^-1*R17_sel),
    dH12(k(C)) = R12(C) F12.

The exact origin matrix is computed as the source first-five Schur complement
of Ann_CB restricted to the twelve physical target columns, with the original
row-unit normalization. Its nonzero determinant is in the report. Therefore
R12 is invertible in the actual local coefficient algebra; its inverse uses a
known center unit. This proves the all-order row identity and coefficient
transport recipe while preserving all five d*u relations after the later
regular-eight adapter.

## 4. Representation boundary

The evidence contains all source inputs, constant matrices, both unit minors
and the fully specified nonlinear matrix-typed unit operator. The source9 map
and exact row12 identity therefore give an actual ambient regular-ring
coordinate isomorphism, not merely an abstract quotient isomorphism. A literal
inverse exists as its pointed nine-unit Hensel system. It does not yet expand
the new group30, coordinate9 and target-row12 arithmetic into the old scalar
coefficient DAG or extract their four R0 components. That remaining
finite serialization/evaluation gate is explicit; no 128-column R0 syzygy
matrix, comparison kernel, ell-torsion or D0 depth is claimed here.

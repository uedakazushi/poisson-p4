# The T37 finite wall is a quartic field plus actual fifth-root scalings

2026-11-07. Hand algebra and actual diagonal coordinate descent. No new matrix-kernel calculation, source point, numerical root or CAS job is used. Existing runner/jobs/checkpoints remain frozen.

## 1. The small quartic algebra

In the actual finite wall presentation, put `v=a5`, `w=v5`, `Z=a3/v²`, and `U=34992w`. The constant term of the saved v-polynomial is nonzero, so v is a unit. The two saved equations become exactly

\[
3U^2+11U+9=0,\qquad
Z^2-(10U+15)Z+\frac{33}{4}U+\frac{99}{2}=0.       \tag{1}
\]

The constants used are `34992=2^4*3^7` and `408146688=2^8*3^13=34992²/3`. Thus the twenty-dimensional algebra is a four-dimensional algebra with `v5=w=U/34992` adjoined; a3 is recovered as Zv².

This quartic algebra is a field K. Indeed `F=Q(U)=Q(sqrt(13))`, with

\[
U=\frac{-11\pm\sqrt{13}}6.
\]

If `J=2Z-(10U+15)`, the second equation gives

\[
J^2=100U^2+267U+27
   =-\frac{13}{3}(23U+63)
   =-\frac{13}{18}(125\pm23\sqrt{13}).           \tag{2}
\]

Both real embeddings of the right side are negative (`125>23sqrt(13)`). It therefore cannot be a square in the real field F. Equation (1) is quadratic irreducible over F, giving `[K:Q]=4`. It is a totally imaginary quadratic extension of the totally real quadratic field F. No Galois-normality assertion is needed.

## 2. The whole20 algebra is a field

Since `[K:F]=2`,

\[
N_{K/Q}(w)=N_{F/Q}(w)^2
 =\left(\frac1{408146688}\right)^2
 =2^{-16}3^{-26}.                              \tag{3}
\]

There is a place of K above 2 where the normalized valuation of w is not divisible by five: otherwise its norm's 2-adic valuation, the residue-degree-weighted sum, would be divisible by five, contrary to −16.

Let L be K with one root of `v5=w` adjoined. At such a place, extending the valuation to L gives

\[
5\,\operatorname{val}_L(v)=e\,\operatorname{val}_K(w).
\]

Hence 5 divides the ramification index e. Since `[L:K]≤5`, the extension has degree exactly five. Thus `X5-w` is irreducible over K. The saved20 algebra is consequently the field L of degree20, under the inverse formulas `U=34992v5`, `a3=Zv²`. In characteristic zero it is reduced and has twenty geometric embeddings. This establishes field/reducedness rather than assuming it from the vector-space dimension.

## 3. Actual diagonal descent of the accepted parent tensor

Use the original T37 homogeneous coordinates X0,...,X4 and weights

\[
(d_0,d_1,d_2,d_3,d_4)=(0,3,2,4,1).
\]

Make the actual coordinate substitution `Xi=v^di*Yi`, i.e. pull the tensor back by `C(v)=diag(1,v³,v²,v⁴,v)`. Its determinant is `v10=w²`, a unit. Every monomial's weight difference `sum(di*exponent_i)-di-dj` is divisible by five. Consequently the actual transformed tensor P_w has coefficients in Q(w), and its nonzero upper brackets are

```
P01 = w*Y3²
P02 = Y4²
P12 = Y1*Y2 + Y3*Y4/18 - Y0²/(48*w)
P13 = 3*Y1*Y3 + Y0*Y2/w - Y4²/(3*w)
P14 = -Y1*Y4 + Y0*Y3/3
P23 = Y2*Y3 + Y0*Y4/(3*w)
P24 = -3*Y2*Y4 + Y0*Y1 + (w/3)*Y3²
P34 = -2*Y3*Y4 - Y0²/(2*w).
```

Each displayed coefficient follows monomial by monomial from the frozen accepted T37 tensor with all five source coefficients equal to one. This does not introduce independently varying source periods. For every nonzero v the tensor is literally the pullback of that same known point. As a consistency check, w=1/31104 gives precisely the earlier rational tensor in the separate algebraic-plane receipt.

Replacing v by ζv with ζ5=1 is an actual symmetry of the original tensor, given by `diag(ζ^di)`. The five root choices at fixed (U,Z) are therefore actual PGL-equivalent parent/plane pairs, rather than five new cases requiring P4 anchors.

## 4. The entire plane descends to K

The original cell plane is

```
X0=A, X3=B, X4=C,
X1=a0*A+a1*B+a2*C,
X2=a3*A+a4*B+a5*C.
```

In the new coordinates put `A_old=A_new`, `B_old=v4*B_new`, `C_old=v*C_new`. Its coefficients become

\[
(b_0,b_1,b_2,b_3,b_4,b_5)
=(a_0/v^3,\;va_1,\;a_2/v^2,\;a_3/v^2,\;v^2a_4,\;1).
\]

They have the following small expressions over K:

\[
\begin{aligned}
b_0&=\frac{306U+522+12(U+4)Z}{13},\\
b_1&=\frac{(11U+18)Z-12U-48}{702},\\
b_2&=-\frac{12}{13}\big((U+4)Z+6U+11\big),\\
b_3&=Z,\qquad b_4=-\frac{2U+3}{648},\qquad b_5=1.
\end{aligned}                                  \tag{4}
\]

For a direct derivation using only the original six equations, set `S=a4/v³=1889568w+90=18(3U+5)`. Their three linear relations give

\[
b_0=S/3-b_2,\qquad Z=b_1S-b_2/3,
\qquad b_1-3wSb_2+1/18=0.
\]

The first equation in (1) implies `wS²=(U+1)/12` and `wS=-(2U+3)/648`. Hence

\[
b_2=\frac{12(Z+3U+5)}{3U-1},
\qquad(3U-1)^{-1}=-(U+4)/13,
\]

which gives (4) and then b0,b1. All remaining verifications can now be done in degree-four K, using the small literal parent P_w and rank-three plane with rows `(1,0,0),(b0,b1,b2),(Z,b4,1),(0,1,0),(0,0,1)`. No denominator depending on a new invariant has been introduced; the only variable denominator w is already a unit.

The original finite20 checkpoint imposed every original Lagrangian equation and all Gröbner relations. Formula (4) is a hand reduction of that exact source map. A native batch should directly replay all six actual Lagrangian coefficients and the quintic in K before using the reduction for a positive recognition theorem.

## 5. What conjugacy proves, and what it does not

Because K is a field, all four embeddings form one arithmetic Galois orbit. Any rank, nonvanishing, intersection-multiplicity or algebraically defined singularity-type assertion with rational coefficients has the same truth value at all four conjugates. In particular an exact nonzero plane-Jacobian determinant in K is a unit and is nonzero at all twenty original embeddings. Conversely a zero determinant is zero everywhere. This avoids treating the twenty algebraic solutions as separate trial anchors.

Galois conjugacy **alone does not prove PGL equivalence** of the four marked divisors: a continuous modulus can have distinct conjugate values. To conclude one actual source orbit, first recognize the intrinsic marked quartic type and apply an already proved whole-marked-pair uniqueness theorem for that type. The accepted finite lists give precisely such uniqueness for the E6 smooth-flex and A5 smooth-hyperflex candidates. Full Saito-frame and magnetic/source identification remain necessary; the cofactor name or quartic singularity label alone is insufficient.

The useful exact native targets are therefore the small six-by-six Lagrangian-plane derivative determinant and the actual restricted quintic/marked-factor geometry over K. If a relevant plane is an isolated reduced Lagrangian zero, the moving-plane/normal-scaling construction in `T37_OPEN_LAGRANGIAN_CURVE_FRONTIER.md` gives a regular actual parent curve whose magnetic endpoint is the full derivative image. If that derivative is invertible, all six magnetic coefficients are covered. None of these geometric recognition or determinant conclusions has been computed in this note.

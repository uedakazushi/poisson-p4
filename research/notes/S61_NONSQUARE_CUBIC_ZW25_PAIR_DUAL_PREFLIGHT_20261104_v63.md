# The old pair dual after the eight projected-cubic unit pivots

The fixed pair-dual preflight is exact over the FULLU coefficient ring
\(R=\mathbf Q[a,e,d,1/(adFBCJ)]\). It does not add a pivot. The old
pair dual, restricted to the surviving three pair columns, has rank at most
one. All three tested two-by-two determinants vanish identically. This is
a limitation of this particular restriction, not a rank-drop result for
the actual cubic matrix.

The source is the immutable eight-step checkpoint
`S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json`; the old
eight-step Gaussian is not rerun. The new producer is
`code/s61_nonsquare_cubic_zw25_pair_dual_preflight_20261104_v63.py` and its
report and exact scalar certificate are in the matching report directory.
The job took 0.564 seconds, evaluated only the pair-dual scalar ancestors,
four dual identities, six fixed dot products and three fixed determinants,
and used no new point, field specialization, Gaussian, GB or determinant
search. Independent review is pending.

## The source and normalization binding

In the initial projected 25-by-32 matrix, the two rows numbered 8 and 9
are respectively \(q_0L\) and \(1600q_5L\). Their columns 20 through 25
are \(L\operatorname{Sym}^2(y_1,y_2,y_3)\), in the L-frame order.
They are literally the saved `pair2x6.txt`, coefficient by coefficient.
The accepted extractor dual instead uses recursive ascending monomial
order and the unscaled \(q_5\). Thus its transported dual is

\[
 D_L=\operatorname{reverse}_6(D_{\rm actual})
       \begin{pmatrix}1&0\\0&1/1600\end{pmatrix}.
\]

All four entries of \(C_LD_L=I_2\) are checked exactly. The native
historical pair source has an additional \(Cq_5\) normalization; that
normalization is already undone in the accepted `D_actual` roots. An
additional factor of \(C\) must not be inserted here.

After the eight saved pivots, the remaining quadratic source rows are
4, 5 and 9: \(1600q_5x_0\), \(1600q_5y_1\), and
\(1600q_5L\). The surviving original pair columns are 21, 22 and 25,
at residual target positions 16, 17 and 18. Let \(Q_s\) be their
literal three-by-three checkpoint submatrix, and let \(D_s\) be these
three rows of \(D_L\). The certificate saves every coefficient of
\(Q_s,D_s\), and \(F=Q_sD_s\). Row 4 of \(F\) is zero; each entry
in its other two rows is nonzero. All three two-by-two minors of \(F\)
are zero.

## Why this restricted dual cannot give two new pivots

The construction of the accepted pair dual starts with the \(q_0\)
coefficient \(5a\) in L-frame pair column 23. Outside that distinguished
column the two columns of \(D_L\) satisfy

\[
 (D_L)_{j,0}=\gamma(D_L)_{j,1},\qquad
 \gamma=-\frac{(C_L)_{1,3}}{(C_L)_{0,3}}
       =-\frac{200a}{3}(35ad-74e-2).
\]

This follows directly from the saved pair-dual formula: the first column
is its distinguished pivot vector minus the pivot coefficient times the
second column. Column 23 is among the seven saved quadratic pivots.
The remaining rows of \(D_L\), in particular 21, 22 and 25, therefore
have proportional columns. Restriction cannot preserve a two-dimensional
dual image. The literal six-dot check independently confirms the resulting
three zero determinants. A different completion or a dual for the whole
ten-row quadratic block is needed; no Bézout identity for two new pivots
is supplied by this computation.

For comparison, an actual unit two-minor of \(F\) would have sufficed
without assuming a free projective complement. For such two chosen rows
put \(L=F_I^{-1}(Q_s)_I\). Then \(LD_s=I_2\), and the three columns
\([D_s,\,L_0\mathbin\times L_1]\) have determinant
\(\det(LD_s)=1\). Its inverse has first two rows \(L\) and last row
the cross product of the columns of \(D_s\). This is a valid explicit
GL3 completion criterion, but its hypothesis fails for the fixed dual
above. A general split injection by itself is not an explicit free
target complement or a computed unimodular change of coordinates.

## The actual quadratic and heavy quotients remain separate

The accepted original Z=W quadratic-product map has rank 19 everywhere
on FULLU. Removing its nine monomial unit pivots and its one literal
zero overlap leaves the actual ten-row quadratic block with rank 10
everywhere. Seven of the saved eight pivots are quadratic pivots; the
eighth uses the first heavy coordinate of genuine source G1. Every
quadratic row has zero heavy coordinates, so that last step does not
modify the three remaining quadratic rows. They still have rank 3 on
FULLU. The failed restricted dual does not contradict this fact. It
does not compute a global three-row right inverse or a free basis of
the complementary target bundle.

The five heavy coordinates retained by Z=W projection, together with
\(\psi=-wH_Z+zH_W\), are a complete heavy quotient modulo the known
source line \(C_\delta=\delta(G_{61}-G_{20})\). In the accepted unit-pair
packet, \(c_Zz+c_Ww=\eta\) is a FULLU unit; consequently

\[
 \begin{pmatrix}c_Z/\eta&c_W/\eta\\-w&z\end{pmatrix}
\]

has determinant one. Thus this six-coordinate heavy quotient omits no
heavy direction other than the explicitly split line. After the saved
G1 pivot it has five remaining coordinates. Any further attempt should
distinguish the rank-three residual pure quadratic subbundle from the
fourteen genuine source rows; after a justified quotient by that
subbundle, their target has rank 21 (pure rank 16 plus heavy rank 5).
An explicit coefficient dual or a local atlas is needed to perform
that quotient constructively. No global basis is inferred from rank.

Whole-FULLU cubic rank 57, full quartic rank 239, and the known uniform
kernel count 38 with gap 42 are unchanged. The accepted cubic-rank-57
result still applies on the additional open D(Delta34). Unknown factors,
including Delta34, are not inverted in this preflight.

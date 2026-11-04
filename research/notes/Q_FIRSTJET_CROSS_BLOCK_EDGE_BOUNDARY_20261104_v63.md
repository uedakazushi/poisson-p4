# cross block の edge 境界の小さい普遍恒等式

2026-11-04。root 提案の ordered normal form における edge 条件を独立に検証した。**元の NEW25 source への代入も、新しい点の評価も行っていない。** 本稿は幾何的条件の普遍式だけで、NEW25 方程式の主張ではない。

前稿の D=diag(x,−x,y,−y)、b=e=1 の開 chart を使う。ここで cross block は

\[
 \begin{pmatrix}g&h\\i&j\end{pmatrix},\qquad
 u_e=gj,\quad v_e=hi,\quad m=u_e-v_e.
\]

uₑ,vₑ は block の積であり、元の三変数 Poisson 族の source u,v とは別である。

S=x²+y²=F₄/20、Δ=x²−y²、D₀=Δ²=p₀=S²−4p₄ とする。

\[
 p_3=-2ax\Delta-2x^2S,\quad A=p_3+S^2,
\]
\[
 B_2=-4a^2x^2-2ax^3-6axy^2+3x^4+y^4.
\]

全係数で

\[
 A=-\Delta(2ax+S),\qquad
 B_2=S^2-p_3+D_0-A^2/D_0
\]

を確認した。κ=p₂−B₂ と置き、

\[
 K_{16}=D_0(p_2-S^2+p_3-D_0)+A^2
\]

とすると κ=K₁₆/D₀。前稿の二つの積の式から

\[
 u_ev_e=\frac{\kappa^2-2S\kappa m+D_0m^2}{16p_4}.
\]

さらに普遍的な m=N₁₂/(2D₀) を使うと

\[
 G_{32}=D_0N_{12}^2-4SK_{16}N_{12}+4K_{16}^2
\]

について、正確な cofactor は

\[
 \boxed{G_{32}=64p_4D_0^2u_ev_e.}
\]

これを元の x,y,a,p₂,m の全有理式として検査した。従って **D₀p₄≠0 の表示 chart 上**で G₃₂=0 は uₑvₑ=ghij=0 と同値である。すなわち四つの cross-entry edge の和集合を表す。

S=F₄/20、D₀=F₄²/400−4p₄を戻すと G₃₂は F₄,p₂,p₃,p₄,N₁₂ の36項多項式で、βについて斉次32次。N₁₂² の係数 D₀ は非零である。六つの universal generator の代数的独立性により、全W上で零な恒等式ではない。

一方、実際の β 空間へ引き戻した後の factorization、D₀p₄=0上の境界成分、NEW25 との一致は本稿で決めない。特に degree24の41部分空間の負の結果から degree32候補が正しいとは推論しない。N₁₂や pᵢ の NEW25 source 値・既存点の値も今回代入していない。

保存：`code/q_firstjet_cross_block_edge_boundary_20261104_v63.py`、
`reports/Q_FIRSTJET_CROSS_BLOCK_EDGE_BOUNDARY_20261104_v63/{report.json,certificate.json}`。

証明書は全恒等式と primitive integer にした36項を保持する。実行約0.30秒、60秒・2GB上限。次の source 消滅検査は未実行である。

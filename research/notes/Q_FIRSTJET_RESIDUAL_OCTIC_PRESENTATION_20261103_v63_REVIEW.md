# 固定 base 八次表示の独立レビュー

2026-11-04。`Q_FIRSTJET_RESIDUAL_OCTIC_PRESENTATION_20261103_v63` を独立受理する。選ばれた実像成分の正確な次数・種数や NEW25 multisection の次数は受理範囲に含まない。

著者の P₈ の因数式を出発点とせず、二つの T 二次式 Q₁,Q₂ の**4×4 Sylvester 行列の全24置換 determinant**を別の QQ 多項式環で計算した。m=u−v は最初から消去する。結果は厳密に

\[
 \operatorname{Res}_T(Q_1,Q_2)=-cP_8.
\]

ここから作り直した P₈ の九係数は保存済みの全係数と一致し、leading=x²v⁴、constant=x²f⁴u⁴も確認した。次数は一般係数体で正確に8である。

同じ独立多項式環で、E=LT+M、a₁=−cyR として次を全係数で確認した。

\[
 E=yRQ_2+AQ_1,
\]
\[
 yRP_8=L^2Q_1-a_1E^2-(bL-2a_1M)E,
\]
\[
 a_1M^2-bML+c_1L^2=yRP_8,
\]
\[
 yR(a_2M^2-b_2ML+c_2L^2)=-AyRP_8.
\]

これにより T=−M/L の両方向の ideal identity を検証した。実際の固定 base の二つの (2,2) 式への全係数 binding も新たに比較した。

一般係数体での単元証明も受理する。P₈ の定数項は非零なので R は単元。P₈=L=0 なら q=0。q は K₀[R] に属するので、その根は K₀ 上代数的である。一方 K₁ は K₀ 上超越的で、L=L₀−K₁RA。したがって共通根では A=0となる。しかし独立に計算した

\[
 \operatorname{Res}_R(A,q)=f^2(u+v)^2+4d^2fuv
\]

は、実係数へ戻すと m の4次式で leading=(x²+y²)²/(4x²y²)≠0。従って A,q は一般係数体で共通根を持たず、gcd(P₈,L)=1。q も同じ議論で P₈ と互いに素だから単元で、T=−ymRq/L も単元となる。これは一般の関数体についての証明であり、全特殊化での非消失ではない。

まず K[R]/(P₈) 内で L,q,R の逆を作り、上の恒等式から Laurent 商からこの環への写像を得る。逆写像では ideal identity から P₈=0となり、得られた L の逆が写るため T の回収も成立する。この順序なら、元の Laurent 商で先に L が単元だと仮定する循環はない。

したがって

\[
 K[R^{\pm1},T^{\pm1}]/(Q_1,Q_2)\cong K[R]/(P_8)
\]

を全スキーム代数として受理する。右辺の長さ8は全交点代数の値であり、被約性、既約性、八つの幾何点、特定の実成分が全長さを占めることを与えない。NEW25 の base 上の multisection 次数とも別である。

独立コード `code/q_firstjet_residual_octic_presentation_review_20261103_v63.py`、報告と証明書は `reports/Q_FIRSTJET_RESIDUAL_OCTIC_PRESENTATION_20261103_v63_REVIEW/`。実行約0.76秒、追加点・GB・高次数不変式評価はなし。著者のコード・ノート・全依存と証明書の hash を照合した。

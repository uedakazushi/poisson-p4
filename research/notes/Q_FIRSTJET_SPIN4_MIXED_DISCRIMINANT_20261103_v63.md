# Spin(4) 切片の二形式と五交点の混合判別式

2026-11-03。前 packet `Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63` の gamma-traceless vector-spinor 模型を使う。**Spin(4) 切片の厳密な同定、混合判別式の170-node回路、有理的な Sp₄ 降下まで**を証明する。NEW25 上での消滅、24次の多項式への延長、既存の曲線判別式との同一性は証明していない。NEW25 の点への評価は行っていない。

## 1. 全12次元切片と二形式

前 packet の Γ₀=diag(1,1,−1,−1) に対し、spinor を A⊕B、dim A=dim B=2 と分解する。ベクトル空間 E₅ の最初の基底 e₀ の直交補空間は A⊗B である。Ψ の最初の行を0とした切片を L とする。

cross 2形式の行を (ij), i,j=0,1 と書き、残りの4行を Ψᵢⱼ,ₖ とする。Clifford trace 条件はちょうど

\[
 \Psi_{1j,0}=\Psi_{0j,1},\qquad
 \Psi_{i1,2}=\Psi_{i0,3}.
\]

従って

\[
 U(a,b)=\sum_{i,j,k=0}^1\Psi_{ij,k}\,a_i a_k b_j,
 \qquad
 V(a,b)=\sum_{i,j,k=0}^1\Psi_{ij,k+2}\,a_i b_j b_k
\]

により

\[
 L\simeq (\operatorname{Sym}^2A\otimes B)
       \oplus(A\otimes\operatorname{Sym}^2B),
 \qquad \dim L=6+6=12.
\]

U の a₀a₁、V の b₀b₁ の係数は Ψ の対応する係数の2倍である。計算コードはこの係数を1/2で戻す。

`certificate.json` は12基底の**元の24係数への全埋め込み**と受理済み16基底への全埋め込みを保存する。block diagonal な全6個の sp₂⊕sp₂ 生成元で、元の物理的 Lie 微分と U,V の変数置換微分が全係数で一致する。従って、単に同じ次元の表現へ名前を付けた同定ではない。

## 2. V が滑らかな場合の完全な binary quintic

\[
 U=u_0(b)a_0^2+u_1(b)a_0a_1+u_2(b)a_1^2,
 \quad V=a_0v_0(b)+a_1v_1(b)
\]

と書く。uᵢ は binary linear、vᵢ は binary quadratic である。R=Res₂(v₀,v₁)≠0 なら、V=0 は

\[
 b\longmapsto (a,b)=((-v_1(b),v_0(b)),b)
\]

で表される滑らかな (1,2) 曲線で、P¹ と同型である。U の制限は

\[
 F_5(b)=u_0(b)v_1(b)^2-u_1(b)v_0(b)v_1(b)+u_2(b)v_0(b)^2.
\]

これは (U,V) に関して二重次数 (1,2) の binary quintic である。6個の u 係数から6個の quintic 係数への線形写像の行列式は、**全係数計算で R³** と確認した。従って R≠0 では全 binary quintic を得る。幾何的にも、制限列の核と余核は O(1,−1) の H⁰,H¹ であり、両方0である。

交点数は (2,1)·(1,2)=2·2+1·1=5。R≠0 のとき U,V の共通零点が reduced でないことは、F₅ に重根があることと同値である。

## 3. 余分な resultant 平方因子と多項式への延長

binary quintic の判別式は係数について8次なので

\[
 \operatorname{Disc}_5(F_5)
\]

の二重次数は (8,16) である。しかし V が基点を持つ divisor R=0 上では、v₀,v₁ が共通零点を持ち、F₅ は自動的に重根を持つ。この余分な因子は**ちょうど R²**である。

証明する。quadratic resultant R は既約である。これは P¹ の共通零点を記した incidence の既約な hypersurface 像である。その一般点の DVR で、共通零点の座標を t=0、uniformizer を ε とする。

\[
 v_0=t a(t)+\varepsilon p(t),\qquad
 v_1=t b(t)+\varepsilon q(t)
\]

と書けるから、F₅ の t=0 付近の二次 Weierstrass 因子は

\[
 t^2+\varepsilon B(\varepsilon)t+\varepsilon^2C(\varepsilon)
\]

の形になる。残りの三根は一般には単純である。従って判別式の R-adic valuation は少なくとも2である。これは U,V 全係数の多項式環が UFD であることにより R² の多項式としての整除性を与える。

次数が正確に2であることには、小さい厳密な fixture を使う。

\[
 v_0=t,\ v_1=t^2+\varepsilon,\quad
 (u_0,u_1,u_2)=(1+t,3+2t,7+5t).
\]

homogeneous quadratic の resultant は ε であり、全 quintic 判別式は

\[
 \varepsilon^2(4096\varepsilon^6+26368\varepsilon^5
 +144080\varepsilon^4+674816\varepsilon^3
 +1554120\varepsilon^2+1825600\varepsilon+1466325).
\]

最後の定数は非零だから R³ は整除しない。従って

\[
 \boxed{D(U,V)=\operatorname{Disc}_5(F_5)/R^2}
\]

は二重次数 (8,8)、合計16の多項式である。

R を可逆にすると U→F₅ は線形同型で、binary quintic 判別式は既約である。上の exact valuation により D は R 因子を持たないから、D は全係数環でも既約である。この divisor は、P¹×P¹ 上の五交点の非横断性の閉包である。V 自体が特異なすべての点で同じ単純な幾何的説明を主張してはいない。

回路は、6個の quintic 係数、多項式 R、一般 quintic 判別式の全59項、R⁻²を組み合わせる。**12変数の D の全展開は行わない。** 回路は R≠0 上の式を保存し、R=0 への多項式延長はこの証明による。単なる有理回路の数値的な整除性を主張しているのではない。

## 4. 実 normalizer と有理的な Sp₄ 降下

Sp₄ の元

\[
 g=\begin{pmatrix}0&I_2\\-I_2&0\end{pmatrix}
\]

は e₀ の normal line を保ち、その向きを反転する。元の24係数への全作用を計算すると、切片上では

\[
 (U(a,b),V(a,b))\longmapsto
 (-V(b,a),U(b,a)).
\]

δ=det(Ψ の cross 4×4 行列) と置くと δ↦−δ である。一方 D の零点条件は U,V と二つの P¹ を入れ替えても同じである。D の既約性と同じ次数から、g*D=cD。固定した整数係数の fixture で両辺を厳密に評価し、c=1を定めた。元の R と変換後の R は50,37で、D の値は両方53268630696となる。この fixture は既約性から既に分かる比例定数の規格化だけに使い、標本から多項式恒等式を推測していない。

したがって D は e₀ line の stabilizer N 全体で不変である。N の恒等成分は Spin(4) で、上の g がもう一方の成分を代表する。

前 packet の C に対し p₀=det C とする。Ψ の4×4 maximal minors を交代符号付きの列 χ として、normal vector を n=Q⁻¹χ と置く。Cauchy–Binet により

\[
 Q(n,n)=16\det(\Psi^{\mathsf T}Q\Psi)=16p_0.
\]

切片上では n=δe₀、従って

\[
 p_0|_L=\delta^2/16.
\]

p₀≠0 では Ψ は rank4 で、その像の直交補 normal line は一意な nonisotropic line である。Spin(5) はこれらの line に推移的に作用する。line を記した incidence を使うと

\[
 W_{p_0\ne0}\simeq G\times^N L_{\delta\ne0}.
\]

忘却写像の逆は maximal minors による normal line なので、この主張は単なる一般点の軌道同定ではなく、この開集合上の同型である。従って N 不変な D は **W の p₀≠0 上の regular Sp₄ 不変式**に降下する。β の scalar 倍と降下が可換なので、その斉次次数は16である。

W は affine space、p₀ は8次不変式だから、この有理不変式は

\[
 \mathcal D=P/p_0^k,\qquad
 P\in\mathbb C[W]^{Sp_4},\quad \deg P=16+8k
\]

と書ける。P の不変性は分母 p₀ の不変性から従う。したがって **k=1 の24次 lift は次数と normalizer の両方に整合する。**

ただし、今回の証明は k≤1 を与えない。極の位数を generic p₀=0 上で計算するか、24次の多項式 P を別に構成して切片への制限 P|L=(δ²/16)D を証明する必要がある。k=0 の可能性も、この周囲の構成だけからは除外していない。NEW25 上の消滅はさらに別の未証明条件である。

## 保存と次の有限な問題

* `code/q_firstjet_spin4_mixed_discriminant_20261103_v63.py`
* `reports/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63/` の全 physical embedding、6作用行列、normalizer、restriction determinant、59項公式、170-node回路、全 fixture。

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_firstjet_spin4_mixed_discriminant_20261103_v63.py \
  --output work/spin4-mixed-replay
```

60秒・2GB上限。現在残る具体的な問題は、(i) p₀ に沿う極の位数、(ii) NEW25 の実 source と五交点非横断条件の対応、(iii) 既存24次曲線判別式との一致または相違である。未証明の段階で NEW25 の点を使った高次数探索や判別式候補の認定を行わない。

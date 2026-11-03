# NEW25 は全8次 spectral base を優越する

2026-11-03。**既に使った同じ QQ 点[0] 一点で、一回の exact source 微分**を実行した。五つの8次不変式の NEW25 source 上の微分 rank は5である。従って NEW25 は canonical adjoint pencil の5次元係数空間を優越し、**その base だけから作る非自明な divisor は NEW25 を定義できない。** 24次以上の多項式評価や追加点は使わない。

## actual source からの全微分を固定する

対象は受理済み `Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz` の6079-node mapである。入力は n₀,n₁,n₂、全 scalar、評価点4変数、Sp₄ big cell10変数の18個。点は `Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63` と同一の既存点であり、新しく選んだ点ではない。

全 node の Fraction forward AD により、全61出力の値と61×18 source 微分を再構成した。value pivot と Pfaffian は既存 QQ 記録と一致し、全逆数 node は定義される。S61 の primitive 一次部分から受理済み16基底への**全線形変換**を作り直し、β₁₆の全値が既存記録と一致することを確認した。これにより16×18 source 微分を原 Poisson 族へ結び付けている。

前 packet の C,D の20×16微分を使う。全8次不変式の一つの基底は

\[
 I=(I_0,I_1,I_2,I_3,I_4)
  =\bigl(\operatorname{tr}C^4,\operatorname{tr}(C^3D),
      \operatorname{tr}(C^2D^2),\operatorname{tr}(CDCD),
      \operatorname{tr}(CD^3)\bigr).
\]

各 word の trace product rule を4×4行列で直接計算し、5×20 pair gradient、5×16 β gradient、5×18 actual source 微分を得た。10本の Sp₄ source 列が全て厳密に0になることも確認した。

## 非零 QQ minor と優越性

source の列

\[
 (n_0,n_1,n_2,\mathrm{scale},\mathrm{base1})
\]

を取る5×5 minor は QQ 上非零である。完全な有理数を証明書に保存した。分子461 bit、分母329 bit、補助的な32003還元値は29835である。rank の根拠は QQ minor であり、modular rank のみではない。

この18変数 source は actual NEW25 の first-jet 像を与え、その像は既約である。標数0の微分判定により、I₀,…,I₄ はその像上で代数的独立であり、

\[
 \boxed{B_{\rm NEW}\longrightarrow\mathbb A^5_I
        \text{ は dominant}.}
\]

低 rank を一般的な関係式と読み替える操作はしていない。今回は上限5に達する下界 minor による優越性の証明である。

## canonical pencil 係数への全係数変換

\[
 \det(tC+uD)=\sum_{k=0}^4p_k t^{4-k}u^k
\]

とする。四式と Newton 恒等式から

\[
 \begin{pmatrix}p_0\\p_1\\p_2\\p_3\\p_4\end{pmatrix}
 =\begin{pmatrix}
 0&0&-1/2&0&-1\\
 -1/2&-1&1&0&2\\
 1/2&0&-2&-1/2&-2\\
 -1/4&0&1/2&0&0\\
 1/32&0&1/16&0&1/8
 \end{pmatrix}
 \begin{pmatrix}I_0\\I_1\\I_2\\I_3\\I_4\end{pmatrix}.
\]

この定数行列の行列式は **−1/128**。保存済みの全8次係数辞書を比較し、五つ全ての恒等式を全単項式で確認した。この変換により

\[
 B_{\rm NEW}\longrightarrow\mathbb A^5_{p_0,\ldots,p_4}
\]

も dominant である。従って p₀,…,p₄ は NEW25 上でも代数的独立である。

特性多項式の全 base は F₄ も含む。

\[
 \mathcal B=\{(F_4,p_0,\ldots,p_4):
                  F_4^2=400p_0+1600p_4\}.
\]

𝓑 は p₀を消去して A⁵と同型で、p-spaceへの射影は有限二重写像である。NEW25 の像の閉包は既約で、p-spaceを優越するため次元5を持つ。既約な5次元の𝓑内ではその閉包は𝓑全体となる。従って

\[
 \boxed{B_{\rm NEW}\longrightarrow\mathcal B
             \text{ も dominant}.}
\]

これは F₄ と4個の独立な8次生成元が、周囲だけでなく NEW25 上でも多項式部分環を作ることを与える。全 invariant 環がこの部分環に等しいとは述べない。

## 除外できる候補のクラスと残る一方向

𝓑 の任意の真の閉部分集合の逆像は、NEW25 全体を含めない。従って adjoint pencil の行列式 binary quartic、特性多項式、二重被覆の branch quartic などの**係数だけ**から作る条件で、𝓑上の真の divisor を定めるものは NEW25 方程式にならない。全 base 上で恒等的に零となる式はそもそも divisor ではなく、この結論と区別する。

これは一つの高次数判別式を評価して非零と確認する方法ではなく、全 base への優越性により、そのような base-only divisor のクラスを一度に除外する。

前 packet により全 W の Sp₄ rational quotient は6次元で、一般軌道は二共変量対の同時共役で分類できる。従ってその𝓑への一般 fiber は1次元である。一方 NEW25 の rational quotient は5次元で𝓑を優越するから、一般関数体の拡大は有限となる。

すなわち NEW25 の未知 divisor は、この残る1次元方向での**generically finite multisection**として現れる。multisection の次数、一般 fiber の種数、商の有理性、全域での有限性や固有性は未決である。七座標正規形の周囲 hypersurface だけから NEW25 条件が決まるとも主張しない。

もし NEW25 方程式が最小可能次数24で存在しても、F₄と全8次不変式だけの式にはならない。追加の invariant・fiber coordinate に依存する必要がある。今回24次の基底や方程式を評価・探索してはいない。

## 保存と再現

* `code/q_new_family_octic_base_dominance_20261103_v63.py`
* `reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63/{report.json,certificate.json.gz}`

証明書は同じ point、全61出力、61×18 actual AD、primitive 変換、16×18 β source 微分、5×20/5×16/5×18の全微分、非零 minor、全係数で確認した5×5変換を保持する。

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_octic_base_dominance_20261103_v63.py \
  --output work/new25-octic-base-replay
```

一回の同点微分検査、上限60秒・2GB、実行約2.49秒。新規 β 点、24次以上の多項式評価、GB、新しい NEW25 方程式は含まない。

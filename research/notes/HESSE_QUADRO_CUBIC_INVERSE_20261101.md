# Hesse 五二次式の普遍的な三次逆写像

2026-11-01。以下は楕円二族を比較するための birational map の厳密な証明書である。
**primitive 一次ジェットの同変性はこのノートの主張に含めない。**

## 1. 二次写像と三次逆写像

添字は mod 5、λ≠0 とする。

\[
 Q_i(x)=x_i^2+\lambda x_{i-1}x_{i+1}-\lambda^{-1}x_{i-2}x_{i+2}.
\]

三次式 R0 を

\[
\begin{aligned}
 R_0(y)={}&\lambda^2y_0^3-\lambda^4y_0y_1y_4-y_0y_2y_3\\
 &-\lambda^3(y_1^2y_3+y_2y_4^2)
 +\lambda(y_1y_2^2+y_3^2y_4)
\end{aligned}
\]

とし、Ri は変数の添字を i だけ巡回させた式とする。
すると五次斉次式 gλ(x),hλ(y) があって、QQ(λ) 上の多項式恒等式

\[
 R(Q(x))=g_\lambda(x)x,\qquad Q(R(y))=h_\lambda(y)y \tag{1}
\]

が成立する。`gλ(e0)=λ²`, `hλ(e0)=λ⁴` である。
従って `gλ≠0`, `hλ≠0` の非空な稠密開集合上で [Q] と [R] は互いに逆であり、
[Q] は birational である。全展開した二つの五次式は証明書に保存した。

この結論は λ≠0 なら恒等式として成立し、楕円基底曲線の滑らかさを追加の前提としない。
Poisson 族への応用では、それぞれの族に必要な滑らかさ等の開条件を別に課す。
(1)から `h(Qx)=g(x)²`、`g(Ry)=h(y)³` も得られるため、開集合の対応は両方向で成立する。
これにより、順方向写像だけを用いて解を失ったり増やしたりする問題を避けられる。

## 2. 逆写像の導出

`t=(0,−λ,1,−1,λ)` として、交代行列

\[
 \Phi_{ij}(x)=t_{j-i}x_{3(i+j)}
\]

を作る。直接展開すると `Φ(x)Q(x)=0` となる。
`Φ(x)y=A(y)x` によって A を定めると

\[
 A(y)=\begin{pmatrix}
0&y_2&\lambda y_4&-\lambda y_1&-y_3\\
-y_4&0&y_3&\lambda y_0&-\lambda y_2\\
-\lambda y_3&-y_0&0&y_4&\lambda y_1\\
\lambda y_2&-\lambda y_4&-y_1&0&y_0\\
y_1&\lambda y_3&-\lambda y_0&-y_2&0
\end{pmatrix}.
\]

交代性から `yᵀA(y)=0`。e0 における第0行・列を除いた小行列式は `−λ²` なので、
一般の A の階数は4である。`−adj(A)` の第0列を y0 で割ると、上の三次ベクトル R を得る。
割り算の結果が全て多項式であること、`A R=0`、(1) を記号計算で厳密に確認した。

この導出は一般点での核の一意性を説明し、(1)はさらに全パラメータの恒等式を与える。
単なる λ の標本や Jacobian の標本から birationality を推論していない。

## 3. 文献との関係と応用上の境界

楕円正規五次曲線の五二次式による quadro-cubic Cremona 変換は既知である。
Polishchuk, *Feigin–Odesskii brackets, syzygies, and Cremona transformations*,
Theorem B と §3.2 は一般の Pfaffian 表示から逆写像を構成し、n=5では
q(5,1) と q(5,2) の Poisson birational map を与える。
[確認した一次資料本文](https://arxiv.org/html/2204.09152)。
今回の式は保存済み Hesse 規約における短い具体的証明書であり、この既知の変換を新規とは主張しない。

非線形な Poisson 同値であることだけから、S61 の primitive 一次部分を保存するとはいえない。
楕円二成分の次数比較に用いるには、正規化の倍率と線形微分を含む一次ジェット恒等式を
別に証明する必要がある。その恒等式と一般開集合の対応が確立された場合にのみ、
この逆写像証明書をパラメータファイバーの比較へ適用できる。

## 4. 再現

```sh
.venv/bin/python code/hesse_quadro_cubic_inverse_20261101.py \
  --output reports/new-hesse-inverse.json
```

保存済み結果は `reports/HESSE_QUADRO_CUBIC_INVERSE_20261101.json`。
Q、Φ、A、R、g、h、階数4小行列式と入力ハッシュを保持する。

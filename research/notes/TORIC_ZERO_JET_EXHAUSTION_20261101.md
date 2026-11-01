# 実際の T01 の零一次部分ファイバーの被約な尽くし

2026年11月1日。基礎体は \(\mathbb C\)。

## 定理

実際の既約トーリック成分 T01 の S61 射影における零一次部分ファイバーについて

\[
 \boxed{\ (T01\cap\{\beta=0\})_{\mathrm{red}}=\overline{AB}\ }
\]

が成り立つ。右辺は旧 compatible quadratic 分類の nil22 AB 成分である。
従って被約ファイバーは既約で、次元は正確に十である。
これはファイバーのスキーム構造、冪零部分、定義イデアルの生成系を決める定理ではない。
全 S61 の一般ファイバーの尽くしとも別の結果である。

証明の残っていた部分は、零モジュラーの C0、L0 の境界を含む全点の処理である。
本ノートでは、その全パラメータ族について strict 安定化次元が三以下になることを示す。
従って T01 に必要な四次元可換 strict 部分代数を持てず、R0 以外の点は除外される。

## 1. 既に確定した部分と全点被覆

以下を使う。

1. [対数型退化](TORIC_ZERO_JET_FIBER_20261101.md)は
   \(\overline{AB}\subset F_0\) を任意の AB 係数で示した。
2. [非零トレースの全点判定](TORIC_ZERO_JET_NONZERO_TRACE_20261101.md)は
   \(F_0\cap\{s\ne0\}\subset\overline{AB}\) を示した。
   最小型では全固定トレース表示の Pfaffian から M0、M1 に戻し、
   [二族の全点判定](TORIC_ZERO_JET_MINIMAL_20261101.md)を使う。
   nil22 型では完全分解条件が AB 正規形を直接強制する。
3. [R0 の判定](TORIC_ZERO_JET_R0_20261101.md)は
   \(F_0\cap R0\subset\overline{AB}\) を示した。
4. 旧 v07 の [被約零モジュラー分類](../../legacy/v07/proofs/zero_reduced_components_theorem.txt)
   は全点について \(Q_0=R0\cup C0\cup L0\) を与える。

従って C0、L0 の R0 外の点を除けばよい。
その全境界を保持する proper incidence と正常形は
[C0 ノート](S61_C0_INFINITY_20261101.md)、[L0 ノート](S61_L0_INFINITY_20261101.md)
に証明され、[独立点検](S61_ZERO_MODULAR_INFINITY_REVIEW_20261101.md)も済んでいる。
ここではそれらの**点集合の全被覆**を使い、一般 \(\beta\) の無限遠に関する
結論を零一次部分のアフィンファイバーへ転用するのではない。

## 2. 零モジュラー型における affine strict stabilizer

R0、C0、L0 はすべて \(\operatorname{Pf}(Q)=0\) を満たす。
従って \(P=J+Q\) の射影退化因子は \(x_0^5\)。
strict 射影自己同型の Lie 代数はその唯一の超平面を保ち、アフィン座標では
\(v+A\mathbf y\) である。次数ごとの条件は

\[
A\in\mathfrak{sp}_4,\qquad \partial_vQ=0,\qquad A\cdot Q=0.
\]

従って

\[
 \dim\mathfrak{aut}_{\mathrm{strict}}(J+Q)
 =10-\operatorname{rank}(\mathfrak{sp}_4\to Q_{45})
   +\dim\ker(v\mapsto\partial_vQ).
\]

非零スカラー倍 \(Q\mapsto cQ\) はこの Lie 代数を変えない。
以下の射影正常形のスカラーを一にできるのは、この直接の次数別条件による。
スカラーを変える射影同型が存在すると仮定してはいない。

## 3. L0 の全境界と二つの正常形

L0 の全非零点は、少なくとも一つの非零行列 \(U\) について

\[
 U^2=0,\qquad U^tJU=0,\qquad (U\mathbf y)\wedge Q=0
\]

を満たす。これは射影化した行列 U の空間からの proper incidence なので、
密な族だけでなく全境界を含む。
rank 一の U は Q に定数因子を与え、既に R0 に入る。
rank 二では \(\operatorname{im}U=\ker U\) は Lagrangian 二平面。
可逆二次双線形形式の合同分類により、二つの型

\[
 U=z\partial_y+t x\partial_w\quad(t\ne0),\qquad
 U=(x+z)\partial_y-x\partial_w
\]

だけを考えればよい。
零トレース・compatibility・接条件の全線形解は四次元である。
そのうち一つの係数が零なら \(Q=F(x,z)\partial_y\wedge\partial_w\) で R0 に入る。
それ以外の係数を一にして、任意の二次式 \(F=f_0x^2+f_1xz+f_2z^2\) を残すと、
全正常形は以下の通り。

非対角型：

\[
\begin{aligned}
Q^{xy}&=-xz,&Q^{xz}&=0,&Q^{xw}&=-t x^2,\\
Q^{yz}&=-z^2/t,&Q^{zw}&=xz,&
Q^{yw}&=(2t-1)xy+(2/t-1)zw+F(x,z).
\end{aligned}
\]

例外型：

\[
\begin{aligned}
Q^{xy}&=-x(x+z),&Q^{xz}&=0,&Q^{xw}&=x^2,\\
Q^{yz}&=(x+z)^2,&Q^{zw}&=x(x+z),&
Q^{yw}&=-3xy-4xw-3zw+F(x,z).
\end{aligned}
\]

ここで接条件だけを使った全線形核は、L0 の必要な包含先である。
全核の点が L0 に属するという逆向きの主張は不要である。

## 4. 全 t、全 F における strict 次元上界

位置座標 \((x,z)\)、運動量 \((y,w)\) の順で

\[
\mathfrak{sp}_4=
\left\{\begin{pmatrix}a&b\\c&-a^t\end{pmatrix}:b=b^t,c=c^t\right\}
\]

と分ける。列群を b の三方向、a の四方向、c の三方向の順とする。
係数行を選ぶと、作用行列はこの列群に関してブロック下三角になる。
各対角ブロックには F が現れない。証明書は全 \(60\times10\) 行列と
選んだ行・列を保存し、後の列群の消滅も検査する。

### 非対角型

- b ブロックには定数 \(-2\) の三次 minor があり、常に rank 三。
- a ブロックには \((t+1)^2/t\) の三次 minor がある。
  \(t\ne-1\) では rank は少なくとも三。\(t=-1\) でも成分一があり rank は少なくとも一。
- c ブロックは正確に

\[
\operatorname{diag}\left(1-3t,-\frac{2(t-1)^2}{t},\frac{t-3}{t}\right).
\]

三つの零点 \(1/3,1,3\) は異なり、任意の \(t\ne0\) で rank は少なくとも二。
\(t=-1\) では rank 三である。

従って \(\mathfrak{sp}_4\) 作用の rank は \(t\ne-1\) で少なくとも八、
\(t=-1\) では少なくとも七。

平行移動核はさらに直接読める。\(Q^{xw}\)、\(Q^{yz}\) から x、z 成分は零。
残る条件は

\[
(2t-1)v_y=0,\qquad (2/t-1)v_w=0.
\]

両係数が同時に零になることはない。従って核次元は高々一で、
例外値は \(t=1/2,2\) だけ。\(t=-1\) では核は零である。
全ての t、F について strict 安定化次元は高々三になる。

### 例外型

b、a、c のブロックはそれぞれ rank 三で、c の行列式は \(128\)。
従って \(\mathfrak{sp}_4\) 作用の rank は少なくとも九。
平行移動核は零であり、証明書には F に依存しない四次 minor がある。
strict 安定化次元は高々一。

したがって L0 の R0 外の全点は、T01 の strict 安定化次元四以上という
必要条件に反する。この段階では可換性をさらに使う必要さえない。

## 5. C0 の全境界

C0 ノートの基点を持たない可換二平面の射影写像により、C0 の全境界も保持される。
その初等的 Jordan 場合分けから、R0 外の点はスカラーを除いて高々二つの
\(\mathrm{Sp}_4\) 軌道に入る。

\[
Q_c=(x\partial_x-y\partial_y)\wedge(z\partial_z-w\partial_w),
\]

\[
Q_j=(-x\partial_x+y\partial_y+z\partial_z-w\partial_w)
     \wedge(z\partial_y+x\partial_w).
\]

\(Q_j\) は上の非対角型の \(t=1,F=0\) と正確に一致し、既に除外した。
\(Q_c\) は \(\mathfrak{sp}_4\) 軌道 rank 八、平行移動作用 rank 四で、
strict 安定化次元は二。これも T01 に入らない。
非零スカラー倍にも同じ計算が適用される。

## 6. 結論と限界

零トレースの全点は R0、C0、L0 で覆われる。
C0、L0 の R0 外は上の一様な strict 条件で除外され、R0 内は既に AB に入る。
非零トレースも全点が AB に入ることを別途示した。
逆の包含は AB の普遍的対数型退化によるので、定理を得る。

これは実際の T01 の零一次部分の**被約台の完全決定**である。
従来の有限可換代数模型はこのファイバーで原点しか与えないため、十次元の
全非零部分を欠くことも明確になる。ファイバーの冪零構造や、全 S61 の一般
ファイバーに未知の延長が存在しないことは、本定理からは従わない。

## 再現

```
.venv/bin/python code/toric_zero_jet_exhaustion_20261101.py --output /tmp/toric-zero-exhaustion-review.json
```

[証明書](../../reports/TORIC_ZERO_JET_EXHAUSTION_20261101.json)は全行列、
一様な三ブロック、平行移動核の境界値、Cartan 軌道と Jordan 代表の一致を保存する。
既存の全境界 incidence と旧被約零モジュラー分類は、明示した依存定理である。

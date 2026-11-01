# 零一次部分の R0 型は AB に含まれる

2026年11月1日。複素数体上の被約な集合についての主張。

## 結論と依存関係

実際のトーリック成分の零一次部分を
\(F_0=(T01\cap\{\beta=0\})_{\mathrm{red}}\)、旧分類の閉 AB 成分を
\(\overline{AB}\) と書く。このノートは

\[
 F_0\cap R0\ \subset\ \overline{AB}
\]

を証明する。既に \(\overline{AB}\subset F_0\) は
[対数型退化の証明](TORIC_ZERO_JET_FIBER_20261101.md)で示した。
従って R0 内の集合としては \(F_0\cap R0=\overline{AB}\cap R0\) である。
全零一次部分の尽くし、交叉のスキーム構造は主張しない。

使う必要条件は、T01 の各点の**テンソル自体を固定する**射影 Lie 代数が
四次元の可換部分代数を含むことである。
[strict stabilizer のノート](TORIC_STRICT_STABILIZER_20261101.md)の
Grassmannian 上の閉 incidence による証明を使う。
テンソルの直線だけを固定する安定化群では不十分である。
また極限部分代数の二次外積像にテンソルが入る、という主張は使わない。

## 1. 固定定数因子における strict action

\(J=\partial_x\wedge\partial_y+\partial_z\wedge\partial_w\) とし、R0 を

\[
 P=J+Q,
 \qquad Q=\partial_y\wedge(-G_w\partial_z+G_z\partial_w),
\]

\[
G=h_0x^2z+h_1x^2w+h_2xz^2+h_3xzw+h_4xw^2
 +h_5z^3+h_6z^2w+h_7zw^2+h_8w^3
\]

と表す。\(G\) は \(x^3\) を法としている。
退化因子は \(x_0^5\) なので、射影自己同型の Lie 代数はアフィン座標を保つ。
アフィンベクトル場 \(a+A(x,y,z,w)^t\) の strict 条件は、次数ごとに

\[
 A\in\mathfrak{sp}_4,\qquad \partial_aQ=0,\qquad A\cdot Q=0.
\]

もし平行移動核が \(\mathbb C\partial_y\) だけなら、定数因子線も唯一である。
実際、二つの独立な定数因子を持つ場合は \(Q=F u\wedge v\) であり、
\([J,Q]=0\) から \(J^\sharp dF\in\langle u,v\rangle\)。
従って \(F\) は高々二つの線形座標に依存し、平行移動核は少なくとも二次元になる。

因子線 \(\mathbb C\partial_y\) を保つ \(\mathrm{Sp}_4\) の放物型 Lie 代数は
七次元で、そのうち \(y\mapsto y+cx\) は全 \(G\) に自明に作用する。
残りの作用は証明書の \(9\times6\) 行列 \(M(h)\) である。
従って平行移動核が一次元なら

\[
 \dim\mathfrak{aut}_{\mathrm{strict}}(P)=1+7-\operatorname{rank}M.
\]

必要条件は \(\operatorname{rank}M\le4\) を与える。
二元三次式 \(H=G(0,z,w)\) の根型を分ける。

- 三つの相異なる根：\(H=z^3+w^3\) とでき、下位五係数に依存しない
  六次小行列式が \(486\)。従って rank は少なくとも六。
- 二重根：\(H=z^2w\) とでき、五次小行列式が定数 \(4\)。
- 三重根：\(H=z^3\) とする。定数 \(9\) の三次 pivot を消去すると、
  残余行列の rank は高々一でなければならない。その二次小行列式から順に
  \(16h_4^2=0\)、\(-2h_3^3=0\)、\(10h_1^2=0\) を得る。
  従って \(G=z^3+h_2xz^2+h_0x^2z\) で、\(\partial_w\) も平行移動核に入る。
- \(H=0\)：\(G=xH_2(z,w)+x^2H_1(z,w)\)。\(H_2\) が rank 二なら
  \(H_2=zw\) とでき、定数 \(2\) の五次小行列式がある。
  rank 高々一なら \(H_2=L^2\) または零である。
  \(L\) を殺す \(z,w\) 方向の微分は \(G\) を \(\mathbb Cx^2\) に移すので、
  やはり追加の平行移動がある。

根型の正規化には因子線を保つ \(\mathrm{Sp}_4\) の作用を使える。
\(z,w\) の \(\mathrm{SL}_2\) と \(x,y\) の逆数スケールが、二元式の
通常の \(\mathrm{GL}_2\) 根型分類に必要なスカラーも与える。
従って strict 安定化次元が四以上なら、追加の定数平行移動が必ず存在する。

## 2. 追加の平行移動が x 成分を持たない場合

\(\partial_y\) を引き、\(z,w\) のシンプレクティック変換を施すと、追加の
平行移動を \(\partial_w\) とできる。このとき

\[
G=kx^2w+H(x,z)\pmod{x^3},\qquad
Q=-kx^2\partial_y\wedge\partial_z+H_z(x,z)\partial_y\wedge\partial_w.
\]

この族は係数を制限せず \(\overline{AB}\) に入る。次の明示的退化を使う。
AB の既知の正規形は

\[
P_{AB}=J-R^2\partial_R\wedge\partial_U-Z^2\partial_Z\wedge\partial_V
       +g(R,Z)\partial_U\wedge\partial_V,
\]

で、\(g\) は任意の二次式である。\(t d\ne0\) に対して

\[
\binom RZ=A_t\binom xz,\quad
\binom UV=A_t^{-t}\binom yw,\qquad
A_t=\begin{pmatrix}t&0\\2t&dt^3\end{pmatrix}
\]

という余接持ち上げを使う。これは \(J\) を固定する。
混合ブロックは

\[
A_t^{-1}\operatorname{diag}(-R^2,-Z^2)A_t
=\begin{pmatrix}
-t^2x^2&0\\
-6x^2/d-8t^2xz-2dt^4z^2&-(2tx+dt^3z)^2
\end{pmatrix}.
\]

任意の二次式 \(F(x,z)\) に対し
\(g_t(R,Z)=F(A_t^{-1}(R,Z)^t)/\det A_t\) と取れば、垂直ブロックは
すべての \(t\ne0\) でちょうど \(F(x,z)\partial_y\wedge\partial_w\) になる。
従って極限は

\[
J+\frac6d x^2\partial_y\wedge\partial_z
 +F(x,z)\partial_y\wedge\partial_w.
\]

\(d\) により任意の非零係数が得られ、零係数も閉包に入る。
\(g_t\) の係数が極を持つことは問題でない。引き戻したテンソルの全係数は
上式の通り \(t=0\) で正則であり、\(t\ne0\) の各点は AB に属する。

## 3. 追加の平行移動が x 成分を持つ場合

その成分を一にし、\(\partial_y\) を引いて
\(d=\partial_x+a\partial_z+b\partial_w\) とする。
シンプレクティック座標

\[
X=x,\quad Z=z-ax,\quad W=w-bx,\quad Y=y-bz+aw
\]

で \(d=\partial_X\)、定数因子は \(\partial_Y\) のままである。
従って \(G=H(Z,W)\pmod{X^3}\) とできる。二元三次式の四根型で十分である。

相異なる三根では定数因子は唯一、平行移動核は二次元、線形作用 rank は六なので、
strict 安定化次元は \(2+7-6=3\) で除外される。
三重根は \(G=z^3\) で、前節の \(k=0\) の族に入る。零式も同様である。

二重根 \(G=z^2w\) では四次元という次元条件だけでは除外できない。
しかし全射影 strict 安定化 Lie 代数の基底は、斉次座標順序
\((x_0,x,y,z,w)\) で

\[
E_{10},\quad E_{20},\quad E_{21},\quad
\tfrac12\operatorname{diag}(0,-1,1,1,-1).
\]

証明書は全 \(25\) 行列方向の作用の rank が \(20\) であることと、
この四方向に Euler 方向を加えた五方向が核を尽くすことを有理数で確認する。
\([E_{21},E_{10}]=E_{20}\ne0\) なのでこの四次元 Lie 代数は非可換である。
四次元可換部分代数を含めず、T01 に属せない。
これは二重根の全軌道の排除であり、一点の一般性からの推測ではない。

以上で R0 の全場合を覆った。

## 証明書と再現

```
.venv/bin/python code/toric_zero_jet_r0_20261101.py --output /tmp/toric-r0-review.json
```

[証明書](../../reports/TORIC_ZERO_JET_R0_20261101.json)は strict 作用行列を
\(\mathfrak{sp}_4\) の作用から再構成し、定数小行列式、三重根の Schur 消去、
全射影安定化核と非零括弧、AB からのシンプレクティック退化を確認する。
旧 Q45 の全二十二成分分類の再検証はこの R0 内の主張には不要である。
他の modular 型が AB で尽くされるかは別途残っている。

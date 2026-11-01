# M0・M1 の実際のトーリック零一次部分は AB に入る

2026年11月1日。複素数体上の被約集合についての部分分類。

## 結論

\(F_0=(T01\cap\{\beta=0\})_{\mathrm{red}}\) とし、旧最小冪零トレースの
二族を M0、M1 と書く。このノートは

\[
 F_0\cap(M0\cup M1)\subset\overline{AB}
\]

を証明する。R0 内の同様の結果は [別ノート](TORIC_ZERO_JET_R0_20261101.md)にある。
必要条件は、strict Poisson 安定化 Lie 代数が四次元可換部分代数を含むこと。
十分性には [AB の対数型退化](TORIC_ZERO_JET_FIBER_20261101.md)を使う。
第三の最小トレース成分、その他の modular 型の全特殊点はここでは扱わない。

## 1. 正規形と作用行列

\(J=\partial_x\wedge\partial_y+\partial_z\wedge\partial_w\) とし、
\(G\) と \(h_0,\ldots,h_8\) は R0 ノートと同じ三次式である。

\[
\begin{aligned}
P_0&=(1+x^2)\partial_x\wedge\partial_y+\partial_z\wedge\partial_w
 +\partial_y\wedge(-G_w\partial_z+G_z\partial_w),\\
P_1&=(1+x^2)(\partial_x\wedge\partial_y+\partial_z\wedge\partial_w)
 +\partial_y\wedge((-G_w-xz)\partial_z+(G_z-xw)\partial_w).
\end{aligned}
\]

退化因子はそれぞれ \(x_0^3(x_0^2+x^2)\)、\(x_0(x_0^2+x^2)^2\)。
従って連結自己同型の Lie 代数は三つの各超平面を保つ。
Euler 方向を除くと、斉次行列の \(x_0,x\) 行は零である。
定数 Poisson テンソルを固定する条件を解くと九次元。
この中には二つの普遍的な \(y\) shear、二つの \(z,w\) 平行移動、
五つの有効な斉次作用方向がある。

平行移動条件の rank を \(b\in\{0,1,2\}\)、五方向の二次作用行列の
rank を \(r\) とすると strict 安定化次元は \(9-b-r\) である。
一次項と二次項は独立なので、この次元式には異次数間の相殺がない。
証明書は全作用から \(12\times5\) 行列 \(M_0,M_1\) を再構成する。

まず追加の平行移動がない場合 \(b=2\) を除外する。この場合、必要条件は
\(r\le3\) である。

## 2. 非零二元三次式の排除

\(H=G(0,z,w)\) の根型を分ける。
三つの相異なる根なら \(H=\lambda(z^3+w^3)\)、\(\lambda\ne0\) とでき、
下位五係数によらない四次小行列式は \(972\lambda^4\)。
二重根 \(H=z^2w\) では定数 \(24\) の四次小行列式がある。
いずれも \(r\le3\) に反する。

三重根 \(H=z^3\) では行 \((3,7,9)\)、列 \((1,0,4)\) の pivot が
定数 \(-162\)。その Schur 残差が零であるためには

\[
-h_3-\epsilon=0,\qquad -h_3(h_3-\epsilon)/3=0,
\quad \epsilon=1\ (M0),\ -1\ (M1)
\]

が必要だが、第一式を第二式に代入すると \(-2/3=0\) となる。
従って \(b=2\) で可能なのは \(H=0\) だけである。

## 3. H=0 の全場合

\[
G=x(h_2z^2+h_3zw+h_4w^2)+x^2(h_0z+h_1w)
 =xH_2+x^2H_1.
\]

シンプレクティック shear

\[
(X,Y,Z,W)=(x,y-bz+aw,z-ax,w-bx)
\]

は \(H_2\) を固定し、係数を正確に

\[
 h'_0=h_0+2ah_2+b(h_3+\epsilon),\qquad
 h'_1=h_1+a(h_3-\epsilon)+2bh_4
\]

と変える。以下、\(z,w\) の \(\mathrm{SL}_2\) 作用で \(H_2\) を分類する。

### M0

\(H_2\) が rank 高々一なら、\(H_2=L^2\) または零。
\(L\) を殺す平行移動は \(Q\) を固定するため \(b<2\) である。

rank 二なら \(H_2=\kappa zw\)、\(\kappa\ne0\) とできる。
\(\kappa\ne\pm1\) では shear により \(h_0=h_1=0\) とでき、
四次小行列式 \(4\kappa^2(\kappa-1)(\kappa+1)\ne0\) がある。
\(z,w\) の交換により共鳴値は \(\kappa=1\) だけでよい。
shear で \(h_0=0\) とし、\(h_1=k\) とすると四次小行列式は \(-8k\)。
従って \(k\ne0\) も除外される。

残る \(G=xzw\) の全 strict 安定化は正確に四次元で、基底は

\[
 E_{20},\ E_{21},\ E_{31}-E_{24},\ \tfrac12\operatorname{diag}(0,0,0,1,-1).
\]

最後の要素と第三の要素の括弧は第三の要素の \(1/2\) 倍である。
従って非可換で、四次元可換部分代数を含まない。
以上から M0 に属する T01 点は必ず \(b<2\) である。

### M1

\(H_2\) が非零 rank 一なら \(H_2=z^2\) とでき、
任意の \(h_0,h_1\) に対して四次小行列式が定数 \(8\) である。
rank 二で \(\kappa\ne\pm1\) なら、M0 と同じ shear と非零 minor で除外される。
\(\kappa=1\) では shear により
\(G=xzw+kx^2z\) とできるが、このテンソルは \(w\) に依存せず、既に \(b<2\)。
\(\kappa=-1\) は交換により同じ場合となる。

\(H_2=0\) では shear により \(H_1=0\)。\(G=0\) の全 strict 安定化は

\[
 \langle E_{20},E_{21}\rangle\ \oplus\
 \langle \tfrac12\operatorname{diag}(0,0,0,1,-1),E_{34},E_{43}\rangle
 \simeq\mathbb C^2\oplus\mathfrak{sl}_2.
\]

\(\mathfrak{sl}_2\) の可換部分代数は高々一次元なので、全体の可換部分代数も
高々三次元。従ってこの場合も T01 ではない。
M1 でも T01 点は必ず \(b<2\) である。

ここでの二つの全安定化代数は、全 \(25\) 射影行列方向の rank と核を
有理数で確認している。単に表示した部分代数だけを調べたのではない。

## 4. b<2 の族はすべて AB の極限

追加の \(z,w\) 平行移動を \(\partial_z\) に取れる。この場合の全ポテンシャルは

\[
 M0:\quad G=\alpha x^2z+H(x,w),\qquad
 M1:\quad G=-xzw+\alpha x^2z+H(x,w),
\]

である。\(H\) は \(x^3\) を法として任意の三次式。
位置座標を \((x,w)\)、運動量を \((y,-z)\) とする。
AB の既知の正規形に、次の位置変換とその余接持ち上げを施す。

\[
 A_t^{(0)}=\begin{pmatrix}i&0\\\alpha t&t\end{pmatrix},\qquad
 A_t^{(1)}=\begin{pmatrix}i&0\\i-\alpha t/2&t\end{pmatrix}.
\]

\(\binom RZ=A_t\binom xw\) とおくと、二次混合ブロックの極限はそれぞれ

\[
\lim_{t\to0}A_t^{-1}\operatorname{diag}(-R^2,-Z^2)A_t
=\begin{pmatrix}x^2&0\\-\alpha x^2&0\end{pmatrix},\qquad
\begin{pmatrix}x^2&0\\2xw-\alpha x^2&x^2\end{pmatrix}.
\]

AB の任意二次式の係数を
\(g_t(R,Z)=H_w(A_t^{-1}(R,Z)^t)/\det A_t\) と取れば、垂直ブロックは
\(H_w\partial_y\wedge\partial_{-z}=-H_w\partial_y\wedge\partial_z\) となる。
従って極限は正確に上の M0、M1 の全 \(b<2\) 族である。
引き戻したテンソルの係数は \(t=0\) で正則であり、\(t\ne0\) では AB の点。
閉包を取って包含を得る。

## 証明書と範囲

```
.venv/bin/python code/toric_zero_jet_minimal_20261101.py --output /tmp/toric-minimal-review.json
```

[証明書](../../reports/TORIC_ZERO_JET_MINIMAL_20261101.json)は全作用からの行列の再構成、
定数および一変数 minor、Schur 矛盾、全係数 shear 恒等式、例外点の全 Lie 核と括弧、
任意係数を残した二つの AB 退化を含む。
旧分類の第三最小トレース成分との特殊交叉や、nil22・零モジュラー型の残りは別問題。
全 S61 の一般ファイバーを尽くしたという主張は含まない。

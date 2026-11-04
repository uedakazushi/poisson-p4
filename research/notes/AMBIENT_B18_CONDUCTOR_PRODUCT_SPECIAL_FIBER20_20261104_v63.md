# B18 conductor 積の実 DVR 閉ファイバー：長さ20と乗法表

2026-11-04、v63。元の深さ・socle 問題のための有限入力を一つ確定する。

実共通正則底を (R_0)、両実 conductor ideal を (I_A,I_B) とし、

\[
 H=R_0/(I_A I_B),\qquad \ell=5(s_A-r_A)
\]

と置く。ここで ((z_0,\ldots,z_6,r_A,s_A)) は凍結済みの実 CA F9 座標であり、
型だけを合わせた抽象座標ではない。既受理の conductor DVR 比較定理は
\(\operatorname{length}_{\mathbb Q}H/(\ell)=20\) を与える。
今回、この実閉ファイバーの全八座標乗法表を有理数上で構成した。

結果は

\[
 \mathfrak m^5(H/(\ell))=0,\qquad
 \operatorname{HF}(H/(\ell))=(1,8,8,2,1).
\]

一般 DVR ファイバーの作用・比較 kernel の torsion・元の (D_0) の深さは計算していない。

## 実 source から必要な係数を得る

入力は既存の二つの57変数 source の全二次係数、CA 零平面上の105変数 source
の三次係数、CB 自身の零平面上の三次係数である。新たな source 方程式や高次
Newton 段階は解いていない。入力の全ファイル hash、実 physical rows9、F9 と
conductor の元 node ID を証拠へ保存した。

両二次 physical 座標写像を (x_A(Y_A),x_B(Y_B)) とする。
定数線形逆行列だけを用い、(x_B(Y_B(Y_A))=x_A(Y_A)\pmod{\mathfrak m^3})
から実座標変換の全二次係数を復元し、この合成恒等式も係数ごとに確認した。
その一次部分は既存の凍結 adapter と一致する。最初の七座標と末尾の二座標の間の
一次 off-diagonal block は両方零である。

全次数で成り立つ ordered-prefix Hadamard 恒等式を

\[
 z_B=A(z,r,s)z+b(r,s),\quad A=A_0+A_1+O(\mathfrak m^2)
\]

と定める。(A_0) の行列式は \(-1/1024\) であり単位である。
二次 monomial が (z) を含む場合、最大の添字の (z_j) で割る規約を用いた。
これは prefix による正則商の規約そのものである。

IB の最初の七生成元を (A^{-1}z_B) へ替えることは、実 ideal の単位行列変換である。
その三次までの必要係数は

\[
 z'_B=z+\beta_2+\beta_3+O(\mathfrak m^4),\qquad
 \beta_2=A_0^{-1}b_2,\quad
 \beta_3=A_0^{-1}(b_3-A_1A_0^{-1}b_2).
\]

従って全九変数三次 source は不要である。(A_1) は全二次 cache、
 (b_2,b_3) は実 CA 零平面上の既存三次 cache が与える。
未知の (A_2) は (b_2) との積で四次以降にしか寄与しない。
IA の線形生成元との積ではその誤差は五次へ入る。

両 ideal の最後の生成元も、それぞれ自身の最初の七零座標での Hadamard 差を引いて
実平面関数へ替える。これは actual ideal の合法な生成元変換である。
CB の末尾二座標の一次部分に (z) がないことにより、その cubic も既存 cache だけで足りる。
\(\ell=0\)、すなわち正確に (s=r\) とすると、最後の生成元の必要な初期式は

\[
 g_{A,3}=\frac{125}{24}r^3,\qquad
 g_{B,3}=-\frac{15125}{144}r^3.
\]

両方の actual source 祖先と元 cubic の係数は保存している。
IB の最後の生成元を型だけから (r^3) と仮定していない。

参考に、正規化後の二次差は

\[
 \beta_2=r^2\left(\frac{9835}{864},0,0,0,
 \frac{29875}{3456},-\frac{4565}{864},0\right).
\]

全三次差も証拠に明示した。

## 有限商の係数検算

八変数 (z_0,\ldots,z_6,r) の四次以下の495 monomial を用いた。
上の actual IA8 と IB8 の**全64組の積**を五次未満へ切り、その全ての必要な
monomial 倍を含めた。定数倍だけ異なる行を保存済み係数で同定した結果は
579行・824非零係数であり、有理数上の rank は475であった。
各 pivot について元の積・乗数・過去の pivot を引いた係数・正規化係数を保存した。

得られた20個の monomial 基底は

\[
 1,\quad z_0,\ldots,z_6,r,\quad
 rz_0,\ldots,rz_6,r^2,\quad r^2z_5,r^3,\quad r^4.
\]

全495 monomial の基底表示、八座標の20×20乗法行列を列方向 sparse 形式で保存した。
八行列の可換性、全64 actual 積が全20基底元を殺すこと、全 relation 行が零へ落ちることを
正確に確認した。冪次元は (20,19,11,3,1,0\) である。

この係数計算はまず

\[
 H/(\ell)\twoheadrightarrow H/(\ell,\mathfrak m^5)
\]

の右辺の長さ20を示す。左辺も既受理の actual DVR 定理により長さ20なので、
この全射は同型である。したがって三次 cache から作った乗法表は actual 閉ファイバー
そのものの表であり、五次消滅も実環についての結論となる。

## 再現と証拠

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_conductor_product_special_fiber20_20261104_v63.py
```

計算本体0.162秒、最大 RSS 約28 MB。新しい source solve、有限体、点の選択、GB は使わない。
主要証拠は
`reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json`。
この閉ファイバーの正確な入力は、DVR 上の実比較 module・完全 syzygy・kernel torsion を
今後構成するためのものであり、その未計算の性質を特殊化 kernel から推論しない。

# 実 Cbar の三つの cyclic sector と全最小自由分解

2026-11-04、v63。前 packet の実閉 cokernel
\(\bar C=C/(\ell)\)、\(\ell=5(s_A-r_A)\) の八つの作用行列だけを用いる。
追加 source solve、jet、rank search はない。

\(\bar R=R_0/(\ell)=\mathbb Q[[z_0,\ldots,z_6,r]]\) は実 CA F9 座標である。
凍結済み八基底は、actual normalizer tensor の三 target により正確に

\[
 \tau_A,r\tau_A,r^2\tau_A,
 \tau_B,r\tau_B,r^2\tau_B,
 \tau_A\tau_B,r\tau_A\tau_B
\]

となる。この持上げは source の次数二以下の monomial であり、抽象的な同型を選んだものではない。
raw26 基底では位置 (0,1,2,3,4,5,6,14\)、三 target は位置 (0,3,6\) である。

保存済み実作用行列には混合 sector 成分がない。(r\) は上の三つの chain に作用し、長さは3,3,2。
最初と最後の sector で全 (z_i\) は零、中央 sector では (-c_i r^2\) として作用する。ここで

\[
 c=\left(\frac{9835}{864},0,0,0,
 \frac{29875}{3456},-\frac{4565}{864},0\right)
\]

は actual product20 packet の \(\beta_2/r^2\) そのものである。
全八行列の全成分でこの表示を確認した。従って

\[
 \bar C\cong
 \frac{\bar R}{(z_0,\ldots,z_6,r^3)}\ \oplus\
 \frac{\bar R}{(z_0+c_0r^2,\ldots,z_6+c_6r^2,r^3)}\ \oplus\
 \frac{\bar R}{(z_0,\ldots,z_6,r^2)}.
\]

各商の長さは3,3,2であり、その和が actual 八基底を全単射に送るので、これは全次数の実 module 同型である。
単に Betti 数や Loewy 数の一致から型を推測していない。

三つの八生成 regular sequence の Koszul complex の直和を、全て具体的に保存した。
各自由基底の sector と wedge の添字、全八 differential の全3072非零多項式成分がある。
全 chain 積の零、augmentation の実八作用表との一致を検査した。
中央 sector の座標変換 (z_i\mapsto z_i+c_i r^2\) は三角 automorphism なので exactness は各 regular sequence から従う。
全係数が maximal ideal に入るため最小であり、rank 列は

\[
 3,24,84,168,210,168,84,24,3
\]

である。これは先行の独立 Koszul homology 計算の Betti 列とも一致する。

同じ構成を (S=\mathbb Q[[r]]\) 上から記述すると、

\[
 D_0=\operatorname{diag}(r^3,r^3,r^2),\qquad
 T_i=S_i=\operatorname{diag}(0,-c_i r^2,0),\qquad H_{ij}=0.
\]

全 (T_iD_0=D_0S_i\)、可換子、coherence/Jacobi 恒等式を検査した。
ここでの (S\) は閉 module の nilpotent (r\) 作用を表示するための環であり、
元の比較で用いる (\mathbb Q[[\ell]]\) 上の一般 rank に関する主張ではない。

主要証拠は
`reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/actual-Cbar-minimal-resolution-and-three-targets.json`。
一回の固定表検査は0.076秒、約16 MB、証拠126 KBで完了した。
比較 (T\)、\(K\) の長さ・torsion・深さは計算していない。

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_Cbar_literal_resolution_20261104_v63.py
```
